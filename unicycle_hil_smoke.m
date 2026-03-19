function results = unicycle_hil_smoke(ctrl, cfg)
%UNICYCLE_HIL_SMOKE  Synthetic HIL smoke test against the MATLAB plant.
%  Mirrors the intent of hil_plug_play_smoke.py from the sim-to-real repo.
%
%  results = unicycle_hil_smoke(ctrl, cfg)
%    ctrl    - controller struct from unicycle_design_controllers
%    cfg     - configuration from unicycle_config()
%    results - struct with:
%                .passed
%                .scenarios
%                .report_file

if nargin < 2 || isempty(cfg)
    cfg = unicycle_config();
end
if nargin < 1 || isempty(ctrl)
    [A0, B0, C0, ~] = unicycle_plant(cfg);
    ctrl = unicycle_design_controllers(A0, B0, C0, cfg);
end

if isfield(ctrl, 'lqr_current')
    ctrl_lqr = ctrl.lqr_current;
else
    ctrl_lqr = ctrl;
end

[A, B, C, ~] = unicycle_plant(cfg);
[Ad, Bd] = unicycle_discretize(A, B, cfg.dt);

scenario_defs = [ ...
    struct('name', 'Idle', ...
           'description', 'Idle (zero torque at rest)', ...
           'x0', zeros(5,1), ...
           'mode', 'normal'); ...
    struct('name', 'Pitch forward', ...
           'description', 'Pitch forward -> +tau_base', ...
           'x0', [deg2rad(5); 0; 0; 0; 0], ...
           'mode', 'normal'); ...
    struct('name', 'Pitch backward', ...
           'description', 'Pitch backward -> -tau_base', ...
           'x0', [deg2rad(-5); 0; 0; 0; 0], ...
           'mode', 'normal'); ...
    struct('name', 'Roll right', ...
           'description', 'Roll right -> +tau_rw', ...
           'x0', [0; 0; deg2rad(5); 0; 0], ...
           'mode', 'normal'); ...
    struct('name', 'Roll left', ...
           'description', 'Roll left -> -tau_rw', ...
           'x0', [0; 0; deg2rad(-5); 0; 0], ...
           'mode', 'normal'); ...
    struct('name', 'Tilt ESTOP', ...
           'description', 'Tilt ESTOP (>35 deg)', ...
           'x0', [deg2rad(40); 0; 0; 0; 0], ...
           'mode', 'normal'); ...
    struct('name', 'Timeout', ...
           'description', 'Timeout (frozen sensor)', ...
           'x0', zeros(5,1), ...
           'mode', 'freeze')];

scenario_results = repmat(struct( ...
    'name', '', ...
    'description', '', ...
    'passed', false, ...
    'detail', '', ...
    'fail_reason', '', ...
    'survival', false, ...
    'final_state', zeros(5,1), ...
    'max_u', 0, ...
    'res', struct()), numel(scenario_defs), 1);

report_lines = {};
report_lines{end+1} = sprintf('\n=== HIL SMOKE TEST RESULTS ==='); %#ok<AGROW>
report_lines{end+1} = sprintf('%-6s  %-30s  %-8s  %s', 'Scen', 'Description', 'Result', 'Detail'); %#ok<AGROW>
report_lines{end+1} = repmat('-', 1, 90); %#ok<AGROW>

for si = 1:numel(scenario_defs)
    sc = scenario_defs(si);
    res = run_smoke_episode(Ad, Bd, C, ctrl_lqr, cfg, sc.x0, sc.mode);

    scenario_results(si).name = sc.name;
    scenario_results(si).description = sc.description;
    scenario_results(si).survival = res.survival;
    scenario_results(si).final_state = res.x_true(:, end);
    scenario_results(si).max_u = max(abs(res.u(:)));
    scenario_results(si).res = res;

    [passed, detail, fail_reason] = evaluate_scenario(si, res, cfg);
    scenario_results(si).passed = passed;
    scenario_results(si).detail = detail;
    scenario_results(si).fail_reason = fail_reason;

    report_lines{end+1} = sprintf('%-6s  %-30s  %-8s  %s', ... %#ok<AGROW>
        sprintf('S%d', si), sc.description, pass_fail(passed), detail);
end

sum_passed = sum([scenario_results.passed]);
report_lines{end+1} = repmat('-', 1, 90); %#ok<AGROW>
report_lines{end+1} = sprintf('PASSED: %d / 7', sum_passed); %#ok<AGROW>
if sum_passed == numel(scenario_results)
    report_lines{end+1} = 'ALL SMOKE TESTS PASS -- sim stack sign conventions verified'; %#ok<AGROW>
else
    report_lines{end+1} = 'SMOKE TEST FAILURES -- check sign conventions before hardware'; %#ok<AGROW>
end

fprintf('%s\n', strjoin(report_lines, newline));

if ~exist(cfg.results_dir, 'dir')
    mkdir(cfg.results_dir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
report_file = fullfile(cfg.results_dir, sprintf('hil_smoke_%s.txt', timestamp));
fid = fopen(report_file, 'w');
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', strjoin(report_lines, newline));

results.passed = all([scenario_results.passed]);
results.scenarios = scenario_results;
results.report_file = report_file;
end


function res = run_smoke_episode(Ad, Bd, C, ctrl, cfg, x0, mode)
dt = cfg.dt;
N = round(cfg.hil_smoke_duration / dt);
nx = size(Ad, 1);
ny = size(C, 1);
nu = size(Bd, 2);

x_true = zeros(nx, N + 1);
x_est = zeros(nx, N + 1);
u_hist = zeros(nu, N);
y_hist = zeros(ny, N);

x_true(:,1) = x0;
x_est(:,1) = x0;

freeze_step = max(round(0.5 / dt), 1);
y_hold = [];
freeze_count = 0;
timeout_zero_step = N + 1;
crashed = false;
crash_step = N + 1;
crash_t = cfg.hil_smoke_duration;

if is_crashed_state(x_true(:,1), cfg)
    crashed = true;
    crash_step = 1;
    crash_t = 0.0;
    x_true(:, 2:end) = repmat(x_true(:,1), 1, N);
    x_est(:, 2:end) = repmat(x_est(:,1), 1, N);
    res = pack_episode_result(x_true, x_est, u_hist, y_hist, crashed, crash_step, crash_t);
    return;
end

for k = 1:N
    y_now = C * x_true(:,k);
    if strcmp(mode, 'freeze') && k >= freeze_step
        if isempty(y_hold)
            y_hold = y_now;
        end
        y_used = y_hold;
        freeze_count = freeze_count + 1;
    else
        y_used = y_now;
        if strcmp(mode, 'freeze')
            y_hold = y_used;
        end
    end
    y_hist(:,k) = y_used;

    x_est(:,k) = x_est(:,k) + ctrl.L * (y_used - C * x_est(:,k));
    if strcmp(mode, 'freeze') && freeze_count >= cfg.hil_timeout_steps
        u_k = zeros(nu, 1);
        timeout_zero_step = k;
    else
        u_k = -ctrl.K * x_est(:,k);
        u_k = saturate_control(u_k, cfg);
    end
    u_hist(:,k) = u_k;

    x_true(:,k+1) = x_true(:,k) + dt * local_dynamics(x_true(:,k), u_k, cfg);
    x_est(:,k+1) = Ad * x_est(:,k) + Bd * u_k;

    if is_crashed_state(x_true(:,k+1), cfg)
        crashed = true;
        crash_step = k + 1;
        crash_t = k * dt;
        if k < N
            x_true(:,k+2:end) = repmat(x_true(:,k+1), 1, N - k);
            x_est(:,k+2:end) = repmat(x_est(:,k+1), 1, N - k);
            u_hist(:,k+1:end) = 0;
            y_hist(:,k+1:end) = repmat(y_used, 1, N - k);
        end
        break;
    end
end

res = pack_episode_result(x_true, x_est, u_hist, y_hist, crashed, crash_step, crash_t);
res.freeze_step = freeze_step;
res.timeout_zero_step = timeout_zero_step;
end


function res = pack_episode_result(x_true, x_est, u_hist, y_hist, crashed, crash_step, crash_t)
res.x_true = x_true;
res.x_est = x_est;
res.u = u_hist;
res.y = y_hist;
res.crashed = crashed;
res.crash_step = crash_step;
res.crash_t = crash_t;
res.survival = ~crashed;
end


function [passed, detail, fail_reason] = evaluate_scenario(idx, res, cfg)
passed = false;
detail = '';
fail_reason = '';

switch idx
    case 1
        max_u = max(abs(res.u(:)));
        passed = max_u < 0.01;
        detail = sprintf('max|u|=%.3f N*m', max_u);
        if ~passed
            fail_reason = sprintf('equilibrium torque too large (max|u|=%.4f N*m)', max_u);
        end

    case 2
        u_base = res.u(2,1);
        passed = u_base > 0;
        detail = sprintf('u_base(1)=%+.3f N*m', u_base);
        if ~passed
            fail_reason = sprintf('expected positive base torque, got %.4f N*m', u_base);
        end

    case 3
        u_base = res.u(2,1);
        passed = u_base < 0;
        detail = sprintf('u_base(1)=%+.3f N*m', u_base);
        if ~passed
            fail_reason = sprintf('expected negative base torque, got %.4f N*m', u_base);
        end

    case 4
        u_rw = res.u(1,1);
        passed = u_rw > 0;
        detail = sprintf('u_rw(1)=%+.3f N*m', u_rw);
        if ~passed
            fail_reason = sprintf('expected positive reaction-wheel torque, got %.4f N*m', u_rw);
        end

    case 5
        u_rw = res.u(1,1);
        passed = u_rw < 0;
        detail = sprintf('u_rw(1)=%+.3f N*m', u_rw);
        if ~passed
            fail_reason = sprintf('expected negative reaction-wheel torque, got %.4f N*m', u_rw);
        end

    case 6
        post_start = min(max(res.crash_step, 1), size(res.u, 2) + 1);
        if post_start <= size(res.u, 2)
            post_crash_max = max(abs(res.u(:, post_start:end)), [], 'all');
        else
            post_crash_max = 0.0;
        end
        passed = res.crashed && res.crash_t < 0.1 && post_crash_max < 1e-6;
        detail = sprintf('crash_t=%.3f s, post|max u|=%.3e', res.crash_t, post_crash_max);
        if ~passed
            fail_reason = sprintf('crash=%d crash_t=%.4f post|maxu|=%.3e', ...
                res.crashed, res.crash_t, post_crash_max);
        end

    case 7
        freeze_step = res.freeze_step;
        max_u_freeze = max(abs(res.u(:, freeze_step:end)), [], 'all');
        finite_u = all(isfinite(res.u(:)));
        passed = finite_u && ~res.crashed && max_u_freeze < cfg.tau_rw_max;
        detail = sprintf('max|u|=%.3f N*m', max_u_freeze);
        if ~passed
            fail_reason = sprintf('finite=%d crashed=%d max|u|=%.4f N*m after freeze', ...
                finite_u, res.crashed, max_u_freeze);
        end
end
end


function u = saturate_control(u, cfg)
u(1) = max(min(u(1), cfg.tau_rw_max), -cfg.tau_rw_max);
u(2) = max(min(u(2), cfg.tau_base_max), -cfg.tau_base_max);
end


function tf = is_crashed_state(x, cfg)
tf = abs(x(1)) > cfg.pitch_crash_rad || abs(x(3)) > cfg.roll_crash_rad;
end


function xdot = local_dynamics(x, u, cfg)
m = cfg.m_body + cfg.payload_mass;
L = cfg.L_body;
Ib_p = cfg.I_body + m * L^2;
Ib_r = cfg.I_body_r;

ddtheta_p = (m * cfg.g * L * sin(x(1)) - u(2)) / Ib_p;
ddtheta_r = (m * cfg.g * L * sin(x(3)) - u(1)) / Ib_r;
domega_rw = u(1) / cfg.I_rw;
xdot = [x(2); ddtheta_p; x(4); ddtheta_r; domega_rw];
end


function txt = pass_fail(tf)
if tf
    txt = 'PASS';
else
    txt = 'FAIL';
end
end
