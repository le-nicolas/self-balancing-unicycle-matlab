% VALIDATION DESIGN RULE:
% This system is an open-loop unstable inverted pendulum.
% Any test that integrates the plant for more than ~0.3 s
% from a non-zero IC MUST apply a stabilising controller.
% Open-loop integration is only valid for:
%   (a) one-step impulse-response checks  (Stage 5)
%   (b) confirming divergence  (Stage 6a)
%   (c) zero-IC equilibrium checks  (Stage 2a)
% All other tests run closed-loop.
function unicycle_model_validate()
%UNICYCLE_MODEL_VALIDATE  Validate the unicycle plant and simulator model.
%  Prints PASS/FAIL for each validation stage and throws an error if any
%  stage fails. This file is self-contained and does not modify model code.

clc;
fprintf('=== UNICYCLE MODEL VALIDATION ===\n\n');
prev_rng = rng;
cleanup_rng = onCleanup(@() rng(prev_rng)); %#ok<NASGU>
rng(1234);

cfg = unicycle_config();
[A, B, C, ~] = unicycle_plant(cfg);
[Ad, Bd] = unicycle_discretize(A, B, cfg.dt);
ctrl_all = unicycle_design_controllers(A, B, C, cfg);
L = ctrl_all.lqr_current.L;

stage_names = { ...
    'Stage 1  Symbolic linearisation'; ...
    'Stage 2  Equilibrium & signs'; ...
    'Stage 3  Lin vs nonlin agreement'; ...
    'Stage 4  ZOH discretisation'; ...
    'Stage 5  Step-response signs'; ...
    'Stage 6  OL diverge / CL converge'; ...
    'Stage 7  Axis decoupling'; ...
    'Stage 8  Kalman audit'};
stage_pass = false(8, 1);
failed_stages = {};

%% Stage 1 - Jacobian audit
fprintf('--- STAGE 1: SYMBOLIC LINEARISATION CHECK ------------------------\n');
[stage_pass(1), msgs] = stage1_linearisation(A, B, cfg);
print_stage_messages(msgs);
record_stage_result(1);

%% Stage 2 - Equilibrium and sign checks
fprintf('\n--- STAGE 2: EQUILIBRIUM AND ENERGY CHECK ------------------------\n');
[stage_pass(2), msgs] = stage2_signs(cfg);
print_stage_messages(msgs);
record_stage_result(2);

%% Stage 3 - Linear vs nonlinear agreement
fprintf('\n--- STAGE 3: LINEAR VS NONLINEAR AGREEMENT -----------------------\n');
[stage_pass(3), msgs] = stage3_lin_vs_nl(A, B, ctrl_all.lqr_current.K, cfg);
print_stage_messages(msgs);
record_stage_result(3);

%% Stage 4 - Discretisation audit
fprintf('\n--- STAGE 4: ZOH DISCRETISATION AUDIT ----------------------------\n');
[stage_pass(4), msgs] = stage4_discretisation(A, B, Ad, Bd, cfg);
print_stage_messages(msgs);
record_stage_result(4);

%% Stage 5 - Discrete impulse sign audit
fprintf('\n--- STAGE 5: DISCRETE STEP-RESPONSE SIGN CHECK -------------------\n');
[stage_pass(5), msgs] = stage5_step_signs(Ad, Bd);
print_stage_messages(msgs);
record_stage_result(5);

%% Stage 6 - Open-loop divergence and closed-loop convergence
fprintf('\n--- STAGE 6: OPEN-LOOP / CLOSED-LOOP TRAJECTORIES ----------------\n');
[stage_pass(6), msgs] = stage6_open_closed_loop(ctrl_all.lqr_current.K, cfg);
print_stage_messages(msgs);
record_stage_result(6);

%% Stage 7 - Axis decoupling
fprintf('\n--- STAGE 7: AXIS DECOUPLING CHECK -------------------------------\n');
[stage_pass(7), msgs] = stage7_axis_decoupling(cfg);
print_stage_messages(msgs);
record_stage_result(7);

%% Stage 8 - Kalman audit
fprintf('\n--- STAGE 8: KALMAN BIAS AND OBSERVABILITY AUDIT -----------------\n');
[stage_pass(8), msgs] = stage8_kalman_audit(Ad, Bd, C, ctrl_all.lqr_current.K, L, cfg);
print_stage_messages(msgs);
record_stage_result(8);

%% Summary
fprintf('\n--- SUMMARY ------------------------------------------------------\n');
for i = 1:numel(stage_names)
    fprintf('%-36s %s\n', stage_names{i}, pass_fail_text(stage_pass(i)));
end

n_pass = sum(stage_pass);
if n_pass == numel(stage_pass)
    fprintf('\nOVERALL: PASS (%d/8)\n', n_pass);
else
    failed_list = strjoin(failed_stages, ', ');
    fprintf('\nOVERALL: FAIL (%d/8 -- stages %s failed)\n', n_pass, failed_list);
    error('unicycle_model_validate:failed', ...
        'Model validation failed. Fix stages: %s', failed_list);
end

    function record_stage_result(stage_idx)
        fprintf('%s\n', stage_banner(stage_idx, stage_pass(stage_idx)));
        if ~stage_pass(stage_idx)
            failed_stages{end+1} = num2str(stage_idx); %#ok<AGROW>
        end
    end
end


function [pass, msgs] = stage1_linearisation(A, B, cfg)
eps_fd = 1e-6;
x0 = zeros(5, 1);
u0 = zeros(2, 1);
A_fd = zeros(size(A));
B_fd = zeros(size(B));

for i = 1:size(A, 2)
    xp = x0;
    xp(i) = xp(i) + eps_fd;
    A_fd(:, i) = (local_dynamics(xp, u0, cfg) - local_dynamics(x0, u0, cfg)) / eps_fd;
end

for j = 1:size(B, 2)
    up = u0;
    up(j) = up(j) + eps_fd;
    B_fd(:, j) = (local_dynamics(x0, up, cfg) - local_dynamics(x0, u0, cfg)) / eps_fd;
end

a_err = A - A_fd;
b_err = B - B_fd;
a_norm = norm(a_err, 'fro');
b_norm = norm(b_err, 'fro');

msgs = {
    sprintf('A_fd Frobenius error = %.3e', a_norm)
    sprintf('B_fd Frobenius error = %.3e', b_norm)
    };

if a_norm >= 1e-6
    msgs{end+1} = 'A mismatch details:'; %#ok<AGROW>
    msgs = [msgs; elementwise_diff('A', a_err, 1e-9)]; %#ok<AGROW>
end
if b_norm >= 1e-6
    msgs{end+1} = 'B mismatch details:'; %#ok<AGROW>
    msgs = [msgs; elementwise_diff('B', b_err, 1e-9)]; %#ok<AGROW>
end

pass = a_norm < 1e-6 && b_norm < 1e-6;
end


function [pass, msgs] = stage2_signs(cfg)
msgs = {};
pass = true;

xdot_eq = local_dynamics(zeros(5, 1), zeros(2, 1), cfg);
eq_pass = norm(xdot_eq) < 1e-12;
msgs{end+1} = sprintf('equilibrium xdot norm = %.3e -> %s', norm(xdot_eq), pass_fail_text(eq_pass)); %#ok<AGROW>
pass = pass && eq_pass;

xdot_pitch = local_dynamics([0.1; 0; 0; 0; 0], zeros(2, 1), cfg);
pitch_grav_pass = xdot_pitch(2) > 0;
msgs{end+1} = sprintf('pitch gravity sign: xdot(2)=%.6f -> %s', xdot_pitch(2), pass_fail_text(pitch_grav_pass)); %#ok<AGROW>
pass = pass && pitch_grav_pass;

xdot_roll = local_dynamics([0; 0; 0.1; 0; 0], zeros(2, 1), cfg);
roll_grav_pass = xdot_roll(4) > 0;
msgs{end+1} = sprintf('roll gravity sign:  xdot(4)=%.6f -> %s', xdot_roll(4), pass_fail_text(roll_grav_pass)); %#ok<AGROW>
pass = pass && roll_grav_pass;

xdot_rw = local_dynamics(zeros(5, 1), [0.1; 0], cfg);
rw_roll_pass = xdot_rw(4) < 0;
rw_spin_pass = xdot_rw(5) > 0;
msgs{end+1} = sprintf('tau_rw sign: xdot(4)=%.6f -> %s', xdot_rw(4), pass_fail_text(rw_roll_pass)); %#ok<AGROW>
msgs{end+1} = sprintf('tau_rw wheel spin: xdot(5)=%.6f -> %s', xdot_rw(5), pass_fail_text(rw_spin_pass)); %#ok<AGROW>
pass = pass && rw_roll_pass && rw_spin_pass;

xdot_base = local_dynamics(zeros(5, 1), [0; 0.1], cfg);
base_pass = xdot_base(2) < 0;
msgs{end+1} = sprintf('tau_base sign: xdot(2)=%.6f -> %s', xdot_base(2), pass_fail_text(base_pass)); %#ok<AGROW>
pass = pass && base_pass;
end


function [pass, msgs] = stage3_lin_vs_nl(A, B, K, cfg)
dt = cfg.dt;
T_check = 2.0;
N = round(T_check / dt);
x0 = [0.05; 0; 0.03; 0; 0];
t = (0:N) * dt;
x_lin = zeros(5, N + 1);
x_nl = zeros(5, N + 1);
x_lin(:, 1) = x0;
x_nl(:, 1) = x0;

for k = 1:N
    u_lin = saturate_control(-K * x_lin(:, k), cfg);
    u_nl = saturate_control(-K * x_nl(:, k), cfg);
    x_lin(:, k + 1) = x_lin(:, k) + dt * (A * x_lin(:, k) + B * u_lin);
    x_nl(:, k + 1) = x_nl(:, k) + dt * local_dynamics(x_nl(:, k), u_nl, cfg);
end

pitch_err_deg = max(abs(rad2deg(x_lin(1, :) - x_nl(1, :))));
roll_err_deg = max(abs(rad2deg(x_lin(3, :) - x_nl(3, :))));
sample_times = [0.0, 0.5, 1.0, 2.0];
pitch_lin_samples = sample_values_deg(x_lin(1, :), t, sample_times);
pitch_nl_samples = sample_values_deg(x_nl(1, :), t, sample_times);

pitch_pass = pitch_err_deg < 1.0;
roll_pass = roll_err_deg < 1.0;
pass = pitch_pass && roll_pass;

msgs = {
    sprintf('max pitch error = %.6f deg -> %s', pitch_err_deg, pass_fail_text(pitch_pass))
    sprintf('max roll error  = %.6f deg -> %s', roll_err_deg, pass_fail_text(roll_pass))
    'pitch trajectory samples [deg]:'
    };
for i = 1:numel(sample_times)
    msgs{end+1} = sprintf('  t=%3.1f  linear=%8.3f  nonlinear=%8.3f', ...
        sample_times(i), pitch_lin_samples(i), pitch_nl_samples(i)); %#ok<AGROW>
end
end


function [pass, msgs] = stage4_discretisation(A, B, Ad, Bd, cfg)
dt = cfg.dt;
Ad_ref = expm(A * dt);
Bd_ref = simpson_zoh_integral(A, B, dt, 2001);
ad_err = norm(Ad - Ad_ref, 'fro');
bd_err = norm(Bd - Bd_ref, 'fro');
rho_ol = max(abs(eig(Ad)));

msgs = {
    sprintf('Ad Frobenius error = %.3e -> %s', ad_err, pass_fail_text(ad_err < 1e-10))
    sprintf('Bd Frobenius error = %.3e -> %s', bd_err, pass_fail_text(bd_err < 1e-10))
    sprintf('rcond(A) = %.3e (A is singular, so the A\\\\B shortcut is not used)', rcond(A))
    'Bd ='
    sprintf('  [% .6e % .6e]', Bd(1,1), Bd(1,2))
    sprintf('  [% .6e % .6e]', Bd(2,1), Bd(2,2))
    sprintf('  [% .6e % .6e]', Bd(3,1), Bd(3,2))
    sprintf('  [% .6e % .6e]', Bd(4,1), Bd(4,2))
    sprintf('  [% .6e % .6e]', Bd(5,1), Bd(5,2))
    sprintf('open-loop spectral radius rho(Ad) = %.6f -> %s', rho_ol, pass_fail_text(rho_ol > 1.0))
    };

sign_checks = {
    Bd(2,2) < 0, 'Bd(2,2) < 0 (tau_base reduces pitch rate)'
    Bd(4,1) < 0, 'Bd(4,1) < 0 (tau_rw reduces roll rate)'
    Bd(5,1) > 0, 'Bd(5,1) > 0 (tau_rw spins wheel up)'
    abs(Bd(2,1)) < 1e-12, 'Bd(2,1) ~ 0 (no pitch-rate coupling from tau_rw)'
    abs(Bd(4,2)) < 1e-12, 'Bd(4,2) ~ 0 (no roll-rate coupling from tau_base)'
    abs(Bd(5,2)) < 1e-12, 'Bd(5,2) ~ 0 (no wheel-speed coupling from tau_base)'
    };

pass = ad_err < 1e-10 && bd_err < 1e-10 && rho_ol > 1.0;
for i = 1:size(sign_checks, 1)
    ok = sign_checks{i, 1};
    label = sign_checks{i, 2};
    msgs{end+1} = sprintf('%s -> %s', label, pass_fail_text(ok)); %#ok<AGROW>
    pass = pass && ok;
end

% Same-axis angle entries are expected to be nonzero under ZOH because
% one sample of acceleration integrates into one sample of angle change.
msgs{end+1} = sprintf('Bd(1,2)=%.6e and Bd(3,1)=%.6e are expected same-axis angle terms under ZOH.', ...
    Bd(1,2), Bd(3,1)); %#ok<AGROW>
end


function [pass, msgs] = stage5_step_signs(Ad, Bd)
x_after_rw = Ad * zeros(5,1) + Bd * [1; 0];
x_after_base = Ad * zeros(5,1) + Bd * [0; 1];

msgs = {
    sprintf('x_after_rw   = [% .6e % .6e % .6e % .6e % .6e]', x_after_rw)
    sprintf('x_after_base = [% .6e % .6e % .6e % .6e % .6e]', x_after_base)
    };

checks = {
    x_after_rw(4) < 0, 'tau_rw impulse gives x_after_rw(4) < 0'
    x_after_rw(5) > 0, 'tau_rw impulse gives x_after_rw(5) > 0'
    abs(x_after_rw(2)) < 1e-8, 'tau_rw impulse keeps x_after_rw(2) ~ 0'
    abs(x_after_rw(1)) < 1e-8, 'tau_rw impulse keeps x_after_rw(1) ~ 0'
    x_after_base(2) < 0, 'tau_base impulse gives x_after_base(2) < 0'
    abs(x_after_base(4)) < 1e-8, 'tau_base impulse keeps x_after_base(4) ~ 0'
    abs(x_after_base(3)) < 1e-8, 'tau_base impulse keeps x_after_base(3) ~ 0'
    abs(x_after_base(5)) < 1e-8, 'tau_base impulse keeps x_after_base(5) ~ 0'
    };

pass = true;
for i = 1:size(checks, 1)
    ok = checks{i, 1};
    msgs{end+1} = sprintf('%s -> %s', checks{i, 2}, pass_fail_text(ok)); %#ok<AGROW>
    pass = pass && ok;
end
end


function [pass, msgs] = stage6_open_closed_loop(K, cfg)
dt = cfg.dt;
N = round(3.0 / dt);
t = (0:N) * dt;
x0 = [0.05; 0; 0.03; 0; 0];

x_ol = zeros(5, N + 1);
x_cl = zeros(5, N + 1);
x_ol(:, 1) = x0;
x_cl(:, 1) = x0;
crashed_cl = false;

for k = 1:N
    x_ol(:, k + 1) = x_ol(:, k) + dt * local_dynamics(x_ol(:, k), [0; 0], cfg);
    u = -K * x_cl(:, k);
    x_cl(:, k + 1) = x_cl(:, k) + dt * local_dynamics(x_cl(:, k), u, cfg);
    if abs(x_cl(1, k + 1)) > cfg.pitch_crash_rad || abs(x_cl(3, k + 1)) > cfg.roll_crash_rad
        crashed_cl = true;
    end
end

ol_diverged = max(abs(x_ol(1, :))) > cfg.pitch_crash_rad || max(abs(x_ol(3, :))) > cfg.roll_crash_rad;
cl_small = norm(x_cl([1, 3], end)) < deg2rad(2);
cl_not_large = norm(x_cl([1, 3], end)) <= deg2rad(5);

sample_times = [0.0, 0.5, 1.0, 2.0, 3.0];
ol_samples = sample_values_deg(x_ol(1,:), t, sample_times);
cl_samples = sample_values_deg(x_cl(1,:), t, sample_times);

msgs = {
    sprintf('open-loop diverges past crash threshold -> %s', pass_fail_text(ol_diverged))
    sprintf('closed-loop crashed before 3 s -> %s', pass_fail_text(~crashed_cl))
    sprintf('closed-loop final angle norm = %.6f rad -> %s', norm(x_cl([1, 3], end)), pass_fail_text(cl_small))
    'pitch trajectory samples [deg]:'
    };
for i = 1:numel(sample_times)
    msgs{end+1} = sprintf('  t=%3.1f  open-loop=%8.3f  closed-loop=%8.3f', ...
        sample_times(i), ol_samples(i), cl_samples(i)); %#ok<AGROW>
end

pass = ol_diverged && ~crashed_cl && cl_small && cl_not_large;
end


function [pass, msgs] = stage7_axis_decoupling(cfg)
dt = cfg.dt;
N = round(1.0 / dt);

x_pitch = zeros(5, N + 1);
x_roll = zeros(5, N + 1);
x_pitch(:, 1) = [0.05; 0; 0; 0; 0];
x_roll(:, 1) = [0; 0; 0.05; 0; 0];

for k = 1:N
    x_pitch(:, k + 1) = x_pitch(:, k) + dt * local_dynamics(x_pitch(:, k), [0; 0], cfg);
    x_roll(:, k + 1) = x_roll(:, k) + dt * local_dynamics(x_roll(:, k), [0; 0], cfg);
end

max_roll_from_pitch_deg = max(abs(rad2deg(x_pitch(3, :))));
max_pitch_from_roll_deg = max(abs(rad2deg(x_roll(1, :))));
pitch_only_pass = max_roll_from_pitch_deg < 0.01;
roll_only_pass = max_pitch_from_roll_deg < 0.01;

msgs = {
    sprintf('pitch-only excitation max roll = %.6f deg -> %s', max_roll_from_pitch_deg, pass_fail_text(pitch_only_pass))
    sprintf('roll-only excitation max pitch = %.6f deg -> %s', max_pitch_from_roll_deg, pass_fail_text(roll_only_pass))
    };
pass = pitch_only_pass && roll_only_pass;
end


function [pass, msgs] = stage8_kalman_audit(Ad, Bd, C, K, L, cfg)
obs_poles = eig(Ad - L * C * Ad);
obs_pass = all(abs(obs_poles) < 1.0);

t = 5.0;
x0 = [0.05; 0; 0.03; 0; 0];
nominal = run_closed_loop_kalman(Ad, Bd, C, K, L, cfg, x0, t, 0.0);
bias_mag = 3 * cfg.noise_pitch_std;
biased = run_closed_loop_kalman(Ad, Bd, C, K, L, cfg, x0, t, bias_mag);

idx_after_1s = nominal.t > 1.0;
pitch_rms = rms_manual(nominal.x_true(1, idx_after_1s) - nominal.x_est(1, idx_after_1s));
roll_rms = rms_manual(nominal.x_true(3, idx_after_1s) - nominal.x_est(3, idx_after_1s));
bias_max = max(abs(biased.x_true(1, idx_after_1s) - biased.x_est(1, idx_after_1s)));
sample_times = 0.5:0.5:t;
pitch_err_samples = sample_series(nominal.x_true(1, :) - nominal.x_est(1, :), nominal.t, sample_times);

pitch_pass = pitch_rms < deg2rad(1);
roll_pass = roll_rms < deg2rad(1);
bias_pass = bias_max < deg2rad(3.0);

msgs = {
    sprintf('observer pole magnitudes = [%s] -> %s', sprintf('%.6f ', abs(obs_poles)), pass_fail_text(obs_pass))
    sprintf('pitch RMS after 1 s = %.6f rad -> %s', pitch_rms, pass_fail_text(pitch_pass))
    sprintf('roll RMS after 1 s  = %.6f rad -> %s', roll_rms, pass_fail_text(roll_pass))
    'pitch estimation error samples [rad]:'
    };
for i = 1:numel(sample_times)
    msgs{end+1} = sprintf('  t=%3.1f  err=% .6f', sample_times(i), pitch_err_samples(i)); %#ok<AGROW>
end
msgs{end+1} = sprintf('max pitch bias error after 1 s with +%.6f rad meas bias = %.6f rad -> %s', ...
    bias_mag, bias_max, pass_fail_text(bias_pass)); %#ok<AGROW>
pass = obs_pass && pitch_pass && roll_pass && bias_pass;
end


function xdot = local_dynamics(x, u, cfg)
theta_p = x(1);
dtheta_p = x(2); %#ok<NASGU>
theta_r = x(3);
dtheta_r = x(4); %#ok<NASGU>

tau_rw = u(1);
tau_base = u(2);

m = cfg.m_body + cfg.payload_mass;
L = cfg.L_body;
g = cfg.g;
Ib_p = cfg.I_body + m * L^2;
Ib_r = cfg.I_body_r;
Irw = cfg.I_rw;

ddtheta_p = (m * g * L * sin(theta_p) - tau_base) / Ib_p;
ddtheta_r = (m * g * L * sin(theta_r) - tau_rw) / Ib_r;
domega_rw = tau_rw / Irw;

xdot = [x(2); ddtheta_p; x(4); ddtheta_r; domega_rw];
end


function sim = run_closed_loop_kalman(Ad, Bd, C, K, L, cfg, x0, T_sim, pitch_bias)
dt = cfg.dt;
N = round(T_sim / dt);
t = (0:N) * dt;
ny = size(C, 1);
noise_std = measurement_noise_std(ny, cfg);

x_true = zeros(5, N + 1);
x_est = zeros(5, N + 1);
u_hist = zeros(2, N);
x_true(:, 1) = x0;

for k = 1:N
    u = saturate_control(-K * x_est(:, k), cfg);
    u_hist(:, k) = u;
    x_true(:, k + 1) = x_true(:, k) + dt * local_dynamics(x_true(:, k), u, cfg);

    y = C * x_true(:, k + 1) + noise_std .* randn(ny, 1);
    y(1) = y(1) + pitch_bias;

    x_pred = Ad * x_est(:, k) + Bd * u;
    x_est(:, k + 1) = x_pred + L * (y - C * x_pred);
end

sim.t = t;
sim.x_true = x_true;
sim.x_est = x_est;
sim.u = u_hist;
end


function u_sat = saturate_control(u, cfg)
u_sat = u;
u_sat(1) = max(min(u_sat(1), cfg.tau_rw_max), -cfg.tau_rw_max);
u_sat(2) = max(min(u_sat(2), cfg.tau_base_max), -cfg.tau_base_max);
end


function Bd_ref = simpson_zoh_integral(A, B, dt, n_pts)
if mod(n_pts, 2) == 0
    error('simpson_zoh_integral requires an odd number of points.');
end

tau = linspace(0, dt, n_pts);
h = tau(2) - tau(1);
Bd_ref = zeros(size(B));
for i = 1:n_pts
    weight = 2;
    if i == 1 || i == n_pts
        weight = 1;
    elseif mod(i, 2) == 0
        weight = 4;
    end
    Bd_ref = Bd_ref + weight * (expm(A * tau(i)) * B);
end
Bd_ref = Bd_ref * (h / 3);
end


function values_deg = sample_values_deg(signal, t, sample_times)
values_deg = zeros(1, numel(sample_times));
for i = 1:numel(sample_times)
    [~, idx] = min(abs(t - sample_times(i)));
    values_deg(i) = rad2deg(signal(idx));
end
end


function values = sample_series(signal, t, sample_times)
values = zeros(1, numel(sample_times));
for i = 1:numel(sample_times)
    [~, idx] = min(abs(t - sample_times(i)));
    values(i) = signal(idx);
end
end


function values = measurement_noise_std(ny, cfg)
base = [cfg.noise_pitch_std;
        cfg.noise_roll_std;
        cfg.noise_dpitch_std;
        cfg.noise_droll_std;
        cfg.noise_omega_rw_std];
values = base(1:ny);
end


function r = rms_manual(x)
r = sqrt(mean(x(:).^2));
end


function msgs = elementwise_diff(label, err_mat, tol)
[rows, cols] = find(abs(err_mat) > tol);
msgs = cell(numel(rows), 1);
for i = 1:numel(rows)
    msgs{i} = sprintf('  %s(%d,%d): %.3e', label, rows(i), cols(i), err_mat(rows(i), cols(i)));
end
if isempty(msgs)
    msgs = {sprintf('  no %s elements exceeded %.1e', label, tol)};
end
end


function print_stage_messages(msgs)
for i = 1:numel(msgs)
    fprintf('%s\n', msgs{i});
end
end


function txt = pass_fail_text(tf)
if tf
    txt = 'PASS';
else
    txt = 'FAIL';
end
end


function txt = stage_banner(stage_idx, tf)
txt = sprintf('Stage %d result: %s', stage_idx, pass_fail_text(tf));
end
