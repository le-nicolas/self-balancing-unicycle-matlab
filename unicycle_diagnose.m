function report = unicycle_diagnose()
%UNICYCLE_DIAGNOSE  Structured diagnosis for benchmark parity and stress survival.
%
%  Run this before tuning. It inspects the current controller, estimator,
%  discretization, and disturbance configuration to identify systematic
%  failures that can make every controller family collapse in stress tests.

clc;
fprintf('=== UNICYCLE DIAGNOSE ===\n\n');

cfg = unicycle_config();
[A, B, C, ~] = unicycle_plant(cfg);
ctrl_all = unicycle_design_controllers(A, B, C, cfg);
[Ad, Bd] = unicycle_discretize(A, B, cfg.dt);

report = struct();
flags = {};

fprintf('--- PLANT CHECK --------------------------------------------------\n');
ev_open = eig(A);
disp_labeled('eig(A)', ev_open);
report.open_loop_eigs = ev_open;
report.n_rhp = sum(real(ev_open) > 0);
fprintf('RHP eigenvalues: %d\n', report.n_rhp);
if report.n_rhp ~= 2
    flags{end+1} = sprintf('plant has %d RHP eigenvalues instead of 2', report.n_rhp); %#ok<AGROW>
end

cond_ctrb = cond(ctrb(A, B));
cond_obsv = cond(obsv(A, C));
report.cond_ctrb = cond_ctrb;
report.cond_obsv = cond_obsv;
fprintf('cond(ctrb(A,B)) = %.3e\n', cond_ctrb);
fprintf('cond(obsv(A,C)) = %.3e\n', cond_obsv);

disp_labeled('Ad', Ad);
disp_labeled('Bd', Bd);
Ad_exact = expm(A * cfg.dt);
ad_err = norm(Ad - Ad_exact);
report.ad_exact_error = ad_err;
fprintf('norm(Ad - expm(A*dt)) = %.3e\n', ad_err);
if ad_err >= 1e-8
    flags{end+1} = sprintf('ZOH Ad mismatch is %.3e', ad_err); %#ok<AGROW>
end

families = {'lqr_current', 'lqr_dob', 'hybrid_modern', 'paper_split', 'baseline_mpc', 'robust_hinf_like'};
report.closed_loop = struct();
for i = 1:numel(families)
    fn = families{i};
    ctrl = ctrl_all.(fn);
    [K_ct, Kd_embedded] = effective_gains(ctrl, fn);
    ev_ct = eig(A - B * K_ct);
    sr_dt = max(abs(eig(Ad - Bd * Kd_embedded)));
    report.closed_loop.(fn).eig_ct = ev_ct;
    report.closed_loop.(fn).spectral_radius_dt = sr_dt;

    fprintf('\n[%s]\n', fn);
    disp_labeled('eig(A - B*K)', ev_ct);
    fprintf('spectral_radius(Ad - Bd*K) = %.6f\n', sr_dt);
    if any(real(ev_ct) > -0.1)
        msg = sprintf('%s has weak continuous-time poles (max real part %.3f)', ...
            fn, max(real(ev_ct)));
        fprintf('FLAG: %s\n', msg);
        flags{end+1} = msg; %#ok<AGROW>
    end
    if sr_dt < 0.960
        msg = sprintf('%s is over-aggressive (spectral radius %.6f)', fn, sr_dt);
        fprintf('FLAG: %s\n', msg);
        flags{end+1} = msg; %#ok<AGROW>
    elseif sr_dt > 0.995
        msg = sprintf('%s is weakly stable (spectral radius %.6f)', fn, sr_dt);
        fprintf('FLAG: %s\n', msg);
        flags{end+1} = msg; %#ok<AGROW>
    elseif sr_dt >= 1.0
        msg = sprintf('%s discrete closed loop is unstable (spectral radius %.6f)', fn, sr_dt);
        fprintf('FLAG: %s\n', msg);
        flags{end+1} = msg; %#ok<AGROW>
    end
end

fprintf('\n--- MPC CHECK ----------------------------------------------------\n');
mpc_ctrl = ctrl_all.baseline_mpc;
report.mpc = struct();
report.mpc.license_gate = mpc_ctrl.license_gate;
report.mpc.quadprog_available = mpc_ctrl.quadprog_available;
report.mpc.constrained = mpc_ctrl.constrained;
report.mpc.horizon = cfg.mpc_horizon;
fprintf('Optimization Toolbox present: %d\n', mpc_ctrl.constrained);
fprintf('license gate result:         %d\n', mpc_ctrl.license_gate);
fprintf('quadprog available:          %d\n', mpc_ctrl.quadprog_available);
fprintf('MPC horizon N:                %d steps (%.3f s)\n', ...
    cfg.mpc_horizon, cfg.mpc_horizon * cfg.dt);
if mpc_ctrl.constrained
    ev_H = eig(mpc_ctrl.H);
    min_eig_H = min(real(ev_H));
    report.mpc.min_eig_H = min_eig_H;
    fprintf('min eig(H):                   %.6f\n', min_eig_H);
    if min_eig_H <= 0
        msg = 'MPC Hessian H is not positive definite';
        fprintf('  FLAG: %s\n', msg);
        flags{end+1} = msg; %#ok<AGROW>
    else
        fprintf('  H is positive definite -> OK\n');
    end

    fprintf('lb(1) / tau_rw_max:           %.4f (should be -1.0)\n', ...
        mpc_ctrl.lb(1) / cfg.tau_rw_max);
    fprintf('ub(1) / tau_rw_max:           %.4f (should be +1.0)\n', ...
        mpc_ctrl.ub(1) / cfg.tau_rw_max);

    f_test = mpc_ctrl.F_half * cfg.x0_nominal;
    [~, ~, ef] = quadprog(mpc_ctrl.H, f_test, [], [], [], [], ...
        mpc_ctrl.lb, mpc_ctrl.ub, [], mpc_ctrl.qp_options);
    report.mpc.test_exitflag = ef;
    fprintf('Test QP solve exitflag:       %d (1 = success)\n', ef);
    if ef ~= 1
        msg = 'MPC test QP solve failed at nominal IC';
        fprintf('  FLAG: %s\n', msg);
        flags{end+1} = msg; %#ok<AGROW>
    end

    t_start = tic;
    for i = 1:100
        quadprog(mpc_ctrl.H, f_test, [], [], [], [], ...
            mpc_ctrl.lb, mpc_ctrl.ub, [], mpc_ctrl.qp_options);
    end
    t_per_solve_ms = toc(t_start) * 10;
    report.mpc.solve_time_ms = t_per_solve_ms;
    fprintf('QP solve time (mean):         %.3f ms\n', t_per_solve_ms);
    if t_per_solve_ms > 4.0
        fprintf('  WARNING: QP solve > 4 ms -- consider reducing N or precomputing\n');
    end
else
    fprintf('Constrained MPC inactive — dlqr fallback will be used.\n');
end

fprintf('\n--- KALMAN CHECK -------------------------------------------------\n');
L_lqe = lqe(A, eye(size(A, 1)), C, cfg.Q_kf, cfg.R_kf);
disp_labeled('L returned by lqe()', L_lqe);
report.L_lqe = L_lqe;

Qd = cfg.Q_kf * cfg.dt;
Rd = cfg.R_kf;
[~, ~, L_dare] = dare(Ad', C', Qd, Rd);
L_dare = L_dare';
disp_labeled('L used by discrete DARE design', L_dare);
report.L_dare = L_dare;
fprintf('norm(L_lqe - L_dare) = %.3e\n', norm(L_lqe - L_dare));

L_active = ctrl_all.lqr_current.L;
disp_labeled('L currently used by simulator', L_active);
report.L_active = L_active;

[rms_err, est_diverged] = estimator_only_check(Ad, Bd, C, L_active);
report.estimator_rms_err = rms_err;
report.estimator_diverged = est_diverged;
fprintf('Estimator-only RMS error by state = [%s]\n', sprintf('%.4e ', rms_err));
if est_diverged
    msg = 'estimator-only replay diverges with the active discrete observer gain';
    fprintf('FLAG: %s\n', msg);
    flags{end+1} = msg; %#ok<AGROW>
end

obs_poles = eig(Ad - L_active * C * Ad);
report.observer_poles = obs_poles;
disp_labeled('eig(Ad - L*C*Ad)', obs_poles);
if any(abs(obs_poles) >= 1.0)
    msg = sprintf('observer poles are unstable (max magnitude %.6f)', max(abs(obs_poles)));
    fprintf('FLAG: %s\n', msg);
    flags{end+1} = msg; %#ok<AGROW>
end

fprintf('\n--- WHEEL BUDGET CHECK -------------------------------------------\n');
cfg_budget = cfg;
cfg_budget.disturbance_enable = false;
cfg_budget.x0_nominal = [0.04; 0; 0.03; 0; 0];
budget = wheel_budget_check(A, B, C, ctrl_all.lqr_current, cfg_budget);
report.wheel_budget = budget;
fprintf('max |omega_rw|            = %.3f rad/s\n', budget.max_abs_omega_rw);
fprintf('soft suppression fraction = %.3f\n', budget.soft_fraction);
fprintf('hard cutoff fraction      = %.3f\n', budget.hard_fraction);
if budget.soft_fraction > 0.10
    msg = sprintf('wheel budget soft suppression fires in %.1f%% of steps', 100 * budget.soft_fraction);
    fprintf('FLAG: %s\n', msg);
    flags{end+1} = msg; %#ok<AGROW>
end

fprintf('\n--- DISTURBANCE CHECK --------------------------------------------\n');
if isfield(cfg, 'dist_kicks') && ~isempty(cfg.dist_kicks)
    kick_times = cfg.dist_kicks(:).';
else
    kick_times = [];
end
kick_signs = ones(size(kick_times));
if ~isempty(kick_signs)
    kick_signs(2:2:end) = -1;
end
active_fraction = numel(kick_times) * cfg.dist_duration / cfg.T_sim;
if isempty(kick_times)
    challenge_fraction = 0;
else
    challenge_fraction = ((kick_times(end) + cfg.dist_duration) - kick_times(1)) / cfg.T_sim;
end
report.benchmark_disturbance = struct( ...
    'kick_times', kick_times, ...
    'kick_signs', kick_signs, ...
    'dist_amp', cfg.dist_amp, ...
    'dist_duration', cfg.dist_duration, ...
    'active_fraction', active_fraction, ...
    'challenge_fraction', challenge_fraction, ...
    'benchmark_seed', cfg.benchmark_seed);
fprintf('benchmark dist_kicks    = [%s] s\n', sprintf('%.2f ', kick_times));
fprintf('benchmark kick signs    = [%s]\n', sprintf('%+d ', kick_signs));
fprintf('benchmark dist_amp      = %.3f N*m\n', cfg.dist_amp);
fprintf('benchmark dist_duration = %.3f s\n', cfg.dist_duration);
fprintf('benchmark seed          = %d\n', cfg.benchmark_seed);
fprintf('active kick fraction    = %.3f (sum of kick durations / T_sim)\n', active_fraction);
fprintf('challenge window frac   = %.3f ((last kick end - first kick) / T_sim)\n', challenge_fraction);
if challenge_fraction <= 0.30
    msg = sprintf('benchmark challenge window is too short (%.1f%% of episode)', ...
        100 * challenge_fraction);
    fprintf('FLAG: %s\n', msg);
    flags{end+1} = msg; %#ok<AGROW>
end

fprintf('\n[Sine window]\n');
if isfield(cfg, 'dist_sine_start') && isfield(cfg, 'dist_sine_end') ...
        && isfield(cfg, 'dist_sine_amp') && isfield(cfg, 'dist_sine_freq')
    omega_sine = 2 * pi * cfg.dist_sine_freq;
    pole_ps = 3.2;
    pole_lqr = 6.6;
    lag_ps = rad2deg(atan(omega_sine / pole_ps));
    lag_lqr = rad2deg(atan(omega_sine / pole_lqr));
    lag_diff = lag_ps - lag_lqr;

    report.benchmark_disturbance.dist_sine_start = cfg.dist_sine_start;
    report.benchmark_disturbance.dist_sine_end = cfg.dist_sine_end;
    report.benchmark_disturbance.dist_sine_amp = cfg.dist_sine_amp;
    report.benchmark_disturbance.dist_sine_freq = cfg.dist_sine_freq;
    report.benchmark_disturbance.sine_duration = cfg.dist_sine_end - cfg.dist_sine_start;
    report.benchmark_disturbance.phase_lag_paper_split = lag_ps;
    report.benchmark_disturbance.phase_lag_lqr_current = lag_lqr;
    report.benchmark_disturbance.phase_lag_difference = lag_diff;

    fprintf('dist_sine_start        = %.1f s\n', cfg.dist_sine_start);
    fprintf('dist_sine_end          = %.1f s\n', cfg.dist_sine_end);
    fprintf('dist_sine_amp          = %.4f N*m\n', cfg.dist_sine_amp);
    fprintf('dist_sine_freq         = %.1f Hz\n', cfg.dist_sine_freq);
    fprintf('sine window duration   = %.1f s\n', cfg.dist_sine_end - cfg.dist_sine_start);
    fprintf('phase lag paper_split  = %.1f deg\n', lag_ps);
    fprintf('phase lag lqr_current  = %.1f deg\n', lag_lqr);
    fprintf('lag difference         = %.1f deg\n', lag_diff);
    if lag_diff < 10
        fprintf('  WARNING: lag difference < 10 deg -- sine freq may be too low\n');
        fprintf('           consider increasing dist_sine_freq to 1.0-1.2 Hz\n');
    end
    if cfg.dist_sine_amp > 0.08
        fprintf('  WARNING: sine amp > 0.08 N*m -- survival risk for paper_split\n');
    end
end

fprintf('\n[Stress benchmark reference]\n');
stress_cfg = unicycle_config();
stress_cfg.dist_amp = 0.30;
stress_cfg.dist_duration = 0.05;
dist_impulse = stress_cfg.dist_amp * stress_cfg.dist_duration;
dist_amp_max = 1.3 * stress_cfg.dist_amp;
dist_impulse_max = dist_amp_max * stress_cfg.dist_duration;
roll_step_impulse = stress_cfg.tau_rw_max * stress_cfg.dt;
base_step_impulse = stress_cfg.tau_base_max * stress_cfg.dt;
report.stress_disturbance = struct( ...
    'dist_amp', stress_cfg.dist_amp, ...
    'dist_duration', stress_cfg.dist_duration, ...
    'impulse', dist_impulse, ...
    'dist_amp_max', dist_amp_max, ...
    'impulse_max', dist_impulse_max, ...
    'tau_rw_max', stress_cfg.tau_rw_max, ...
    'tau_base_max', stress_cfg.tau_base_max, ...
    'rw_step_impulse', roll_step_impulse, ...
    'base_step_impulse', base_step_impulse);
fprintf('dist_amp               = %.3f N*m\n', stress_cfg.dist_amp);
fprintf('dist_duration          = %.3f s\n', stress_cfg.dist_duration);
fprintf('impulse magnitude      = %.6f N*m*s\n', dist_impulse);
fprintf('max randomized amp     = %.3f N*m\n', dist_amp_max);
fprintf('max randomized impulse = %.6f N*m*s\n', dist_impulse_max);
fprintf('tau_rw_max             = %.3f N*m\n', stress_cfg.tau_rw_max);
fprintf('tau_base_max           = %.3f N*m\n', stress_cfg.tau_base_max);
fprintf('tau_rw_max * dt        = %.6f N*m*s per step\n', roll_step_impulse);
fprintf('tau_base_max * dt      = %.6f N*m*s per step\n', base_step_impulse);
if dist_impulse_max > 20 * roll_step_impulse
    msg = sprintf('stress disturbance impulse is %.1fx the reaction-wheel max per-step impulse', ...
        dist_impulse_max / roll_step_impulse);
    fprintf('FLAG: %s\n', msg);
    flags{end+1} = msg; %#ok<AGROW>
end

fprintf('\n--- NOISE CHECK --------------------------------------------------\n');
cfg_noise = cfg;
cfg_noise.disturbance_enable = false;
cfg_noise.x0_nominal = [0.04; 0; 0.03; 0; 0];
noise_stats = noise_control_check(A, B, C, ctrl_all.lqr_current, cfg_noise);
report.noise = noise_stats;
fprintf('max |u(1)|              = %.4f N*m\n', noise_stats.max_u1);
fprintf('max |u(2)|              = %.4f N*m\n', noise_stats.max_u2);
fprintf('reaction-wheel sat frac = %.3f\n', noise_stats.sat_frac_rw);
fprintf('base-wheel sat frac     = %.3f\n', noise_stats.sat_frac_base);
if noise_stats.sat_frac_rw > 0.05 || noise_stats.sat_frac_base > 0.05
    msg = sprintf('noise-only run saturates control too often (rw %.1f%%, base %.1f%%)', ...
        100 * noise_stats.sat_frac_rw, 100 * noise_stats.sat_frac_base);
    fprintf('FLAG: %s\n', msg);
    flags{end+1} = msg; %#ok<AGROW>
end

fprintf('\n--- SUMMARY ------------------------------------------------------\n');
if isempty(flags)
    likely_root_cause = 'No diagnostic flags were raised.';
else
    likely_root_cause = flags{1};
end
report.flags = flags;
report.likely_root_cause = likely_root_cause;
fprintf('LIKELY ROOT CAUSE: %s\n', likely_root_cause);
end


function [K_ct, Kd_embedded] = effective_gains(ctrl, family_name)
switch family_name
    case 'paper_split'
        K_ct = zeros(2, 5);
        K_ct(1, 3:5) = ctrl.K_roll;
        K_ct(2, 1:2) = ctrl.K_pitch;
        Kd_embedded = K_ct;
    case 'hybrid_modern'
        K_ct = ctrl.K(:, 1:5);
        Kd_embedded = ctrl.K(:, 1:5);
    otherwise
        K_ct = ctrl.K;
        Kd_embedded = ctrl.K;
end
end


function [rms_err, diverged] = estimator_only_check(Ad, ~, C, L)
nx = size(Ad, 1);
steps = 200;
x_true = zeros(nx, steps + 1);
x_est = zeros(nx, steps + 1);
x_true(:,1) = [0.03; 0; 0.02; 0; 0];
err_hist = zeros(nx, steps);

diverged = false;
for k = 1:steps
    x_true(:,k+1) = Ad * x_true(:,k);
    y = C * x_true(:,k+1);
    x_pred = Ad * x_est(:,k);
    innov = y - C * x_pred;
    x_est(:,k+1) = x_pred + L * innov;
    err_hist(:,k) = x_true(:,k+1) - x_est(:,k+1);
    if abs(err_hist(1,k)) > 0.5 || abs(err_hist(3,k)) > 0.5
        diverged = true;
    end
end

rms_err = sqrt(mean(err_hist.^2, 2));
end


function stats = wheel_budget_check(A, B, C, ctrl, cfg)
dt = cfg.dt;
N = cfg.N;
[Ad, Bd] = unicycle_discretize(A, B, dt);
x_true = zeros(5, N + 1);
x_est = zeros(5, N + 1);
x_true(:,1) = cfg.x0_nominal;

Qd = cfg.Q_kf * dt;
Rd = cfg.R_kf;
[~, ~, Lkf] = dare(Ad', C', Qd, Rd);
Lkf = Lkf';

soft_count = 0;
hard_count = 0;
max_abs_omega = 0;
u_prev = zeros(2, 1);

for k = 1:N
    y = C * x_true(:,k);
    x_pred = Ad * x_est(:,k) + Bd * u_prev;
    x_est(:,k+1) = x_pred + Lkf * (y - C * x_pred);
    xe = x_est(:,k+1);

    u = -ctrl.K * xe;
    omega_true = x_true(5,k);
    max_abs_omega = max(max_abs_omega, abs(omega_true));
    if abs(omega_true) > cfg.omega_rw_hard
        u(1) = 0;
        hard_count = hard_count + 1;
    elseif abs(omega_true) > cfg.omega_rw_max
        scale = 1 - (abs(omega_true) - cfg.omega_rw_max) / ...
            (cfg.omega_rw_hard - cfg.omega_rw_max);
        u(1) = u(1) * max(scale, 0);
        soft_count = soft_count + 1;
    end

    u(1) = max(min(u(1), cfg.tau_rw_max), -cfg.tau_rw_max);
    u(2) = max(min(u(2), cfg.tau_base_max), -cfg.tau_base_max);

    xdot = nonlinear_dynamics(x_true(:,k), u, [0; 0], cfg);
    x_true(:,k+1) = x_true(:,k) + dt * xdot;
    u_prev = u;
end

stats.max_abs_omega_rw = max_abs_omega;
stats.soft_fraction = soft_count / N;
stats.hard_fraction = hard_count / N;
end


function stats = noise_control_check(A, B, C, ctrl, cfg)
res = unicycle_simulate(A, B, C, ctrl, cfg);
u = res.u;
stats.max_u1 = max(abs(u(1,:)));
stats.max_u2 = max(abs(u(2,:)));
stats.sat_frac_rw = mean(abs(u(1,:)) >= cfg.tau_rw_max - 1e-12);
stats.sat_frac_base = mean(abs(u(2,:)) >= cfg.tau_base_max - 1e-12);
end


function xdot = nonlinear_dynamics(x, u, disturbance, cfg)
theta_p = x(1);
theta_r = x(3);
m = cfg.m_body + cfg.payload_mass;
L = cfg.L_body;
Ib_p = cfg.I_body + m * L^2;
Ib_r = cfg.I_body_r;
Irw = cfg.I_rw;

ddtheta_p = (m * cfg.g * L * sin(theta_p) - u(2) + disturbance(2)) / Ib_p;
ddtheta_r = (m * cfg.g * L * sin(theta_r) - u(1) + disturbance(1)) / Ib_r;
domega_rw = u(1) / Irw;
xdot = [x(2); ddtheta_p; x(4); ddtheta_r; domega_rw];
end


function disp_labeled(label, value)
fprintf('%s =\n', label);
disp(value);
end
