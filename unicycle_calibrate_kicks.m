function unicycle_calibrate_kicks()
%UNICYCLE_CALIBRATE_KICKS  Sweep kick amplitude using true-state feedback.
%  This isolates disturbance margin from estimator noise by bypassing the
%  Kalman filter and driving the nonlinear plant directly with true state.

clc;
fprintf('=== UNICYCLE KICK CALIBRATION ===\n\n');

cfg_sweep = unicycle_config();
cfg_sweep.x0_nominal = [0.04; 0; 0.03; 0; 0];
cfg_sweep.export_telemetry = false;
if isfield(cfg_sweep, 'dist_sine_amp')
    cfg_sweep.dist_sine_amp = 0.0;
end

[A, B, C, ~] = unicycle_plant(cfg_sweep); %#ok<ASGLU>
ctrl = unicycle_design_controllers(A, B, C, cfg_sweep);

amps = 0.10:0.02:0.35;
fprintf('kick schedule = [%s] s\n', sprintf('%.1f ', cfg_sweep.dist_kicks));
fprintf('kick duration = %.3f s\n\n', cfg_sweep.dist_duration);
fprintf('%-8s  %-12s  %-12s  %-10s\n', 'amp', 'ps_survival', 'mpc_score', 'ps_score');

amp_max_safe = NaN;
for amp = amps
    cfg_sweep.dist_amp = amp;

    res_ps = sweep_one(ctrl.paper_split, cfg_sweep);
    res_mpc = sweep_one(ctrl.baseline_mpc, cfg_sweep);

    fprintf('%.2f      %.3f        %.3f        %.3f\n', ...
        amp, res_ps.survival_frac, res_mpc.composite, res_ps.composite);

    if abs(res_ps.survival_frac - 1.0) < 1e-12
        amp_max_safe = amp;
    end
end

if isnan(amp_max_safe)
    fprintf('\namp_max_safe = none found in sweep range %.2f:%.2f:%.2f N*m\n', ...
        amps(1), amps(2) - amps(1), amps(end));
else
    fprintf('\namp_max_safe = %.2f N*m\n', amp_max_safe);
end
end


function res = sweep_one(ctrl, cfg)
dt = cfg.dt;
N = round(cfg.T_sim / dt);
x_true = zeros(5, N + 1);
u_hist = zeros(2, N);
disturbance = zeros(2, N);
x_true(:, 1) = cfg.x0_nominal;

xi_int = zeros(2, 1);
crashed = false;
crash_step = N;

kick_times = cfg.dist_kicks(:).';
kick_signs = ones(size(kick_times));
kick_signs(2:2:end) = -1;
kick_len = max(round(cfg.dist_duration / dt), 1);

for i = 1:numel(kick_times)
    kick_start = max(round(kick_times(i) / dt), 1);
    kick_end = min(kick_start + kick_len - 1, N);
    disturbance(1, kick_start:kick_end) = disturbance(1, kick_start:kick_end) ...
        + kick_signs(i) * cfg.dist_amp;
end

for k = 1:N
    xk = x_true(:, k);
    if strcmp(ctrl.name, 'Paper Split Baseline')
        u_pitch = -ctrl.K_pitch * xk([1, 2]);
        u_roll = -ctrl.K_roll * xk([3, 4, 5]);
        u = [u_roll; u_pitch];
    elseif isfield(ctrl, 'integrator') && ctrl.integrator
        xi_int = (1 - cfg.integrator_leak * dt) * xi_int + [xk(1); xk(3)] * dt;
        u = -ctrl.K * [xk; xi_int];
    else
        u = -ctrl.K * xk;
    end

    omega_rw = x_true(5, k);
    if abs(omega_rw) > cfg.omega_rw_hard
        u(1) = 0;
    elseif abs(omega_rw) > cfg.omega_rw_max
        scale = 1 - (abs(omega_rw) - cfg.omega_rw_max) / ...
            (cfg.omega_rw_hard - cfg.omega_rw_max);
        u(1) = u(1) * max(scale, 0);
    end

    u(1) = max(min(u(1), cfg.tau_rw_max), -cfg.tau_rw_max);
    u(2) = max(min(u(2), cfg.tau_base_max), -cfg.tau_base_max);
    u_hist(:, k) = u;

    xdot = nonlinear_dynamics(x_true(:, k), u, disturbance(:, k), cfg);
    x_true(:, k + 1) = x_true(:, k) + dt * xdot;

    if abs(x_true(1, k + 1)) > cfg.pitch_crash_rad || abs(x_true(3, k + 1)) > cfg.roll_crash_rad
        crashed = true;
        crash_step = k + 1;
        break;
    end
end

T_full = N * dt;
T_alive = crash_step * dt;
idx_u = max(crash_step - 1, 1);

res.survival_frac = T_alive / T_full;
res.pitch_rms = rms(x_true(1, 1:crash_step));
res.roll_rms = rms(x_true(3, 1:crash_step));
res.ctrl_rms = rms(u_hist(1, 1:idx_u));
res.composite = 100 * res.survival_frac ...
              - 50 * res.pitch_rms ...
              - 40 * res.roll_rms ...
              - 5 * res.ctrl_rms;
res.crashed = crashed;
end


function xdot = nonlinear_dynamics(x, u, disturbance, cfg)
m = cfg.m_body + cfg.payload_mass;
L = cfg.L_body;
Ib_p = cfg.I_body + m * L^2;
Ib_r = cfg.I_body_r;

ddtheta_p = (m * cfg.g * L * sin(x(1)) - u(2) + disturbance(2)) / Ib_p;
ddtheta_r = (m * cfg.g * L * sin(x(3)) - u(1) + disturbance(1)) / Ib_r;
domega_rw = u(1) / cfg.I_rw;

xdot = [x(2); ddtheta_p; x(4); ddtheta_r; domega_rw];
end
