function res = unicycle_simulate(A, B, C, ctrl, cfg)
%UNICYCLE_SIMULATE  Closed-loop simulation of the unicycle robot.
%
%  Features:
%    - Discrete-time Kalman filter (state estimator)
%    - Sensor noise + IMU delay buffer
%    - Disturbance observer (if ctrl.dob = true)
%    - Integrator augmentation (if ctrl.integrator = true)
%    - Actuator saturation
%    - Reaction-wheel speed budget / hard cutoff
%    - Crash detection (fall > threshold)
%    - Disturbance impulse injection

dt = cfg.dt;
nx = 5;
ny = size(C, 1);
nu = 2;

[kick_times, kick_signs, T_sim] = disturbance_profile(cfg);
N = round(T_sim / dt);
single_kick_mode = ~isempty(kick_times) && all(kick_signs == 1) ...
    && isfield(cfg, 'dist_t_kick') && ~isempty(cfg.dist_t_kick);

if ~single_kick_mode && isfield(cfg, 'benchmark_seed') && ~isempty(cfg.benchmark_seed)
    prev_rng = rng; %#ok<NASGU>
    cleanup_rng = onCleanup(@() rng(prev_rng)); %#ok<NASGU>
    rng(cfg.benchmark_seed);
end

%% ── Discrete-time plant ─────────────────────────────────────
[Ad, Bd] = unicycle_discretize(A, B, dt);

% Discrete Kalman
if isfield(ctrl, 'L') && ~isempty(ctrl.L)
    Lkf = ctrl.L;
else
    Qd = cfg.Q_kf * dt;
    Rd = cfg.R_kf;
    [~, ~, Lkf] = dare(Ad', C', Qd, Rd);
    Lkf = Lkf';
end

%% ── Pre-allocate ─────────────────────────────────────────────
x_true = zeros(nx, N+1);
x_est  = zeros(nx, N+1);
u_hist = zeros(nu, N);
y_hist = zeros(ny, N);
y_meas = zeros(ny, N);
disturbance = zeros(nu, N);

x_true(:,1) = cfg.x0_nominal;
x_est(:,1)  = zeros(nx,1);   % estimator starts at zero

% Integrator state (pitch and roll integrals)
xi_int = zeros(2,1);

% DOB state
d_est   = zeros(nu,1);
u_prev  = zeros(nu,1);

% Optional online pitch identification for hybrid_modern
online_id_active = strcmp(ctrl.name, 'Hybrid Modern (LQR+I)') && ...
    isfield(cfg, 'enable_online_id') && cfg.enable_online_id;
if online_id_active
    rls = unicycle_rls('init', A, B, cfg);
    xe_rls_prev = x_est(:,1);
else
    rls = struct('theta_init', [], 'theta', [], ...
        'update_count', 0, 'gain_update_count', 0);
    xe_rls_prev = zeros(nx,1);
end

% Wheel speed (tracked separately for budget logic)
omega_rw = 0;

% Crash flag
crashed = false;
crash_step = N;

%% ── Disturbance schedule ─────────────────────────────────────
if cfg.disturbance_enable
    kick_len = max(round(cfg.dist_duration / dt), 1);
    for i = 1:numel(kick_times)
        kick_start = max(round(kick_times(i) / dt), 1);
        kick_end = min(kick_start + kick_len - 1, N);
        disturbance(1, kick_start:kick_end) = disturbance(1, kick_start:kick_end) ...
            + kick_signs(i) * cfg.dist_amp;
    end

    if isfield(cfg, 'dist_sine_amp') && isfield(cfg, 'dist_sine_start') ...
            && isfield(cfg, 'dist_sine_end') && isfield(cfg, 'dist_sine_freq') ...
            && cfg.dist_sine_amp ~= 0
        for k = 1:N
            t_k = k * dt;
            if t_k >= cfg.dist_sine_start && t_k < cfg.dist_sine_end
                disturbance(1, k) = disturbance(1, k) + ...
                    cfg.dist_sine_amp * sin(2 * pi * cfg.dist_sine_freq * t_k);
            end
        end
    end
end

%% ── Simulation loop ──────────────────────────────────────────
for k = 1:N

    %% (a) Noisy sensor measurement
    v_std = [cfg.noise_pitch_std;
             cfg.noise_roll_std;
             cfg.noise_dpitch_std;
             cfg.noise_droll_std];
    if ny >= 5
        v_std(5,1) = cfg.noise_omega_rw_std;
    end
    v_noise = v_std(1:ny) .* randn(ny,1);
    y_clean = C * x_true(:,k);
    y_meas(:,k) = y_clean + v_noise;

    %% (b) Fixed-lag update for delayed measurements
    meas_idx = k - cfg.imu_delay_samples;
    if meas_idx >= 1
        y_delayed = y_meas(:, meas_idx);
        y_hist(:,k) = y_delayed;

        x_replay = x_est(:, meas_idx) + Lkf * (y_delayed - C * x_est(:, meas_idx));
        x_est(:, meas_idx) = x_replay;
        for j = meas_idx:k-1
            x_replay = Ad * x_replay + Bd * u_hist(:,j);
            x_est(:, j+1) = x_replay;
        end
    end
    xe = x_est(:,k);

    if online_id_active
        if k > 1
            phi_rls = [xe_rls_prev(1); -u_prev(2)];
            y_rls = (xe(2) - xe_rls_prev(2)) / dt;
            rls = unicycle_rls('update', rls, phi_rls, y_rls, cfg);

            if rls.update_count >= cfg.online_id_min_updates && ...
                    mod(rls.update_count, cfg.online_id_recompute_every) == 0
                [A_est, B_est] = unicycle_rls('plant', rls, A, B);
                A_aug_est = zeros(7);
                A_aug_est(1:5,1:5) = A_est;
                A_aug_est(6,1) = 1;
                A_aug_est(6,6) = -cfg.integrator_leak;
                A_aug_est(7,3) = 1;
                A_aug_est(7,7) = -cfg.integrator_leak;
                B_aug_est = [B_est; zeros(2,2)];
                try
                    [Ad_aug_est, Bd_aug_est] = unicycle_discretize(A_aug_est, B_aug_est, dt);
                    ctrl.K = dlqr(Ad_aug_est, Bd_aug_est, ctrl.Q_aug, ctrl.R_aug);
                    rls.gain_update_count = rls.gain_update_count + 1;
                catch
                    % Keep the current gain if the estimated model is ill-conditioned.
                end
            end
        end
        xe_rls_prev = xe;
    end

    %% (d) Compute control
    if strcmp(ctrl.name, 'Paper Split Baseline')
        % Decoupled pitch / roll
        u_pitch = -ctrl.K_pitch * xe([1,2]);
        u_roll  = -ctrl.K_roll  * xe([3,4,5]);
        u = [u_roll; u_pitch];
    elseif strcmp(ctrl.name, 'Baseline MPC (unconstrained)') ...
            && isfield(ctrl, 'constrained') && ctrl.constrained
        u = solve_constrained_mpc(xe, ctrl, cfg);
    elseif ctrl.integrator
        % Augmented state
        xi_int = (1 - cfg.integrator_leak * dt) * xi_int + [xe(1); xe(3)] * dt;
        xe_aug = [xe; xi_int];
        u = -ctrl.K * xe_aug;
    else
        u = -ctrl.K * xe;
    end

    %% (e) Disturbance observer correction
    if isfield(ctrl, 'dob') && ctrl.dob
        alpha = ctrl.dob_alpha;
        xe_prev = x_est(:, max(k - 1, 1));
        m = cfg.m_body + cfg.payload_mass;
        roll_acc_nom = (m * cfg.g * cfg.L_body * sin(xe_prev(3)) - u_prev(1)) / cfg.I_body_r;
        droll_pred = xe_prev(4) + dt * roll_acc_nom;
        tau_roll_residual = cfg.I_body_r * (xe(4) - droll_pred) / dt;
        d_est(1) = (1 - alpha) * d_est(1) + alpha * tau_roll_residual;
        u(1) = u(1) + d_est(1);     % roll-axis disturbance compensation
    end

    %% (f) Wheel-speed budget logic (mirrors MuJoCo config)
    omega_rw = x_true(5,k);   % use true wheel speed for safety logic
    if abs(omega_rw) > cfg.omega_rw_hard
        u(1) = 0;       % hard cutoff: suppress torque same direction
    elseif abs(omega_rw) > cfg.omega_rw_max
        scale = 1 - (abs(omega_rw) - cfg.omega_rw_max) / ...
                    (cfg.omega_rw_hard - cfg.omega_rw_max);
        u(1) = u(1) * max(scale, 0);
    end

    %% (g) Actuator saturation
    u(1) = max(min(u(1), cfg.tau_rw_max),   -cfg.tau_rw_max);
    u(2) = max(min(u(2), cfg.tau_base_max), -cfg.tau_base_max);

    u_hist(:,k) = u;
    u_prev = u;

    %% (h) Plant dynamics (nonlinear Euler integration)
    w_proc = sqrt(diag(cfg.Q_kf)) .* randn(nx,1) * sqrt(dt);
    xdot = unicycle_dynamics(x_true(:,k), u, disturbance(:,k), cfg);
    x_true(:,k+1) = x_true(:,k) + xdot * dt + w_proc;
    x_est(:,k+1) = Ad * x_est(:,k) + Bd * u;

    %% (i) Crash detection
    if ~crashed
        theta_p = x_true(1,k+1);
        theta_r = x_true(3,k+1);
        if abs(theta_p) > cfg.pitch_crash_rad || abs(theta_r) > cfg.roll_crash_rad
            crashed  = true;
            crash_step = k+1;
            fprintf('  [!] CRASH at t=%.2f s  (θ_p=%.1f°, θ_r=%.1f°)\n', ...
                k*dt, rad2deg(theta_p), rad2deg(theta_r));
            break;
        end
    end
end

%% ── Pack results ─────────────────────────────────────────────
t = (0:N) * dt;
res.t         = t;
res.x_true    = x_true;
res.x_est     = x_est;
res.u         = u_hist;
res.y         = y_hist;
res.crashed   = crashed;
res.crash_step = crash_step;
res.crash_t   = crash_step * dt;
res.survival  = ~crashed;
res.name      = ctrl.name;
res.T_sim     = T_sim;
res.online_id_enabled = online_id_active;
res.online_id_gain_updates = rls.gain_update_count;
res.online_id_update_count = rls.update_count;
res.rls_theta_init = rls.theta_init;
res.rls_theta_final = rls.theta;

%% ── Composite score (mirrors benchmark.py metric) ────────────
T_full   = N * dt;
T_alive  = crash_step * dt;
pitch_rms = rms(x_true(1, 1:crash_step));
roll_rms  = rms(x_true(3, 1:crash_step));
u_rms     = rms(u_hist(1, 1:crash_step-1));

res.survival_frac = T_alive / T_full;
res.pitch_rms     = pitch_rms;
res.roll_rms      = roll_rms;
res.ctrl_rms      = u_rms;
res.composite     = 100 * res.survival_frac ...
                  - 50 * pitch_rms ...
                  - 40 * roll_rms ...
                  - 5  * u_rms;

if isfield(cfg, 'export_telemetry') && cfg.export_telemetry
    write_telemetry_csv(res, ctrl, cfg);
end

fprintf('  Survival: %.0f%%  Composite: %.2f  (pitch RMS=%.4f rad)\n', ...
    res.survival_frac*100, res.composite, pitch_rms);
end


function [kick_times, kick_signs, T_sim] = disturbance_profile(cfg)
%DISTURBANCE_PROFILE  Resolve the disturbance schedule for this simulation.
%  Single-kick callers set cfg.dist_t_kick explicitly. The deterministic
%  benchmark uses cfg.dist_kicks with alternating signs.

has_single_kick = isfield(cfg, 'dist_t_kick') && ~isempty(cfg.dist_t_kick);
if has_single_kick
    kick_times = cfg.dist_t_kick(:).';
    kick_signs = ones(size(kick_times));
    if isfield(cfg, 'single_kick_T_sim') && ~isempty(cfg.single_kick_T_sim)
        T_sim = cfg.single_kick_T_sim;
    else
        T_sim = cfg.T_sim;
    end
else
    if isfield(cfg, 'dist_kicks') && ~isempty(cfg.dist_kicks)
        kick_times = cfg.dist_kicks(:).';
    else
        kick_times = [];
    end
    kick_signs = ones(size(kick_times));
    if ~isempty(kick_signs)
        kick_signs(2:2:end) = -1;
    end
    T_sim = cfg.T_sim;
end
end


function u = solve_constrained_mpc(xe_mpc, ctrl, cfg)
%SOLVE_CONSTRAINED_MPC  Solve the condensed finite-horizon QP.
u = -ctrl.K_fallback * xe_mpc;

R_step = ctrl.R_bar_block;
if abs(xe_mpc(5)) > cfg.omega_rw_max
    R_step(1,1) = R_step(1,1) + ctrl.rho_w * ...
        (abs(xe_mpc(5)) - cfg.omega_rw_max) / cfg.omega_rw_max;
end

delta_rw_penalty = R_step(1,1) - ctrl.R_bar_block(1,1);
H_qp = ctrl.H;
if delta_rw_penalty ~= 0
    H_qp(1,1) = H_qp(1,1) + 2 * delta_rw_penalty;
end
f_vec = ctrl.F_half * xe_mpc;

try
    [U_opt, ~, exitflag] = quadprog(H_qp, f_vec, [], [], [], [], ...
        ctrl.lb, ctrl.ub, [], ctrl.qp_options);
    if exitflag > 0 && ~isempty(U_opt)
        u = U_opt(1:2);
    end
catch
    % Keep the dlqr fallback when the QP solver is unavailable or fails.
end
end


function xdot = unicycle_dynamics(x, u, disturbance, cfg)
%UNICYCLE_DYNAMICS  Nonlinear equations of motion (sin model).
theta_p  = x(1);
dtheta_p = x(2);
theta_r  = x(3);
dtheta_r = x(4);

tau_rw   = u(1);
tau_base = u(2);
tau_roll_dist  = disturbance(1);
tau_pitch_dist = disturbance(2);

m  = cfg.m_body + cfg.payload_mass;
L  = cfg.L_body;
g  = cfg.g;
Ib_p = cfg.I_body   + m * L^2;
Ib_r = cfg.I_body_r;
Irw  = cfg.I_rw;

ddtheta_p = (m*g*L*sin(theta_p) - tau_base + tau_pitch_dist) / Ib_p;
ddtheta_r = (m*g*L*sin(theta_r) - tau_rw + tau_roll_dist)    / Ib_r;
domega_rw = tau_rw / Irw;

xdot = [dtheta_p; ddtheta_p; dtheta_r; ddtheta_r; domega_rw];
end


function write_telemetry_csv(res, ctrl, cfg)
%WRITE_TELEMETRY_CSV  Export HIL-friendly telemetry for the current run.
Nv = min(res.crash_step, size(res.x_true, 2));
Nu = min(max(Nv - 1, 1), size(res.u, 2));

t = res.t(1:Nv).';
pitch_deg = rad2deg(res.x_true(1,1:Nv)).';
roll_deg  = rad2deg(res.x_true(3,1:Nv)).';
dpitch    = res.x_true(2,1:Nv).';
droll     = res.x_true(4,1:Nv).';
omega_rw  = res.x_true(5,1:Nv).';

u_pad = zeros(2, Nv);
u_pad(:,1:Nu) = res.u(:,1:Nu);
if Nv > 1
    u_pad(:,Nv) = u_pad(:,max(Nu, 1));
end
tau_rw   = u_pad(1,:).';
tau_base = u_pad(2,:).';

telemetry = table(t, pitch_deg, roll_deg, dpitch, droll, omega_rw, tau_rw, tau_base);
if isfield(ctrl, 'key')
    fname = fullfile(cfg.results_dir, sprintf('telemetry_%s.csv', ctrl.key));
else
    fname = fullfile(cfg.results_dir, 'telemetry_controller.csv');
end
writetable(telemetry, fname);
end
