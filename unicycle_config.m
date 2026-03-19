function cfg = unicycle_config()
%UNICYCLE_CONFIG  All physical and simulation parameters.
%  Derived from final/final.xml and final/config.yaml in the
%  le-nicolas/Self-Balancing-Unicycle-Robot MuJoCo project.

%% ── Physical parameters ─────────────────────────────────────
% Robot body (stick / chassis)
cfg.m_body   = 0.85;        % [kg]   body mass
cfg.payload_mass = 0.0;     % [kg]   optional payload added at body CoM
cfg.L_body   = 0.28;        % [m]    CoM height above wheel axle
cfg.I_body   = 0.022;       % [kg·m²] body pitch moment of inertia
cfg.I_body_r = 0.018;       % [kg·m²] body roll moment of inertia

% Reaction wheel (flywheel)
cfg.m_rw     = 0.12;        % [kg]   reaction wheel mass
cfg.r_rw     = 0.055;       % [m]    reaction wheel radius
cfg.I_rw     = 0.5 * cfg.m_rw * cfg.r_rw^2;  % thin disk

% Base drive wheel
cfg.m_wheel  = 0.08;        % [kg]
cfg.r_wheel  = 0.045;       % [m]    base wheel radius
cfg.I_wheel  = 0.5 * cfg.m_wheel * cfg.r_wheel^2;

% Gravity
cfg.g = 9.81;               % [m/s²]

%% ── Actuator limits ─────────────────────────────────────────
cfg.tau_rw_max   = 0.25;    % [N·m]  max reaction wheel torque
cfg.tau_base_max = 0.40;    % [N·m]  max base wheel torque
cfg.omega_rw_max = 800;     % [rad/s] wheel speed budget
cfg.omega_rw_hard = 950;    % [rad/s] hard cutoff

%% ── Sensor noise (Kalman) ───────────────────────────────────
cfg.noise_pitch_std  = 0.005;   % [rad]
cfg.noise_roll_std   = 0.005;   % [rad]
cfg.noise_dpitch_std = 0.008;   % [rad/s]
cfg.noise_droll_std  = 0.008;   % [rad/s]
cfg.noise_omega_rw_std = 0.02;  % [rad/s] reaction-wheel encoder noise
cfg.noise_bias_std   = 0.002;   % [rad/s]  gyro bias drift std

% IMU delay (samples)
cfg.imu_delay_samples = 2;

%% ── Simulation settings ──────────────────────────────────────
cfg.dt       = 0.005;       % [s]  control period (200 Hz)
cfg.T_sim    = 20.0;        % [s]  deterministic benchmark duration
cfg.N        = round(cfg.T_sim / cfg.dt);
cfg.single_kick_T_sim = 10.0;  % [s] preserve legacy single-kick callers
cfg.benchmark_seed = 1234;     % deterministic nominal benchmark noise

%% ── Initial conditions ───────────────────────────────────────
% State vector: [pitch, dpitch, roll, droll, omega_rw]
cfg.x0_nominal = [0.04; 0; 0.03; 0; 0];   % small push

%% ── Disturbance profile ──────────────────────────────────────
cfg.disturbance_enable = true;
cfg.dist_t_kick   = [];                 % single-kick override for callers
cfg.dist_kicks    = [3.0, 6.0, 10.0];  % [s] deterministic reference kicks
cfg.dist_amp      = 0.34;              % [N·m] torque kick magnitude
cfg.dist_duration = 0.05;              % [s] kick duration
cfg.dist_sine_start = 13.0;            % [s] sustained roll-axis loading start
cfg.dist_sine_end   = 20.0;            % [s] sustained roll-axis loading end
cfg.dist_sine_amp   = 0.04;            % [N·m] sinusoidal disturbance amplitude
cfg.dist_sine_freq  = 0.8;             % [Hz] sinusoidal disturbance frequency

%% ── LQR weights ──────────────────────────────────────────────
% State:   [pitch, dpitch, roll, droll, omega_rw]
% Control: [tau_rw, tau_base]
cfg.Q_lqr = diag([260, 6, 174, 5, 0.0035]);
cfg.R_lqr = diag([0.14, 0.32]);

%% ── Kalman filter noise matrices ─────────────────────────────
cfg.Q_kf = diag([1e-4, 1e-3, 1e-4, 1e-3, 1e-5]);  % process noise
cfg.R_kf = diag([cfg.noise_pitch_std^2, ...
                  cfg.noise_roll_std^2, ...
                  cfg.noise_dpitch_std^2, ...
                  cfg.noise_droll_std^2, ...
                  cfg.noise_omega_rw_std^2]);

%% ── DOB (Disturbance Observer) settings ──────────────────────
cfg.dob_enable     = true;
cfg.dob_cutoff_hz  = 5.0;   % [Hz]
cfg.dob_alpha      = 1 - exp(-2*pi*cfg.dob_cutoff_hz*cfg.dt);
cfg.integrator_leak = 0.10; % [1/s] light anti-windup leak for hybrid integrators
cfg.enable_online_id = false;
cfg.rls_lambda = 0.98;
cfg.online_id_min_updates = 60;
cfg.online_id_recompute_every = 25;

%% ── MPC settings (baseline_mpc family) ──────────────────────
cfg.mpc_horizon = 20;       % prediction horizon (steps)
cfg.Q_mpc = diag([800, 30, 600, 20, 0.002]);
cfg.R_mpc = diag([0.08, 0.2]);

%% ── H-inf-like robust weights ────────────────────────────────
cfg.Q_hinf = diag([300, 20, 220, 14, 0.0035]);
cfg.R_hinf = diag([0.16, 0.4]);
cfg.gamma_hinf = 2.0;       % robustness level

%% ── Crash / safety thresholds ────────────────────────────────
cfg.pitch_crash_rad = deg2rad(35);
cfg.roll_crash_rad  = deg2rad(35);

%% ── Output directory ─────────────────────────────────────────
cfg.results_dir = './results';
cfg.export_telemetry = true;
if ~exist(cfg.results_dir, 'dir'), mkdir(cfg.results_dir); end

end
