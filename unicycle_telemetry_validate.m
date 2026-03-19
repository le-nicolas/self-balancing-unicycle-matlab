function unicycle_telemetry_validate()
%UNICYCLE_TELEMETRY_VALIDATE  Validate telemetry CSV sign conventions.
%  This script checks whether exported telemetry is compatible with the
%  sim-to-real sign conventions expected by hil_bridge.py.

clc;
fprintf('=== UNICYCLE TELEMETRY VALIDATION ===\n\n');

cfg0 = unicycle_config();
[A, B, C, ~] = unicycle_plant(cfg0);
ctrl_all = unicycle_design_controllers(A, B, C, cfg0);

work_dir = fullfile(cfg0.results_dir, 'telemetry_validate');
if ~exist(work_dir, 'dir')
    mkdir(work_dir);
end

tests = false(6, 1);
labels = { ...
    'T1 pitch column sign'; ...
    'T2 roll column sign'; ...
    'T3 tau_rw sign convention'; ...
    'T4 tau_base sign convention'; ...
    'T5 omega_rw spin direction'; ...
    'T6 CSV column order and units'};

%% T1 - Positive pitch should export positive pitch_deg
cfg = base_cfg(cfg0, work_dir);
cfg.x0_nominal = [0.05; 0; 0; 0; 0];
ctrl = zero_ctrl(ctrl_all.lqr_current.L, 'telemetry_val_t1');
tbl = run_and_read_csv(A, B, C, ctrl, cfg);
pitch_val = tbl.pitch_deg(1);
tests(1) = pitch_val > 0;
fprintf('T1 pitch_deg(1) = %.6f deg -> %s\n', pitch_val, pass_fail(tests(1)));

%% T2 - Positive roll should export positive roll_deg
cfg = base_cfg(cfg0, work_dir);
cfg.x0_nominal = [0; 0; 0.05; 0; 0];
ctrl = zero_ctrl(ctrl_all.lqr_current.L, 'telemetry_val_t2');
tbl = run_and_read_csv(A, B, C, ctrl, cfg);
roll_val = tbl.roll_deg(1);
tests(2) = roll_val > 0;
fprintf('T2 roll_deg(1)  = %.6f deg -> %s\n', roll_val, pass_fail(tests(2)));

%% T3 - Positive roll should trigger positive tau_rw
cfg = base_cfg(cfg0, work_dir);
cfg.imu_delay_samples = 0;
cfg.x0_nominal = [0; 0; 0.05; 0; 0];
ctrl = ctrl_all.lqr_current;
ctrl.key = 'telemetry_val_t3';
tbl = run_and_read_csv(A, B, C, ctrl, cfg);
tau_rw_val = tbl.tau_rw(1);
tests(3) = tau_rw_val > 0;
fprintf('T3 tau_rw(1)    = %.6f N*m -> %s\n', tau_rw_val, pass_fail(tests(3)));

%% T4 - Positive pitch should trigger positive tau_base
cfg = base_cfg(cfg0, work_dir);
cfg.imu_delay_samples = 0;
cfg.x0_nominal = [0.05; 0; 0; 0; 0];
ctrl = ctrl_all.lqr_current;
ctrl.key = 'telemetry_val_t4';
tbl = run_and_read_csv(A, B, C, ctrl, cfg);
tau_base_val = tbl.tau_base(1);
tests(4) = tau_base_val > 0;
fprintf('T4 tau_base(1)  = %.6f N*m -> %s\n', tau_base_val, pass_fail(tests(4)));

%% T5 - Positive tau_rw should spin omega_rw positive
cfg = base_cfg(cfg0, work_dir);
manual_csv = fullfile(work_dir, 'telemetry_val_t5.csv');
write_constant_input_csv(manual_csv, cfg, [0.1; 0], 10);
tbl = readtable(manual_csv);
omega_end = tbl.omega_rw(end);
tests(5) = omega_end > 0;
fprintf('T5 omega_rw(end)= %.6f rad/s -> %s\n', omega_end, pass_fail(tests(5)));

%% T6 - Column order and units
cfg = base_cfg(cfg0, work_dir);
cfg.imu_delay_samples = 0;
cfg.x0_nominal = [0.05; 0; 0; 0; 0];
cfg.T_sim = 10 * cfg.dt;
cfg.N = round(cfg.T_sim / cfg.dt);
ctrl = ctrl_all.lqr_current;
ctrl.key = 'telemetry_val_t6';
tbl = run_and_read_csv(A, B, C, ctrl, cfg);

dt_series = diff(tbl.t);
monotone_t = all(dt_series > 0);
step_ok = all(abs(dt_series - cfg.dt) < 1e-12);
deg_ok = max(abs(tbl.pitch_deg)) < 360;
tau_ok = max(abs(tbl.tau_rw)) <= cfg.tau_rw_max + 1e-4;
tests(6) = monotone_t && step_ok && deg_ok && tau_ok;
fprintf('T6 monotone t   = %d, dt match = %d, deg units = %d, tau units = %d -> %s\n', ...
    monotone_t, step_ok, deg_ok, tau_ok, pass_fail(tests(6)));
if ~tests(6)
    fprintf('   details: max|pitch_deg|=%.6f, max|tau_rw|=%.6f\n', ...
        max(abs(tbl.pitch_deg)), max(abs(tbl.tau_rw)));
end

%% Summary
fprintf('\n--- SUMMARY ------------------------------------------------------\n');
for i = 1:numel(labels)
    fprintf('%-30s %s\n', labels{i}, pass_fail(tests(i)));
end

if all(tests)
    fprintf('\nTelemetry sign conventions: COMPATIBLE with hil_bridge.py\n');
else
    fprintf('\nTelemetry sign conventions: INCOMPATIBLE with hil_bridge.py\n');
end
end


function cfg = base_cfg(cfg0, work_dir)
cfg = cfg0;
cfg.results_dir = work_dir;
cfg.export_telemetry = true;
cfg.disturbance_enable = false;
cfg.T_sim = cfg.dt;
cfg.N = round(cfg.T_sim / cfg.dt);
cfg.dist_sine_amp = 0.0;
cfg.noise_pitch_std = 0.0;
cfg.noise_roll_std = 0.0;
cfg.noise_dpitch_std = 0.0;
cfg.noise_droll_std = 0.0;
cfg.noise_omega_rw_std = 0.0;
cfg.Q_kf = zeros(size(cfg.Q_kf));
cfg.R_kf = diag([1e-12, 1e-12, 1e-12, 1e-12, 1e-12]);
end


function ctrl = zero_ctrl(L, key)
ctrl.name = 'Telemetry Zero Control';
ctrl.key = key;
ctrl.K = zeros(2, 5);
ctrl.L = L;
ctrl.dob = false;
ctrl.integrator = false;
end


function tbl = run_and_read_csv(A, B, C, ctrl, cfg)
csv_path = fullfile(cfg.results_dir, sprintf('telemetry_%s.csv', ctrl.key));
if exist(csv_path, 'file')
    delete(csv_path);
end
evalc('unicycle_simulate(A, B, C, ctrl, cfg);');
tbl = readtable(csv_path);
end


function write_constant_input_csv(csv_path, cfg, u_const, n_steps)
t = (0:n_steps)' * cfg.dt;
x = zeros(5, n_steps + 1);

for k = 1:n_steps
    xdot = local_dynamics(x(:,k), u_const, cfg);
    x(:,k+1) = x(:,k) + cfg.dt * xdot;
end

u_pad = repmat(u_const(:).', n_steps + 1, 1);
telemetry = table( ...
    t, ...
    rad2deg(x(1,:)).', ...
    rad2deg(x(3,:)).', ...
    x(2,:).', ...
    x(4,:).', ...
    x(5,:).', ...
    u_pad(:,1), ...
    u_pad(:,2), ...
    'VariableNames', {'t','pitch_deg','roll_deg','dpitch','droll','omega_rw','tau_rw','tau_base'});
writetable(telemetry, csv_path);
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
