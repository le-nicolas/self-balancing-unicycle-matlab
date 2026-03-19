function results = unicycle_test_firmware_parity(out_dir)
%UNICYCLE_TEST_FIRMWARE_PARITY  Verify exported firmware params match MATLAB control.
if nargin < 1 || isempty(out_dir)
    out_dir = tempname();
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

cfg = unicycle_config();
[A, B, C, ~] = unicycle_plant(cfg);
ctrl_all = unicycle_design_controllers(A, B, C, cfg);

cfg_par = cfg;
cfg_par.disturbance_enable = false;
cfg_par.export_telemetry = false;
[Ad, Bd] = unicycle_discretize(A, B, cfg_par.dt);

unicycle_export_firmware(ctrl_all, cfg_par, out_dir);
json_file = fullfile(out_dir, 'firmware_params.json');
raw = fileread(json_file);
data = jsondecode(raw);

K_loaded = reshape_row_major(data.controllers.lqr_current.K, ...
    data.controllers.lqr_current.K_rows, data.controllers.lqr_current.K_cols);
Ad_loaded = reshape_row_major(data.kalman.Ad, data.kalman.Ad_rows, data.kalman.Ad_cols);
Bd_loaded = reshape_row_major(data.kalman.Bd, data.kalman.Bd_rows, data.kalman.Bd_cols);

N_par = round(2.0 / cfg_par.dt);
u_matlab = zeros(2, N_par);
u_firmware = zeros(2, N_par);
x_ml = cfg_par.x0_nominal;
x_fw = single(cfg_par.x0_nominal);

K_fw = single(K_loaded);
Ad_fw = single(Ad_loaded);
Bd_fw = single(Bd_loaded);
tau_rw_s = single(cfg_par.tau_rw_max);
tau_base_s = single(cfg_par.tau_base_max);

for k = 1:N_par
    u_k_ml = -ctrl_all.lqr_current.K * x_ml;
    u_k_ml(1) = max(min(u_k_ml(1), cfg_par.tau_rw_max), -cfg_par.tau_rw_max);
    u_k_ml(2) = max(min(u_k_ml(2), cfg_par.tau_base_max), -cfg_par.tau_base_max);
    u_matlab(:,k) = u_k_ml;
    x_ml = Ad * x_ml + Bd * u_k_ml;

    u_k_fw = -K_fw * x_fw;
    u_k_fw(1) = max(min(u_k_fw(1), tau_rw_s), -tau_rw_s);
    u_k_fw(2) = max(min(u_k_fw(2), tau_base_s), -tau_base_s);
    u_firmware(:,k) = double(u_k_fw);
    x_fw = Ad_fw * x_fw + Bd_fw * u_k_fw;
end

rw_err = abs(u_matlab(1,:) - u_firmware(1,:));
base_err = abs(u_matlab(2,:) - u_firmware(2,:));
[max_rw_err, idx_rw] = max(rw_err);
[max_base_err, idx_base] = max(base_err);
json_roundtrip_ok = abs(K_loaded(1,1) - ctrl_all.lqr_current.K(1,1)) < 1e-10;

rw_pass = max_rw_err < 1e-4;
base_pass = max_base_err < 1e-4;

fprintf('--- FIRMWARE PARITY RESULTS ---\n');
fprintf('RW   max error: %.2e N*m  (limit 1e-4) -> %s\n', ...
    max_rw_err, pass_str(rw_pass));
fprintf('Base max error: %.2e N*m  (limit 1e-4) -> %s\n', ...
    max_base_err, pass_str(base_pass));
fprintf('JSON round-trip: %s\n', pass_str(json_roundtrip_ok));

if ~rw_pass
    fprintf('  RW mismatch at step %d: matlab=%.8f firmware=%.8f\n', ...
        idx_rw, u_matlab(1, idx_rw), u_firmware(1, idx_rw));
end
if ~base_pass
    fprintf('  Base mismatch at step %d: matlab=%.8f firmware=%.8f\n', ...
        idx_base, u_matlab(2, idx_base), u_firmware(2, idx_base));
end

results.max_rw_err = max_rw_err;
results.max_base_err = max_base_err;
results.idx_rw = idx_rw;
results.idx_base = idx_base;
results.json_roundtrip_ok = json_roundtrip_ok;
results.rw_pass = rw_pass;
results.base_pass = base_pass;
results.json_file = json_file;
results.header_file = fullfile(out_dir, 'firmware_params.h');
end


function M = reshape_row_major(flat, rows, cols)
M = reshape(flat, cols, rows).';
end


function s = pass_str(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end
