%% ============================================================
%  UNICYCLE VERIFICATION TESTS
%  Mirrors pytest CI checks in the original MuJoCo project.
%
%  Run: unicycle_verify
%  All tests must print PASS. If any print FAIL, inspect output.
% ============================================================
function unicycle_verify()
clc;
fprintf('=== UNICYCLE MATLAB VERIFICATION TESTS ===\n\n');
pass_count = 0;
fail_count = 0;

cfg = unicycle_config();
[A, B, C, ~] = unicycle_plant(cfg);

%% TEST 1: Controllability
rk = rank(ctrb(A,B));
[pass_count, fail_count] = check(rk == size(A,1), ...
    'TEST 1: Controllability full rank', pass_count, fail_count);

%% TEST 2: Observability
rk2 = rank(obsv(A,C));
[pass_count, fail_count] = check(rk2 == size(A,1), ...
    'TEST 2: Observability full rank', pass_count, fail_count);

%% TEST 3: Open-loop unstable (2 RHP poles expected)
ev = eig(A);
n_rhp = sum(real(ev) > 0);
[pass_count, fail_count] = check(n_rhp == 2, ...
    'TEST 3: Plant has exactly 2 unstable eigenvalues', pass_count, fail_count);

%% TEST 4: LQR stabilises the plant
ctrl_all = unicycle_design_controllers(A, B, C, cfg);
K_lqr = ctrl_all.lqr_current.K;
ev_cl = eig(A - B*K_lqr);
[pass_count, fail_count] = check(all(real(ev_cl) < 0), ...
    'TEST 4: LQR closed-loop is stable', pass_count, fail_count);

%% TEST 5: Paper-split controllers each stabilise their subsystem
K_p = ctrl_all.paper_split.K_pitch;
K_r = ctrl_all.paper_split.K_roll;
A_p = A(1:2,1:2); B_p = B(1:2,2);
A_r = A(3:5,3:5); B_r = B(3:5,1);
ev_p = eig(A_p - B_p*K_p);
ev_r = eig(A_r - B_r*K_r);
[pass_count, fail_count] = check(all(real(ev_p) < 0) && all(real(ev_r) < 0), ...
    'TEST 5: Paper-split subsystem poles stable', pass_count, fail_count);

%% TEST 6: Zero-input, zero-state → stays at origin (linearised)
xk = zeros(5,1);
uk = zeros(2,1);
[Ad, Bd] = unicycle_discretize(A, B, cfg.dt);
for k = 1:100
    xk = Ad*xk + Bd*uk;
end
[pass_count, fail_count] = check(norm(xk) < 1e-10, ...
    'TEST 6: Zero equilibrium is stable (lin)', pass_count, fail_count);

%% TEST 7: Small IC converges with LQR (closed-loop linear simulation)
xk = [0.05; 0; 0.03; 0; 0];
for k = 1:1000
    uk = -K_lqr * xk;
    uk = max(min(uk, [cfg.tau_rw_max; cfg.tau_base_max]), ...
                     [-cfg.tau_rw_max; -cfg.tau_base_max]);
    xk = Ad*xk + Bd*uk;
end
[pass_count, fail_count] = check(norm(xk([1,3])) < deg2rad(1), ...
    'TEST 7: LQR drives small IC to near-zero', pass_count, fail_count);

%% TEST 8: Actuator saturation is respected
res_lqr = unicycle_simulate(A, B, C, ctrl_all.lqr_current, cfg);
u_max = max(abs(res_lqr.u), [], 2);
[pass_count, fail_count] = check(u_max(1) <= cfg.tau_rw_max + 1e-9 && ...
                                  u_max(2) <= cfg.tau_base_max + 1e-9, ...
    'TEST 8: Actuator saturation respected', pass_count, fail_count);

%% TEST 9: Kalman estimator tracks truth (RMS error < 2 deg)
Nv = res_lqr.crash_step;
pitch_est_err = rad2deg(rms(res_lqr.x_true(1,1:Nv) - res_lqr.x_est(1,1:Nv)));
[pass_count, fail_count] = check(pitch_est_err < 2.0, ...
    sprintf('TEST 9: Kalman pitch error < 2 deg (got %.3f)', pitch_est_err), ...
    pass_count, fail_count);

%% TEST 10: At least 3 controllers survive full simulation
n_survived = 0;
fns = fieldnames(ctrl_all);
for i = 1:numel(fns)
    r = unicycle_simulate(A, B, C, ctrl_all.(fns{i}), cfg);
    if r.survival, n_survived = n_survived + 1; end
end
[pass_count, fail_count] = check(n_survived >= 3, ...
    sprintf('TEST 10: ≥3 controllers survive (got %d)', n_survived), ...
    pass_count, fail_count);

%% TEST 11: MPC constraint satisfaction (toolbox-conditional)
mpc_active = isfield(ctrl_all.baseline_mpc, 'constrained') && ctrl_all.baseline_mpc.constrained;
if mpc_active
    res_mpc = unicycle_simulate(A, B, C, ctrl_all.baseline_mpc, cfg);
    u_max_mpc = max(abs(res_mpc.u), [], 2);
    pass11 = u_max_mpc(1) <= cfg.tau_rw_max + 1e-4 && ...
             u_max_mpc(2) <= cfg.tau_base_max + 1e-4;
    [pass_count, fail_count] = check(pass11, ...
        sprintf('TEST 11: MPC u within bounds (max rw=%.4f, base=%.4f)', ...
            u_max_mpc(1), u_max_mpc(2)), pass_count, fail_count);
else
    [pass_count, fail_count] = check(true, ...
        'TEST 11: MPC constraint check skipped (no toolbox)', ...
        pass_count, fail_count);
end

%% TEST 12: Constrained MPC respects tightened bounds
if mpc_active
    cfg_tight = cfg;
    cfg_tight.tau_rw_max = 0.10;
    cfg_tight.tau_base_max = 0.15;
    ctrl_tight = unicycle_design_controllers(A, B, C, cfg_tight);
    if isfield(ctrl_tight.baseline_mpc, 'constrained') && ctrl_tight.baseline_mpc.constrained
        res_tight = unicycle_simulate(A, B, C, ctrl_tight.baseline_mpc, cfg_tight);
        u_tight_max = max(abs(res_tight.u(1,:)));
        pass12 = u_tight_max <= cfg_tight.tau_rw_max + 1e-4;
        [pass_count, fail_count] = check(pass12, ...
            sprintf('TEST 12: MPC respects tightened bounds (max=%.4f, lim=%.4f)', ...
                u_tight_max, cfg_tight.tau_rw_max), pass_count, fail_count);
    else
        [pass_count, fail_count] = check(true, ...
            'TEST 12: MPC bounds test skipped (no toolbox)', ...
            pass_count, fail_count);
    end
else
    [pass_count, fail_count] = check(true, ...
        'TEST 12: MPC bounds test skipped (no toolbox)', ...
        pass_count, fail_count);
end

%% TEST 13: Transient metrics are non-zero and MPC outranks paper_split
res_ps = unicycle_simulate(A, B, C, ctrl_all.paper_split, cfg);
res_mpc = unicycle_simulate(A, B, C, ctrl_all.baseline_mpc, cfg);
[mpk_ps, mst_ps] = compute_transient_metrics_verify(res_ps, cfg);
[mpk_mpc, mst_mpc] = compute_transient_metrics_verify(res_mpc, cfg);
score_ps = compute_benchmark_score_verify(res_ps, mpk_ps, mst_ps);
score_mpc = compute_benchmark_score_verify(res_mpc, mpk_mpc, mst_mpc);
pass13 = score_mpc > score_ps && mst_ps > 0 && mpk_ps > 0;
[pass_count, fail_count] = check(pass13, ...
    'TEST 13: Disturbance-axis transient metrics non-zero and MPC > paper_split', ...
    pass_count, fail_count);

%% TEST 14: Hybrid modern survives with online ID under mass perturbation
cfg_t14 = cfg;
cfg_t14.enable_online_id = true;
cfg_t14.I_body = cfg.I_body * 1.20;
cfg_t14.online_id_min_updates = 40;
cfg_t14.online_id_recompute_every = 20;
[A_t14, B_t14, C_t14, ~] = unicycle_plant(cfg_t14);
ctrl_t14 = ctrl_all.hybrid_modern;
res_t14 = unicycle_simulate(A_t14, B_t14, C_t14, ctrl_t14, cfg_t14);
pass14 = res_t14.survival && res_t14.pitch_rms < deg2rad(3.0);
[pass_count, fail_count] = check(pass14, ...
    sprintf('TEST 14: hybrid online ID survives +20%% pitch inertia  (survival=%d, pitchRMS=%.4f rad)', ...
        res_t14.survival, res_t14.pitch_rms), ...
    pass_count, fail_count);

%% TEST 15: Firmware export produces valid JSON and C header
out_dir = tempname();
mkdir(out_dir);
unicycle_export_firmware(ctrl_all, cfg, out_dir);
json_file = fullfile(out_dir, 'firmware_params.json');
h_file = fullfile(out_dir, 'firmware_params.h');

json_info = dir(json_file);
h_info = dir(h_file);
json_ok = exist(json_file, 'file') == 2 && ~isempty(json_info) && json_info.bytes > 100;
h_ok = exist(h_file, 'file') == 2 && ~isempty(h_info) && h_info.bytes > 100;
try
    raw = fileread(json_file);
    data = jsondecode(raw);
    has_physical = isfield(data, 'physical');
    has_kalman = isfield(data, 'kalman');
    has_controllers = isfield(data, 'controllers');
    has_lqr = has_controllers && isfield(data.controllers, 'lqr_current');
    json_valid = has_physical && has_kalman && has_controllers && has_lqr;
catch
    json_valid = false;
end
h_text = fileread(h_file);
h_has_dt = contains(h_text, 'CTRL_DT');
h_has_k = contains(h_text, 'K_LQR_CURRENT');
pass15 = json_ok && h_ok && json_valid && h_has_dt && h_has_k;
[pass_count, fail_count] = check(pass15, ...
    'TEST 15: firmware export produces valid JSON and C header', ...
    pass_count, fail_count);
rmdir(out_dir, 's');

%% TEST 16: Firmware parity within 1e-4 N*m
K_fw = single(ctrl_all.lqr_current.K);
x0_s = single(cfg.x0_nominal);
Ad_s = single(Ad);
Bd_s = single(Bd);

N_par = round(2.0 / cfg.dt);
u_ml = zeros(2, N_par);
u_fw = zeros(2, N_par);
x_ml = cfg.x0_nominal;
x_fw = x0_s;
tau_rw_s = single(cfg.tau_rw_max);
tau_base_s = single(cfg.tau_base_max);

for k = 1:N_par
    u_k_ml = -ctrl_all.lqr_current.K * x_ml;
    u_k_ml(1) = max(min(u_k_ml(1), cfg.tau_rw_max), -cfg.tau_rw_max);
    u_k_ml(2) = max(min(u_k_ml(2), cfg.tau_base_max), -cfg.tau_base_max);
    u_ml(:,k) = u_k_ml;
    x_ml = Ad * x_ml + Bd * u_k_ml;

    u_k_fw = -K_fw * x_fw;
    u_k_fw(1) = max(min(u_k_fw(1), tau_rw_s), -tau_rw_s);
    u_k_fw(2) = max(min(u_k_fw(2), tau_base_s), -tau_base_s);
    u_fw(:,k) = double(u_k_fw);
    x_fw = Ad_s * x_fw + Bd_s * u_k_fw;
end

max_rw_err = max(abs(u_ml(1,:) - u_fw(1,:)));
max_base_err = max(abs(u_ml(2,:) - u_fw(2,:)));
pass16 = max_rw_err < 1e-4 && max_base_err < 1e-4;
[pass_count, fail_count] = check(pass16, ...
    sprintf('TEST 16: firmware parity within 1e-4 N*m  (rw=%.2e, base=%.2e)', ...
        max_rw_err, max_base_err), pass_count, fail_count);

%% TEST 17: HIL smoke test 7/7 scenarios pass
hil_results = unicycle_hil_smoke(ctrl_all, cfg);
pass17 = hil_results.passed;
[pass_count, fail_count] = check(pass17, ...
    'TEST 17: HIL smoke test 7/7 scenarios pass', ...
    pass_count, fail_count);
if ~pass17
    for si = 1:numel(hil_results.scenarios)
        sc = hil_results.scenarios(si);
        if ~sc.passed
            fprintf('    FAIL in scenario %d (%s): %s\n', ...
                si, sc.name, sc.fail_reason);
        end
    end
end

%% ── Summary ──────────────────────────────────────────────────
fprintf('\n====================================\n');
fprintf('PASS: %d / %d\n', pass_count, pass_count + fail_count);
fprintf('FAIL: %d / %d\n', fail_count, pass_count + fail_count);
if fail_count == 0
    fprintf('✓  ALL TESTS PASSED\n');
else
    fprintf('✗  SOME TESTS FAILED — see above\n');
end
fprintf('====================================\n');
end


function [p, f] = check(cond, label, p, f)
if cond
    fprintf('  ✓ PASS  %s\n', label);
    p = p + 1;
else
    fprintf('  ✗ FAIL  %s\n', label);
    f = f + 1;
end
end


function score = compute_benchmark_score_verify(res, mpk, mst)
score = 100 * res.survival_frac ...
      - 50 * res.pitch_rms ...
      - 40 * res.roll_rms ...
      -  5 * res.ctrl_rms ...
      - 20 * mpk ...
      - 15 * mst;
end


function [mpk, mst] = compute_transient_metrics_verify(res, cfg)
t = res.t;
x_pitch = res.x_true(1, :);
x_roll = res.x_true(3, :);
dt = cfg.dt;
thresh = 0.02;
win_s = 1.5;
consec = 10;

if ~isfield(cfg, 'dist_kicks') || isempty(cfg.dist_kicks)
    mpk = 0;
    mst = 0;
    return;
end

kicks = cfg.dist_kicks;
peaks = zeros(1, numel(kicks));
fracs = zeros(1, numel(kicks));

for ki = 1:numel(kicks)
    k_start = round(kicks(ki) / dt) + 1;
    k_end = min(round((kicks(ki) + win_s) / dt) + 1, numel(t));
    if k_start > numel(t)
        peaks(ki) = 0;
        fracs(ki) = 0;
        continue;
    end

    window_pitch = x_pitch(k_start:k_end);
    window_roll = x_roll(k_start:k_end);
    if max(abs(window_roll)) >= max(abs(window_pitch))
        window_angle = window_roll;
    else
        window_angle = window_pitch;
    end
    peaks(ki) = max(abs(window_angle));

    settled = false;
    if numel(window_angle) >= consec
        [~, peak_idx] = max(abs(window_angle));
        last_start = numel(window_angle) - consec + 1;
        for ks = peak_idx:last_start
            if all(abs(window_angle(ks:ks+consec-1)) < thresh)
                fracs(ki) = min(((ks - 1) * dt) / win_s, 1.0);
                settled = true;
                break;
            end
        end
    end
    if ~settled
        fracs(ki) = 1.0;
    end
end

mpk = mean(peaks);
mst = mean(fracs);
end
