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
cfg_rls = cfg;
cfg_rls.enable_online_id = true;
cfg_rls.m_body = cfg.m_body * 1.20;
cfg_rls.online_id_min_updates = 20;
cfg_rls.online_id_recompute_every = 10;
res_rls = unicycle_simulate(A, B, C, ctrl_all.hybrid_modern, cfg_rls);
theta_shift = norm(res_rls.rls_theta_final - res_rls.rls_theta_init);
pass14 = res_rls.survival && res_rls.online_id_gain_updates > 0 && theta_shift > 1e-6;
[pass_count, fail_count] = check(pass14, ...
    sprintf('TEST 14: Hybrid online ID adapts and survives (updates=%d, theta shift=%.3e)', ...
        res_rls.online_id_gain_updates, theta_shift), ...
    pass_count, fail_count);

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
