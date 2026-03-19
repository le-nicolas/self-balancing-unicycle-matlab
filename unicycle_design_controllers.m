function ctrl = unicycle_design_controllers(A, B, C, cfg)
%UNICYCLE_DESIGN_CONTROLLERS  Design all controller families.
%
%  Returns struct with fields:
%    .lqr_current      — delta-u LQR (matches "current" family)
%    .lqr_dob          — LQR + disturbance observer ("current_dob")
%    .hybrid_modern    — LQR with augmented integrator
%    .paper_split      — decoupled pitch/roll LQR ("paper_split_baseline")
%    .baseline_mpc     — unconstrained LQR-MPC baseline
%    .robust_hinf_like — H-inf-weighted LQR ("baseline_robust_hinf_like")

fprintf('[Design] Computing controllers...\n');

[Ad, Bd] = unicycle_discretize(A, B, cfg.dt);
Qd = cfg.Q_kf * cfg.dt;
Rd = cfg.R_kf;
[~, ~, L_kf] = dare(Ad', C', Qd, Rd);
L_kf = L_kf';

%% ── (1) CURRENT — delta-u LQR ───────────────────────────────
K_lqr = dlqr(Ad, Bd, cfg.Q_lqr, cfg.R_lqr);

ctrl.lqr_current.K   = K_lqr;
ctrl.lqr_current.L   = L_kf;
ctrl.lqr_current.dob = false;
ctrl.lqr_current.integrator = false;
ctrl.lqr_current.name = 'LQR Current';
ctrl.lqr_current.key  = 'lqr_current';

%% ── (2) CURRENT_DOB — LQR + disturbance observer ─────────────
ctrl.lqr_dob.K   = K_lqr;
ctrl.lqr_dob.L   = L_kf;
ctrl.lqr_dob.dob = true;
ctrl.lqr_dob.dob_alpha = cfg.dob_alpha;
ctrl.lqr_dob.integrator = false;
ctrl.lqr_dob.name = 'LQR + DOB';
ctrl.lqr_dob.key  = 'lqr_dob';

%% ── (3) HYBRID MODERN — LQR with integrator (augmented) ──────
% Augment state with integral of pitch and roll errors
% x_aug = [x; ∫θ_p; ∫θ_r]   (7 states)
A_pitch_i = [A(1:2,1:2), zeros(2,1);
             1, 0, -cfg.integrator_leak];
B_pitch_i = [B(1:2,2);
             0];
Q_pitch_i = diag([cfg.Q_lqr(1,1), cfg.Q_lqr(2,2), 40]);
[Ad_pitch_i, Bd_pitch_i] = unicycle_discretize(A_pitch_i, B_pitch_i, cfg.dt);
K_pitch_i = dlqr(Ad_pitch_i, Bd_pitch_i, Q_pitch_i, cfg.R_lqr(2,2));

A_roll_i = [A(3:5,3:5), zeros(3,1);
            1, 0, 0, -cfg.integrator_leak];
B_roll_i = [B(3:5,1);
            0];
Q_roll_i = diag([cfg.Q_lqr(3,3), cfg.Q_lqr(4,4), cfg.Q_lqr(5,5), 30]);
 [Ad_roll_i, Bd_roll_i] = unicycle_discretize(A_roll_i, B_roll_i, cfg.dt);
K_roll_i = dlqr(Ad_roll_i, Bd_roll_i, Q_roll_i, cfg.R_lqr(1,1));

A_aug = zeros(7);
A_aug(1:5,1:5) = A;
A_aug(6,1) = 1;
A_aug(6,6) = -cfg.integrator_leak;
A_aug(7,3) = 1;
A_aug(7,7) = -cfg.integrator_leak;
B_aug = [B; zeros(2,2)];
Q_aug = diag([cfg.Q_lqr(1,1), cfg.Q_lqr(2,2), cfg.Q_lqr(3,3), ...
    cfg.Q_lqr(4,4), cfg.Q_lqr(5,5), 40, 30]);
R_aug = cfg.R_lqr;
K_aug = zeros(2,7);
K_aug(1,[3,4,5,7]) = K_roll_i;
K_aug(2,[1,2,6])   = K_pitch_i;

ctrl.hybrid_modern.K     = K_aug;
ctrl.hybrid_modern.L     = L_kf;
ctrl.hybrid_modern.A_aug = A_aug;
ctrl.hybrid_modern.B_aug = B_aug;
ctrl.hybrid_modern.Q_aug = Q_aug;
ctrl.hybrid_modern.R_aug = R_aug;
ctrl.hybrid_modern.dob   = false;
ctrl.hybrid_modern.integrator = true;
ctrl.hybrid_modern.name  = 'Hybrid Modern (LQR+I)';
ctrl.hybrid_modern.key   = 'hybrid_modern';

%% ── (4) PAPER SPLIT BASELINE — decoupled pitch / roll ────────
% Separate LQR for each axis (mirrors the 2013 Lee et al. paper)
% Pitch subsystem:  [θ_p, θ̇_p]
A_p = A(1:2, 1:2);
B_p = B(1:2, 2);            % only τ_base drives pitch
Q_p = diag([80, 8]);
R_p = 1.0;
[Ad_p, Bd_p] = unicycle_discretize(A_p, B_p, cfg.dt);
K_p = dlqr(Ad_p, Bd_p, Q_p, R_p);

% Roll subsystem: [θ_r, θ̇_r, ω_rw]
A_r = A(3:5, 3:5);
B_r = B(3:5, 1);            % only τ_rw drives roll
Q_r = diag([60, 6, 0.01]);
R_r = 0.5;
[Ad_r, Bd_r] = unicycle_discretize(A_r, B_r, cfg.dt);
K_r = dlqr(Ad_r, Bd_r, Q_r, R_r);

ctrl.paper_split.K_pitch = K_p;
ctrl.paper_split.K_roll  = K_r;
ctrl.paper_split.L       = L_kf;
ctrl.paper_split.dob     = false;
ctrl.paper_split.integrator = false;
ctrl.paper_split.name    = 'Paper Split Baseline';
ctrl.paper_split.key     = 'paper_split';

%% ── (5) BASELINE MPC — unconstrained LQR-MPC ────────────────
% For unconstrained MPC the optimal gain equals LQR
% (we keep separate tuning weights as in the original)
K_mpc = dlqr(Ad, Bd, cfg.Q_mpc, cfg.R_mpc);
cfg.mpc_constrained = license('test', 'optimization_toolbox');
quadprog_available = exist('quadprog', 'file') > 0;
mpc_qp_active = cfg.mpc_constrained && quadprog_available;
if mpc_qp_active
    fprintf('[MPC] Optimization Toolbox available — constrained QP active\n');
elseif cfg.mpc_constrained && ~quadprog_available
    fprintf('[MPC] Optimization Toolbox license available but quadprog missing — falling back to dlqr\n');
else
    fprintf('[MPC] Optimization Toolbox not available — falling back to dlqr\n');
end

[P_mpc, ~, ~] = dare(Ad, Bd, cfg.Q_mpc, cfg.R_mpc);
[Psi_mpc, Gamma_mpc] = compute_mpc_prediction(Ad, Bd, cfg.mpc_horizon);
Q_bar = build_mpc_state_weights(cfg.Q_mpc, P_mpc, cfg.mpc_horizon);
R_bar = kron(eye(cfg.mpc_horizon), cfg.R_mpc);
H_mpc = 2 * (Gamma_mpc' * Q_bar * Gamma_mpc + R_bar);
H_mpc = (H_mpc + H_mpc') / 2 + 1e-8 * eye(size(H_mpc));
F_half = 2 * (Gamma_mpc' * Q_bar * Psi_mpc);
lb = repmat([-cfg.tau_rw_max; -cfg.tau_base_max], cfg.mpc_horizon, 1);
ub = repmat([ cfg.tau_rw_max;  cfg.tau_base_max], cfg.mpc_horizon, 1);

ctrl.baseline_mpc.K   = K_mpc;
ctrl.baseline_mpc.K_fallback = K_mpc;
ctrl.baseline_mpc.L   = L_kf;
ctrl.baseline_mpc.N   = cfg.mpc_horizon;
ctrl.baseline_mpc.N_horiz = cfg.mpc_horizon;
ctrl.baseline_mpc.dob = false;
ctrl.baseline_mpc.integrator = false;
ctrl.baseline_mpc.constrained = mpc_qp_active;
ctrl.baseline_mpc.license_gate = cfg.mpc_constrained;
ctrl.baseline_mpc.quadprog_available = quadprog_available;
ctrl.baseline_mpc.name = 'Baseline MPC (unconstrained)';
ctrl.baseline_mpc.key  = 'baseline_mpc';
ctrl.baseline_mpc.Psi = Psi_mpc;
ctrl.baseline_mpc.Phi = Gamma_mpc;
ctrl.baseline_mpc.H = H_mpc;
ctrl.baseline_mpc.F_half = F_half;
ctrl.baseline_mpc.lb = lb;
ctrl.baseline_mpc.ub = ub;
ctrl.baseline_mpc.R_bar = R_bar;
ctrl.baseline_mpc.R_bar_block = cfg.R_mpc;
ctrl.baseline_mpc.P_terminal = P_mpc;
ctrl.baseline_mpc.rho_w = 0.1;
if mpc_qp_active
    ctrl.baseline_mpc.qp_options = optimoptions('quadprog', ...
        'Display', 'off', 'MaxIterations', 200);
else
    ctrl.baseline_mpc.qp_options = [];
end

%% ── (6) ROBUST H-INF-LIKE — heavier penalty LQR ─────────────
K_hinf = dlqr(Ad, Bd, cfg.Q_hinf, cfg.R_hinf);

ctrl.robust_hinf_like.K   = K_hinf;
ctrl.robust_hinf_like.L   = L_kf;
ctrl.robust_hinf_like.dob = false;
ctrl.robust_hinf_like.integrator = false;
ctrl.robust_hinf_like.name = 'Robust H-inf-like';
ctrl.robust_hinf_like.key  = 'robust_hinf_like';

%% ── Summary ──────────────────────────────────────────────────
fprintf('[Design] Controllers designed:\n');
fns = fieldnames(ctrl);
for i = 1:numel(fns)
    fprintf('         • %-25s\n', ctrl.(fns{i}).name);
end
end


function [Psi, Gamma] = compute_mpc_prediction(Ad, Bd, N)
%COMPUTE_MPC_PREDICTION  Build lifted prediction matrices X = Psi*x0 + Gamma*U.
nx = size(Ad, 1);
nu = size(Bd, 2);

Psi = zeros(N * nx, nx);
Gamma = zeros(N * nx, N * nu);
for i = 1:N
    row_idx = (i-1) * nx + (1:nx);
    Psi(row_idx, :) = Ad^i;
    for j = 1:i
        col_idx = (j-1) * nu + (1:nu);
        Gamma(row_idx, col_idx) = Ad^(i-j) * Bd;
    end
end
end


function Q_bar = build_mpc_state_weights(Q, P, N)
%BUILD_MPC_STATE_WEIGHTS  Assemble blkdiag(Q,...,Q,P) over the prediction horizon.
nx = size(Q, 1);
Q_bar = zeros(N * nx);
for i = 1:N
    row_idx = (i-1) * nx + (1:nx);
    if i < N
        Q_bar(row_idx, row_idx) = Q;
    else
        Q_bar(row_idx, row_idx) = P;
    end
end
end
