%% ============================================================
%  SELF-BALANCING UNICYCLE ROBOT — MATLAB SIMULATION
%  Converted from: le-nicolas/Self-Balancing-Unicycle-Robot
%                  le-nicolas/reaction-wheel-balancer-sim2real
%
%  Architecture:
%    - Body (stick) balanced by a reaction wheel (roll axis)
%    - Base drive wheel corrects pitch axis
%    - LQR controller with Kalman state estimation
%    - Disturbance observer (DOB) optional feed-forward
%    - Configurable controller families: LQR, MPC, H-inf-like
%
%  Run this file. Sub-functions are in separate .m files.
%  Results and plots are saved to ./results/
% ============================================================
clear; clc; close all;
addpath(genpath('.'));

%% ── 0. USER CONFIGURATION ──────────────────────────────────
cfg = unicycle_config();          % load all parameters

%% ── 1. LINEARISED PLANT ─────────────────────────────────────
[A, B, C, D] = unicycle_plant(cfg);
sys = ss(A, B, C, D);
fprintf('[Plant] Eigenvalues of open-loop A:\n');
disp(eig(A));

%% ── 2. CONTROLLER DESIGN ────────────────────────────────────
ctrl = unicycle_design_controllers(A, B, C, cfg);

%% ── 3. RUN SIMULATIONS ──────────────────────────────────────
families = fieldnames(ctrl);
results  = struct();

for fi = 1:numel(families)
    fname = families{fi};
    fprintf('\n[SIM] Running controller: %s\n', fname);
    results.(fname) = unicycle_simulate(A, B, C, ctrl.(fname), cfg);
end

%% ── 4. BENCHMARK & PLOTS ────────────────────────────────────
unicycle_benchmark(results, cfg);
unicycle_plot_all(results, cfg);

fprintf('\n[DONE] All simulations complete. Results in ./results/\n');
