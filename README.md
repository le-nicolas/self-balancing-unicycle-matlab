# Self-Balancing Unicycle MATLAB

Standalone MATLAB simulation, verification, and benchmarking suite for a reaction-wheel self-balancing unicycle robot.

[![Model Validation](https://img.shields.io/badge/model%20validation-8%2F8%20PASS-brightgreen)](#current-status)
[![Verification](https://img.shields.io/badge/verification-13%2F14%20PASS-yellow)](#current-status)
[![Telemetry](https://img.shields.io/badge/telemetry-HIL%20compatible-brightgreen)](#telemetry)
[![Constrained MPC](https://img.shields.io/badge/constrained%20MPC-implemented%20with%20fallback-blue)](#controllers)
[![Online ID](https://img.shields.io/badge/online%20ID-WIP-orange)](#known-issues)

## Current Status

| Area | Status | Notes |
|---|---|---|
| Physics and linearisation | PASS | `unicycle_model_validate` = 8/8 |
| Core verification | PASS | `unicycle_verify` = 13/14 |
| Telemetry export | PASS | CSV sign convention matches HIL expectations |
| Nominal benchmark | PASS | All six controller families survive |
| Stress benchmark | PASS | `baseline_mpc` is highest in current Monte Carlo runs |
| Constrained MPC | READY WITH FALLBACK | Active when `quadprog` is available |
| Online adaptive ID | IN PROGRESS | RLS branch adapts numerically but still destabilises `hybrid_modern` under the current perturbation test |

Status bar:

```text
Overall project maturity  [#########---] 75%
Validated foundation      [##########--] 85%
Adaptive control          [#####-------] 40%
```

## What This Repository Contains

- Nonlinear plant simulation using the `sin(theta)` model
- Six controller families:
  - `lqr_current`
  - `lqr_dob`
  - `hybrid_modern`
  - `paper_split`
  - `baseline_mpc`
  - `robust_hinf_like`
- Discrete DARE-based Kalman estimator
- Deterministic benchmark pipeline and Monte Carlo stress benchmark
- Telemetry CSV export for sim-to-real bring-up
- Model-validation and verification scripts

## Transparent Project Snapshot

This repository is in a useful but not fully finished state.

- The physical model, Jacobian, ZOH discretisation, estimator, and telemetry path are validated.
- The constrained MPC path is implemented, but on the current machine it falls back to `dlqr` because `quadprog` is not available on the MATLAB path.
- The new RLS-based online identification path is present behind `cfg.enable_online_id = false`, but the current verification case still fails when it is enabled.
- The benchmark score rebalancing now measures disturbance-axis transients, but nominal score separation is still much smaller than the target spread.

## Quick Start

```matlab
% From this folder in MATLAB:

unicycle_model_validate
unicycle_verify
unicycle_main
unicycle_stress_benchmark(20)
unicycle_sim2real_sensitivity
unicycle_telemetry_validate
```

## Main Files

| File | Purpose |
|---|---|
| `unicycle_main.m` | End-to-end design, simulation, benchmarking, and plots |
| `unicycle_config.m` | Central parameter and tuning configuration |
| `unicycle_plant.m` | Continuous-time linear model and nonlinear dynamics |
| `unicycle_design_controllers.m` | Controller design for all six families |
| `unicycle_simulate.m` | Closed-loop simulation loop |
| `unicycle_benchmark.m` | Benchmark scoring and export |
| `unicycle_verify.m` | Automated verification gate |
| `unicycle_model_validate.m` | Physics and estimator validation gate |
| `unicycle_diagnose.m` | Structured controller and observer diagnostics |
| `unicycle_telemetry_validate.m` | Telemetry sign/unit validation for HIL compatibility |
| `unicycle_stress_benchmark.m` | Monte Carlo disturbance robustness benchmark |
| `unicycle_sim2real_sensitivity.m` | Parameter-fragility sweep |

## Controllers

| Family | Intent |
|---|---|
| `lqr_current` | Baseline full-state LQR |
| `lqr_dob` | LQR with disturbance observer |
| `hybrid_modern` | LQR with pitch and roll integrators |
| `paper_split` | Decoupled Lee-style pitch/roll baseline |
| `baseline_mpc` | Finite-horizon MPC with runtime fallback to `dlqr` when QP support is unavailable |
| `robust_hinf_like` | Heavier-weight robust LQR approximation |

## Telemetry

Exported telemetry format:

```text
t, pitch_deg, roll_deg, dpitch, droll, omega_rw, tau_rw, tau_base
```

Validated sign convention:

- positive pitch = forward tilt
- positive roll = right tilt
- positive `tau_rw` opposes positive roll
- positive `tau_base` opposes positive pitch
- positive `omega_rw` is the wheel spin direction that counters positive roll

## Requirements

- MATLAB
- Control System Toolbox
- Optimization Toolbox only if you want the constrained `baseline_mpc` QP path to activate

## Known Issues

1. `unicycle_verify` currently passes 13/14, not 14/14.
2. The failing item is the online-ID survival test for `hybrid_modern` under a mass perturbation.
3. The nominal benchmark still shows compressed scores even after transient-aware reweighting.
4. The constrained MPC implementation is present, but this machine currently uses the fallback path because `quadprog` is unavailable.

## Upstream Lineage

This MATLAB package is derived from ideas and behavior in:

- `le-nicolas/Self-Balancing-Unicycle-Robot`
- `le-nicolas/reaction-wheel-balancer-sim2real`
