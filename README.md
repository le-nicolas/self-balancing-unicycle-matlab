# Self-Balancing Unicycle MATLAB

## What this is

A MATLAB simulation of a self-balancing unicycle robot with a reaction wheel on the roll axis and a driven base wheel on the pitch axis. This package was converted from two Python/MuJoCo repositories:

- `le-nicolas/Self-Balancing-Unicycle-Robot`
- `le-nicolas/reaction-wheel-balancer-sim2real`

Its main job is control design, verification, benchmarking, and firmware-parameter export for an ESP32-based hardware build. It is not a real-time simulator.

Hardware target:

- ESP32-WROOM dev board
- BMI088 IMU
- AS5600 magnetic encoder on the reaction wheel BLDC
- 2804 hollow-shaft BLDC + DRV8313 driver for the reaction wheel
- 110 RPM geared base motor + BTS7960 H-bridge
- 3S LiPo 11.1 V

## System model

State vector:

`x = [theta_p, dtheta_p, theta_r, dtheta_r, omega_rw]^T`

- `theta_p` = pitch angle, positive when the robot tilts forward
- `theta_r` = roll angle, positive when the robot tilts right
- `omega_rw` = reaction wheel speed

Control inputs:

`u = [tau_rw, tau_base]^T`

- `tau_rw` = reaction wheel torque, positive when it opposes positive roll
- `tau_base` = base wheel torque, positive when it opposes positive pitch

Nonlinear equations of motion used by the simulator:

```text
(I_body + m*L^2) * ddtheta_p = m*g*L*sin(theta_p) - tau_base
I_body_r * ddtheta_r         = m*g*L*sin(theta_r) - tau_rw
I_rw * domega_rw             = tau_rw
```

Linearisation is intended for small angles and is verified in closed loop to `0.000083 deg` maximum pitch error in `unicycle_model_validate.m` Stage 3.

Key parameters from `unicycle_config.m`:

- `m_body = 0.85 kg`
- `L_body = 0.28 m`
- `dt = 0.005 s` (`200 Hz`)
- `tau_rw_max = 0.25 N*m`
- `tau_base_max = 0.40 N*m`
- `omega_rw_max = 800 rad/s`
- `omega_rw_hard = 950 rad/s`
- crash threshold: `abs(theta_p)` or `abs(theta_r)` greater than `35 deg`

## File map

Core simulation

| File | Purpose |
|---|---|
| `unicycle_main.m` | End-to-end entry point: controller design, simulation, benchmark, plots |
| `unicycle_config.m` | Central parameter, tuning, disturbance, and safety configuration |
| `unicycle_plant.m` | Continuous-time linear model and nonlinear plant equations |
| `unicycle_discretize.m` | Exact ZOH discretisation helper |
| `unicycle_design_controllers.m` | Designs all controller families and observer gains |
| `unicycle_simulate.m` | Closed-loop nonlinear simulation loop with estimator, limits, and crash detection |

Benchmarking

| File | Purpose |
|---|---|
| `unicycle_benchmark.m` | Benchmark scoring and result aggregation |
| `unicycle_plot_all.m` | Time histories, phase portrait, comparison plots, and bar charts |
| `unicycle_stress_benchmark.m` | Monte Carlo disturbance robustness benchmark |
| `unicycle_calibrate_kicks.m` | Kick-amplitude sweep utility for deterministic benchmark calibration |

Verification

| File | Purpose |
|---|---|
| `unicycle_verify.m` | Automated verification gate, now `17/17 PASS` |
| `unicycle_model_validate.m` | Physics/model validation gate, now `8/8 PASS` |
| `unicycle_diagnose.m` | Structured diagnostics for plant, observer, disturbance, and MPC checks |

Sim-to-real

| File | Purpose |
|---|---|
| `unicycle_sim2real_sensitivity.m` | Sensitivity sweep for sim-to-real fragility ranking |
| `unicycle_telemetry_validate.m` | Validates telemetry column order, units, and sign conventions |
| `unicycle_hil_smoke.m` | Synthetic HIL smoke test with 7 sign/ESTOP/watchdog scenarios |
| `unicycle_rls.m` | Pitch-axis RLS estimator used by optional online ID |
| `unicycle_export_firmware.m` | Exports firmware-ready JSON and C header files |
| `unicycle_test_firmware_parity.m` | Verifies exported params reproduce MATLAB control to float precision |

Other source files

| File | Purpose |
|---|---|
| `firmware_export/` | Generated export directory containing `firmware_params.json` and `firmware_params.h` |

Note: `firmware_export/firmware_params.json` and `firmware_export/firmware_params.h` are generated outputs, not source files.

## Quick start

Requirements:

- MATLAB R2020b or newer
- Control System Toolbox
- Optional: Optimization Toolbox plus `quadprog` on the MATLAB path to activate constrained MPC

```matlab
% 1. Verify the model and simulation chain
unicycle_model_validate   % must print OVERALL: PASS (8/8)
unicycle_verify           % must print PASS: 17 / 17

% 2. Run the full benchmark
unicycle_main             % benchmark, plots, CSV to ./results/

% 3. Export firmware parameters
cfg = unicycle_config();
[A, B, C, D] = unicycle_plant(cfg);
ctrl = unicycle_design_controllers(A, B, C, cfg);
unicycle_export_firmware(ctrl, cfg, './firmware_export/');

% 4. Run HIL smoke test
results = unicycle_hil_smoke(ctrl, cfg);

% 5. Stress benchmark (seed=42, 20 episodes)
unicycle_stress_benchmark(20)
```

## Controller families

| Name | MATLAB struct field | Description | Exported? |
|---|---|---|---|
| `lqr_current` | `lqr_current` | Delta-u LQR plus discrete Kalman observer | yes |
| `lqr_dob` | `lqr_dob` | LQR plus first-order disturbance observer | yes |
| `hybrid_modern` | `hybrid_modern` | LQR with pitch and roll integrators (7-state) | no |
| `paper_split` | `paper_split` | Decoupled pitch/roll baseline following Lee 2013 structure | yes |
| `baseline_mpc` | `baseline_mpc` | Constrained MPC when `quadprog` is available, otherwise `dlqr` fallback | no |
| `robust_hinf_like` | `robust_hinf_like` | H-inf-weighted LQR approximation | yes |

Firmware-export notes:

- `hybrid_modern` is excluded because its 7-state integrator path needs persistent integrator state management and ESTOP reset logic in firmware.
- `baseline_mpc` is excluded because it requires a QP solve at `200 Hz`.
- Constrained MPC activates only when both `license('test','optimization_toolbox') == 1` and `exist('quadprog','file') == 1`.

## Verification gates

Two files must pass before trusting results.

`unicycle_model_validate.m` (`8 stages`)

- Stage 1: symbolic linearisation check, Jacobian match below `1e-10`
- Stage 2: equilibrium and sign checks for gravity and input directions
- Stage 3: linear versus nonlinear agreement, `0.000083 deg` max pitch error in closed loop
- Stage 4: exact ZOH discretisation and `Bd` sign audit
- Stage 5: step-response sign checks for `tau_rw` and `tau_base`
- Stage 6: open-loop divergence and closed-loop convergence
- Stage 7: pitch/roll axis decoupling
- Stage 8: Kalman observer audit, about `1.2 mrad` RMS tracking and `17.8 mrad` bias error

`unicycle_verify.m` (`17 tests`)

- Tests 1-10: core correctness, controllability, stability, tracking, and saturation
- Tests 11-12: MPC constraint satisfaction, skipped cleanly when `quadprog` is unavailable
- Test 13: transient metrics are non-zero and correctly order MPC above paper_split
- Test 14: RLS online ID survives a `+20%` pitch inertia perturbation
- Tests 15-16: firmware export validity and double-versus-single parity (`6e-7 N*m` scale)
- Test 17: HIL smoke test, `7/7` scenarios pass

Current gate status:

- `unicycle_model_validate`: `OVERALL: PASS (8/8)`
- `unicycle_verify`: `PASS: 17 / 17`

## Benchmark results

These numbers are from the MATLAB nonlinear simulation using the `sin(theta)` plant, not MuJoCo.

Current deterministic nominal benchmark configuration:

- `T_sim = 20 s`
- deterministic kicks at `t = [3.0, 6.0, 10.0] s`, alternating sign
- `dist_amp = 0.34 N*m`
- sustained sine disturbance from `t = 13 s` to `20 s`
- sine amplitude `0.04 N*m`
- sine frequency `0.8 Hz`

Nominal benchmark:

| Controller | Composite score |
|---|---:|
| Baseline MPC (unconstrained) | 94.348 |
| Robust H-inf-like | 94.216 |
| LQR + DOB | 93.928 |
| Hybrid Modern (LQR+I) | 93.927 |
| LQR Current | 93.925 |
| Paper Split Baseline | 93.674 |

Stress benchmark, seed `42`, `20` episodes:

| Controller | Survival fraction |
|---|---:|
| Baseline MPC (unconstrained) | 0.950 |
| Hybrid Modern (LQR+I) | 0.850 |
| LQR Current | 0.850 |
| Paper Split Baseline | 0.900 |
| LQR + DOB | 0.800 |
| Robust H-inf-like | 0.750 |

Note on stress ranking:

`paper_split` looks unusually strong in the current stress run relative to its nominal score. In this model that is physically plausible: its weaker gains drive the reaction wheel less aggressively, so it sometimes avoids wheel-speed saturation that causes stronger controllers to lose a recovery. The stress benchmark measures disturbance survival, not overall control quality.

## Sim-to-real workflow

Step 1 - Run the verification gates

- `unicycle_model_validate`
- `unicycle_verify`

Both must pass before continuing.

Step 2 - Run the HIL smoke test

- `unicycle_hil_smoke(ctrl, cfg)`

All 7 scenarios must pass before trusting sign conventions or ESTOP behavior.

Step 3 - Export firmware parameters

- `unicycle_export_firmware(ctrl, cfg, './firmware_export/')`
- outputs: `firmware_params.json` and `firmware_params.h`
- verify export parity with `unicycle_test_firmware_parity`

Step 4 - Copy to ESP32 firmware

- copy `firmware_params.h` into `esp32_rw_base_platformio/src/`
- select the controller family at compile time via `#define CTRL_FAMILY`
- recommended bring-up order: `lqr_current`, then `lqr_dob`, then `robust_hinf_like`

Step 5 - Hardware bring-up

- power the robot with motors restrained and wheels lifted
- verify IMU signs match `unicycle_telemetry_validate`
- verify roll commands go to the reaction wheel and pitch commands go to the base wheel
- verify ESTOP at the configured crash thresholds
- save the mapping profile before enabling full closed-loop control

Step 6 - Sim-to-real sensitivity

- run `unicycle_sim2real_sensitivity()`
- before touching hardware, run `unicycle_sim2real_sensitivity()` one more time and record the two most fragile parameters
- in the current model those are body mass at about `+/-15%` and CoM height at about `+/-10%`
- if hardware behavior diverges from prediction, update `m_body` and `L_body` first in `unicycle_config.m`
- then rerun `unicycle_export_firmware`, re-flash the ESP32, and compare again
- that is the intended feedback loop between the MATLAB design and the physical robot

## Known limitations

1. Physics fidelity

This simulator uses rigid-body dynamics with `sin(theta)` nonlinearity. It does not model contact forces, wheel slip, rolling friction, backlash, or structural flex. The original project used MuJoCo. Disturbance rejection results here are optimistic compared to hardware, especially under large disturbances.

2. Nominal benchmark score compression

The current composite score formula

```text
100*survival - 50*pitch_RMS - 40*roll_RMS - 5*ctrl_RMS - 20*peak - 15*settling
```

still produces only about `0.67` points of spread across the top five controllers in the nominal run. Those controllers all drive the state close to the estimator-noise floor, so asymptotic RMS terms remain tightly compressed. The stress benchmark is the more meaningful discriminator.

3. Constrained MPC availability

The constrained MPC branch is implemented in `unicycle_design_controllers.m` and `unicycle_simulate.m`, but it requires Optimization Toolbox support and `quadprog` on the MATLAB path. On the current machine, the code falls back to `dlqr` because `quadprog` is unavailable.

4. hybrid_modern firmware exclusion

`hybrid_modern` is not exported to `firmware_params.h`. The ESP32 scaffold does not yet implement persistent integrator state management or integrator reset on ESTOP. Exporting the `2x7` gain matrix without that runtime logic would be unsafe.

5. RLS adaptive ID scope

The RLS estimator identifies only the pitch subsystem. It estimates stiffness and pitch control authority, then updates only the pitch-related gain columns. Roll dynamics are not identified online. That is conservative, but it means large roll-axis payload asymmetry can still degrade performance.

6. No real-time HIL bridge in MATLAB

The Python sim-to-real repository includes a real-time UDP HIL bridge. This MATLAB repo does not. `unicycle_hil_smoke.m` is a synthetic smoke test against the simulation plant, useful for sign-convention and ESTOP sanity checks, but it does not exercise live hardware.

## Hardware bring-up checklist

1. Confirm `unicycle_model_validate` prints `OVERALL: PASS (8/8)`.
2. Confirm `unicycle_verify` prints `PASS: 17 / 17`.
3. Confirm `unicycle_hil_smoke` reports `PASSED: 7 / 7`.
4. Confirm `unicycle_export_firmware` writes both `firmware_params.json` and `firmware_params.h`.
5. Confirm `unicycle_test_firmware_parity` reports both channels below `1e-4 N*m`.
6. Copy `firmware_params.h` into `esp32_rw_base_platformio/src/`.
7. Power the robot with motors restrained and wheels lifted.
8. Verify IMU pitch sign: forward tilt must produce positive `pitch_deg`.
9. Verify IMU roll sign: rightward tilt must produce positive `roll_deg`.
10. Verify reaction wheel direction: positive `tau_rw` command produces the expected positive wheel spin direction for your mapping.
11. Verify base wheel direction: positive `tau_base` command produces forward corrective wheel torque.
12. Force a tilt past `35 deg` and verify zero motor output.
13. Kill the sensor feed and verify timeout/ESTOP behavior.
14. Save `hardware_mapping_profile.json` with the verified signs.
15. Run the first closed-loop hardware test with `lqr_current` and the robot physically restrained.
16. Check that pitch settling is on the order of `150 ms`, matching the current simulation bandwidth.
17. Graduate to `lqr_dob` and then `robust_hinf_like` only after `lqr_current` is stable.

## License

This repository does not currently include a standalone license file. Until one is added, treat it as all rights reserved by default.

The code and workflow here were derived from earlier work in:

- `le-nicolas/Self-Balancing-Unicycle-Robot`
- `le-nicolas/reaction-wheel-balancer-sim2real`

If you plan to redistribute or reuse this code outside your own development environment, add an explicit project license first.
