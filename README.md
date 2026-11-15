# Guardian turret: predictive tracking and elevation control

ME319 Control Systems (Autumn 2026), Assignment 2: *Mumbai 2044: Predictive Tracking and Turret Control for a Maneuvering Aerial Target*.

The repository contains the estimator, predictor, feedback controller and fire logic for the elevation axis of a laser turret, the plant model, the target/radar generator, and the scripts that produce every number and figure in the paper (`paper/main.pdf`).

## Author's verification and design review

The design was developed with AI assistance, documented in the appendix of the paper. The author reviewed the principal design decisions and recorded his reasoning on each (`design-decisions.md`) and re-derived parts of the model and controller analysis by hand (`handwritten/`).

## Quick start

Requirements: MATLAB (developed and tested on R2025a). No toolboxes are needed.

```matlab
>> run_all        % about 30 minutes; writes ./results and ./figures
```

Or step by step:

| Script | What it does | Outputs |
|---|---|---|
| `design_controller.m` | Checks `turret_plant.m` against the transfer function, builds the compensator, computes poles, margins, robustness | `results/design_summary.txt`, `figures/fig_design.pdf` |
| `run_scenarios.m` | Scenarios A, B, C: one plotted run each plus 50 noise realisations each | `results/scenario_metrics.txt`, `results/tab_scenarios.tex`, `figures/fig_scenarios.pdf`, `figures/fig_geometry.pdf` |
| `run_final_evaluation.m` | 200 randomized held-out engagements and four design alternatives | `results/evaluation_metrics.txt`, `results/tab_ablation.tex`, `results/evaluation.mat`, `figures/fig_evaluation.pdf` |
| `run_stress_tests.m` | 670 engagements with an independent evaluator: other manoeuvre types and calling conventions | `results/stress_metrics.txt`, `results/stress.mat` |

A single engagement:

```matlab
>> L = simulate_engagement('C');   % 'A', 'B', 'C', or a case number such as 137
>> L.M                              % metrics: interception time, miss distances, ...
>> plot(L.t, L.miss)                % beam-target distance over time
```

## The three prescribed functions

```matlab
[r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
[V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)
[theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i)
```

`turret_controller.m` and `turret_plant.m` are **self-contained** (no other file from this repository is needed), so they can be dropped next to an external evaluation script.

Conventions the controller relies on:

- It keeps persistent state. A call with `t = 0` starts a new engagement.
- It may be called at any fixed or variable step; it measures the step from `t`. Results in the paper use a 1 ms step; steps from 0.5 ms to 10 ms were checked on 20 cases with the same outcome (every pulse a hit).
- `meas_valid = 1` marks a new radar sample. A sample is used only if its values differ from the last one used, so a caller that holds `meas_valid = 1` between samples is handled too.
- `fire_cmd` rises to 1 for 0.1 s when a pulse is requested and then returns to 0, so each request is one 0-to-1 transition. Requests are spaced at least 1.0 s apart.

## Design in brief

| Block | Choice | Key parameters (all in `turret_controller.m`, function `ctrl_params`) |
|---|---|---|
| Estimator | Kalman filter in Cartesian coordinates, Singer acceleration model, converted radar measurements | `sigma_m = 10 m/s^2`, `tau_m = 2 s` |
| Predictor | Open-loop propagation of the same model; trajectory valid for `t_p = 0.5 s` | `t_p = 0.5` |
| Reference | Elevation of the predicted position due at the present time, clipped 0.5 deg inside the stops | `t_lead = 0` |
| Feedback | `Gc(s) = Kp (s + z_i)/s * (s^2 + 2 zeta_p wn s + wn^2)/(s^2 + 2 zeta_d wn s + wn^2)` | `Kp = 14748 V/rad`, `z_i = 3 rad/s`, `wn = 55.94 rad/s`, `zeta_d = 0.7` |
| Rate aiding | Back-emf voltage of the predicted line-of-sight rate, outside the loop | `use_ff = true` |
| Limits | Limiter upstream of the notch; integrator active for small errors only | `V_lim = 900 V`, `V_max = 1000 V`, `e_int = 2 mrad` |
| Fire logic | Aim gate, dwell check and 3-sigma confidence gate from the filter covariance | `d_aim = 0.5 m`, `kappa = 3`, `tau_fire = 0.02 s` |

Plant: `Gp(s) = 1000 / (320 s^3 + 3350 s^2 + 1001500 s)`, poles at 0 and -5.23 +/- j55.7 rad/s.
Loop: crossover 15 rad/s, phase margin 56.7 deg, gain margin 13.8 dB, `Ka = 44.2 1/s^2`.

The 1000 V drive limit is a design assumption; the assignment does not specify one.

## Scenarios

`scenario_definition.m` defines all cases.

- `'A'`: straight pass at 100 m/s.
- `'B'`: A plus a two-axis sinusoidal maneuver from t = 1 s (peak 19 m/s^2).
- `'C'`: B with the radar lost from 1.0 s to 1.5 s.
- `n` (a number): randomized case `n` with random geometry, a two-axis maneuver at 60-98 % of the acceleration limit, and one or two radar losses of 0.2-0.5 s. Cases 1-40 were used for design decisions; cases 101-300 are the held-out evaluation set.
- `struct('id', 'B', 'seed', 7)`: the same motion with another noise realisation.

Every trajectory is checked against 800 m <= r <= 2000 m, speed <= 150 m/s and acceleration <= 20 m/s^2. Scenario and noise generation use a small portable random generator, so results do not depend on the MATLAB random stream.

## Files

```
turret_controller.m        estimator, predictor, compensator, fire logic (self-contained)
turret_plant.m             elevation-drive dynamics with stops (self-contained)
target_motion_generator.m  radar model: 20 Hz, noise, measurement loss
target_truth.m             true target kinematics
scenario_definition.m      scenarios A, B, C and randomized cases
simulate_engagement.m      closed-loop simulation (RK4, 1 ms) and metrics
loop_analysis.m            margins, poles, step response (no toolbox)
design_controller.m        design calculations and Fig. 2
run_scenarios.m            standardized scenarios
run_final_evaluation.m     held-out evaluation and design alternatives
plot_evaluation.m          evaluation figure from results/evaluation.mat
run_stress_tests.m         stress suite (other manoeuvres, calling conventions)
stress_trajectory.m        stress-test target trajectories
stress_engagement.m        independent evaluator used by the stress suite
run_all.m                  runs the three scripts above
save_fig.m, shade.m, tee.m, p95.m   small helpers
design-decisions.md        author's reasoning on the ten principal design choices
handwritten/               author's hand derivations (photos)
results/                   text summaries and LaTeX table rows
figures/                   PDF and PNG figures
paper/                     IEEE paper source (main.tex) and PDF
```

## Rebuilding the paper

```
cd paper
pdflatex main
pdflatex main
```

`paper/numbers.tex` holds the simulation numbers quoted in the text. After re-running the scripts, update it from `results/*.txt`; the two tables are read from `results/` automatically.

## Results at a glance

| Test | Intercepted | Pulses hit | Mean interception time | Mean miss at pulse |
|---|---|---|---|---|
| Scenarios A, B, C (50 noise runs each) | 150/150 | 1500/1500 | 0.50 s | 0.32-0.39 m |
| 200 held-out randomized cases | 200/200 | 1964/1964 | 0.68 s | 0.35 m |
| Stress: 7 other manoeuvre families (40 each) | 280/280 | 2661/2661 | 0.75-1.54 s | 0.25-0.45 m |
| Stress: 13 calling conventions (30 each) | 390/390 | 3784/3784 | 0.69-1.27 s | 0.36-0.42 m |

Full numbers: `results/scenario_metrics.txt`, `results/evaluation_metrics.txt`, `results/stress_metrics.txt`, `results/design_summary.txt`.

### Stress tests

`run_stress_tests.m` uses a second evaluator (`stress_engagement.m`) that shares no code with `simulate_engagement.m`: its own radar model (MATLAB `randn`), integration loop and hit test.

- Manoeuvre families (`stress_trajectory.m`): straight passes in any direction, bang-bang jinks, sustained near-limit turns with reversals, a full-acceleration reversal at the onset of a 0.5 s radar loss, a 150 m/s crosser at 810-900 m, engagements up to 84 deg elevation, and four or five radar losses per engagement.
- Calling conventions: Euler or RK4, steps of 1-10 ms, the controller called inside every RK4 stage, `meas_valid` held high between samples, NaN measurements during a loss, the turret starting on either stop, a shifted radar clock, and a first call after t = 0.

The first run of this suite found a fault: with `meas_valid` held high, a repeated pre-loss measurement was fused as current after a radar loss (88 of 569 pulses missed in the two affected conventions). The controller now treats a sample as new only if its values differ from the last one used; the table above is the run after that change.

## Notes and limitations

- Octave: the three prescribed functions avoid toolbox and MATLAB-only calls, but the repository has only been run in MATLAB.
- The instructor's evaluation cases are not available before submission. The "final evaluation" and the stress tests here are our own test sets.
- The stress suite uses MATLAB's `RandStream`, so it is MATLAB-only.
- AI use is documented in the appendix of the paper.
