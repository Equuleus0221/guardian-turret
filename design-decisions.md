# Design decision log

ME319: Control Systems, Assignment 2 (Autumn 2026)
Author: Ryan Sammy D'Souza (24B2267)

The design in this repository was first proposed by an AI assistant (see the appendix of `paper/main.pdf`). This log records my own review of its ten principal design decisions: what the alternative was, why the choice holds, and where I agree or disagree.

## 1. Why estimate in Cartesian x–y and not directly in range and angle?

Cartesian coordinates allow target position, velocity and acceleration to be modeled using linear state equations. Estimating directly in polar coordinates introduces nonlinear relationships, requiring a more complicated filter. I agree, although converting radar measurements to Cartesian coordinates requires careful treatment of their transformed noise covariance.

## 2. Why an acceleration model and not constant velocity? What does 0.5 s of radar loss do to each?

A constant-velocity model assumes acceleration is zero, so an unmodeled acceleration of 20 m/s² can produce 2.5 m of position error over a 0.5 s blackout. The Singer model estimates acceleration and propagates it, reducing errors when maneuvers are reasonably predictable. I agree, but neither model can anticipate an unexpected acceleration reversal during a blackout.

## 3. Where do 10 m/s² and 2 s come from, and what happens if each is doubled?

The 10 m/s² value comes from assuming acceleration is uniformly distributed within a disc of radius 20 m/s², giving a per-axis standard deviation of 10 m/s². The 2 s correlation time is a heuristic based on maneuvers lasting several seconds. Doubling the acceleration standard deviation makes the filter more responsive but less smooth; doubling the correlation time makes modeled acceleration more persistent. I partly disagree: the 2 s choice needs a sensitivity analysis rather than relying mainly on intuition.

## 4. Why must the turret point at the present position and not 0.5 s ahead? What is 0.5 s for, then?

Since laser propagation time is negligible, firing at a target's future position would create an artificial lead error, potentially around 75 m at the maximum line-of-sight rate. The 0.5 s horizon instead allows the system to continue predicting target motion when radar measurements disappear. I agree, although predicting ahead by an appropriately calculated actuator-response time could improve tracking, unlike applying a fixed 0.5 s lead.

## 5. Why does the loop have to be type 2?

A type-1 loop has a finite steady-state error when tracking constant angular velocity, potentially requiring an impractically large velocity constant to meet the pointing tolerance. A type-2 loop eliminates this steady-state error and leaves a finite error for constant angular acceleration. I agree with the design choice, but disagree that type 2 is strictly necessary: sufficiently accurate feedforward could also make a type-1 system viable.

## 6. Why can't the inductance be neglected in the plant?

Neglecting motor inductance removes electrical dynamics that interact with the large reflected motor inertia, concealing the lightly damped resonance near 56 rad/s. This would make the simplified plant appear substantially easier to control than the actual system. I agree: retaining inductance is important because that resonance directly influences the achievable bandwidth and stability margins.

## 7. Why a notch, and what is its weakness?

The notch suppresses the resonance's contribution to the loop response, allowing higher low-frequency gain without destabilizing the controller. However, it relies on approximate pole-zero cancellation rather than physically damping the resonance, making it sensitive to model uncertainty. I agree, but would prefer additional robustness testing against changes in motor parameters and resonance frequency.

## 8. Why crossover at 15 rad/s and not 20 or 10?

A 20 rad/s crossover improves responsiveness but reduces gain margin significantly when the resonance shifts, while 10 rad/s sacrifices tracking performance and acceleration-error rejection. The selected 15 rad/s balances both requirements, yielding approximately 57° phase margin and 13.8 dB gain margin for the nominal plant. I agree, although testing the actual engagement performance at all three crossover frequencies would strengthen the justification.

## 9. Why is rate aiding acceptable alongside an LTI compensator, and would you keep it?

Rate aiding provides a feedforward voltage based on estimated line-of-sight velocity, reducing the tracking burden on the LTI feedback compensator without altering its nominal loop stability margins. The report shows that removing it increases mean interception time from approximately 0.69 s to 1.51 s in the comparison cases. I would keep it, provided noisy velocity estimates are adequately filtered and the combined system is validated under saturation.

## 10. Why gate firing on the filter covariance and not on a fixed threshold?

A fixed pointing-error threshold treats all estimates as equally reliable, even when uncertainty increases during radar loss. Covariance-based gating accounts for this uncertainty and prevents firing when the estimated miss distance is small but confidence is poor. I agree, although its reliability depends on covariance calibration; an overconfident filter could still authorize inaccurate shots.

## Related material

- Hand derivations of the plant model, error budget, the 0.5 s lead error and the compensator: `handwritten/`
- Design numbers and robustness sweep: `results/design_summary.txt`
- Evaluation and stress-test results: `results/evaluation_metrics.txt`, `results/stress_metrics.txt`
