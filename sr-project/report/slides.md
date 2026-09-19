# Slides

One figure per slide. Each title is the takeaway, so the deck reads as an
argument rather than a tour. Numbers match `report/numbers.md`.

---

## 1. A signal too weak to see

No figure.

The drive is A = 0.3 against a threshold of A_c = 0.3849, so without noise the
particle never leaves its well and the signal produces nothing. The usual
instinct is to remove noise.

---

## 2. Adding noise makes the signal appear

`figures/e1_regimes.png`

Three noise levels. Too little: 18.88 dB. Optimal: 28.11 dB. Too much:
20.68 dB. The middle panel tracks the drive; the others do not.

---

## 3. The effect needs no dynamics, only a threshold

`figures/e2_threshold.png`

A hard threshold with a sub-threshold input fires not at all at low noise,
peaks at 39.50 dB at sigma = 0.621, then degrades. Noise is what carries the
signal across.

---

## 4. The optimum lands where theory says

`figures/e3_snr_vs_noise.png`

Measured peak D = 0.1153 at 28.19 dB against a predicted dU/2 = 0.1250. A 7.7 %
error, smaller than the 14.6 % grid step.

---

## 5. The result is not an artefact of the solver

`figures/e3_numerical_checks.png`

Halving the step and switching to Heun both leave the peak at D = 0.1153, with
maximum deviations of 0.69 dB and 0.66 dB.

---

## 6. The optimum ignores amplitude and follows frequency

`figures/e4_robustness.png`

SR peak stays near D = 0.10 to 0.13 for A = 0.1, 0.2 and 0.3, but rises from
0.0879 to 0.1986 as f0 goes from 0.005 to 0.05. Faster signal, faster hopping,
more noise.

---

## 7. Tuning the system instead of the noise does not work here

`figures/e5_parameter_tuned.png`

Peaks sit at a = 0.25 to 0.39 regardless of input noise and do not follow the
prediction. The plan's a = 8 D_in is also a D-sweep condition misapplied; the
correct stationary point is a = 4 D_in.

---

## 8. Noise turns a useless detector into a perfect one

`figures/e6_detection.png`

The threshold detector goes from AUC 0.500 and Pd 0.000, through AUC 1.000 and
Pd 1.000 at sigma = 0.156, back to AUC 0.495 and Pd 0.046. Chance to perfect to
chance, on noise alone.

---

## 9. But a linear filter still wins

`figures/e7_baselines.png`

Adding noise is worth +8.40 dB to the bistable system. The matched linear
band-pass beats that by a further +13.62 dB.

---

## 10. What the low-noise end really shows

`figures/e4_robustness.png`, left panel.

The two-state model has no intra-well motion, so it predicts SNR vanishing as
noise goes to zero. A real particle still responds linearly inside one well,
which is why four separate experiments show their largest value at the lowest
noise rather than at the resonance.

---

## 11. Limits

No figure.

Stochastic resonance helps threshold-type nonlinear detectors, not linear
filtering of a known sine. The optimum needs tuning per signal. The two-state
theory is quantitative only where hopping dominates, and its prefactor is
unverified, leaving a 2.30 dB offset unresolved.
