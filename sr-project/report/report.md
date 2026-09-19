# Stochastic resonance in a bistable system

Every number in this report is reproduced in `report/numbers.md`, which
`report/extract_numbers.py` generates from `results/`. Figures come from
`figures/`, each drawn by the experiment script named in the caption.

## 1. Introduction

A signal too weak to cross a detector's threshold is invisible, and the usual
instinct is to remove noise. Stochastic resonance is the observation that the
opposite can help: in a nonlinear system with a threshold, adding noise can
make a sub-threshold signal easier to see, up to an optimal noise level, beyond
which performance falls again.

This report tests that claim on the standard model, an overdamped particle in a
double-well potential driven by a weak sine. The scope is deliberately narrow.
We measure where the optimum sits, how it moves with the signal, whether it
improves detection rather than just a spectral ratio, and how it compares with
doing nothing and with a linear filter. We do not claim stochastic resonance is
a good engineering choice; section 8 argues it usually is not for this problem.

## 2. Theory

The model is

    dx/dt = a x - b x^3 + A cos(2 pi f0 t) + sqrt(2 D) xi(t)

with xi unit white Gaussian noise. The potential U(x) = -a x^2 / 2 + b x^4 / 4
has minima at x_m = +-sqrt(a/b) separated by a barrier dU = a^2 / (4 b). At the
baseline a = b = 1 the barrier is **0.2500** and the static threshold amplitude
is **A_c = 0.3849**. The drive amplitude is **A = 0.3**, which is **0.779** of
A_c, so without noise the particle can never leave its well.

Noise lets it escape at the Kramers rate
r_K = a / (sqrt(2) pi) exp(-dU / D). The two-state model of McNamara and
Wiesenfeld gives an output SNR proportional to (A x_m / D)^2 r_K, which peaks
where d/dD of its logarithm vanishes, at **D = dU/2 = 0.1250**. That value is
the prediction this project tests, and it is independent of the prefactor.

Two caveats, both of which matter later. The two-state model describes only
well-to-well hopping and contains no intra-well motion, so it says nothing
about what happens when hopping is rare. And its prefactor is unverified: see
section 8.

## 3. Method

The Langevin equation is integrated with Euler-Maruyama at dt = 0.01, with the
noise increment scaled as sqrt(2 D dt); scaling by dt instead of sqrt(dt) is
the classic error and is guarded by a unit test. The solver is vectorised over
noise intensity and trial, so a whole sweep steps as one array.

Output is decimated to fs = 10 and recorded over an exact integer number of
signal periods, which puts f0 on FFT bin 128 with no leakage. The SNR estimator
averages periodograms across trials, then takes the ratio of the f0 bin to the
median of the surrounding background bins. Error bars are bootstrap 16th and
84th percentiles over trials.

Two numerical checks back the solver. Halving dt moves the E3 peak not at all,
still **D = 0.1153**, with a maximum deviation of **0.69 dB** across the whole
curve. Switching to Heun, a stochastic Runge-Kutta scheme, likewise leaves the
peak at **D = 0.1153** with a maximum deviation of **0.66 dB**. Figure
`e3_numerical_checks` shows all three curves.

## 4. Results: stochastic resonance exists

**Figure `e1_regimes`.** Three regimes in the time domain. At D = 0.02 the
particle hops rarely and out of step with the drive, giving **18.88 dB**. At
D = 0.12 it tracks the drive, giving **28.11 dB**. At D = 1.0 hopping is faster
than the signal and coherence is lost, giving **20.68 dB**. The middle panel is
the effect in one picture.

**Figure `e2_threshold`.** The same effect without any dynamics. A hard
threshold at theta = 1 with a sub-threshold drive A = 0.7 produces nothing at
all for the **4** lowest noise levels, rises to **39.50 dB** at
**sigma = 0.621**, and falls to **30.58 dB** by sigma = 3. Noise alone carries
the signal across.

**Figure `e3_snr_vs_noise`, the headline.** Output SNR against noise intensity
peaks at **D = 0.1153** with **28.19 dB**, bootstrap band
**[28.04, 28.34] dB**, against the predicted 0.1250. That is an error of
**7.7 %**, comfortably inside the 30 % the project set as its criterion, and
smaller than the sweep can resolve: the 40-point log grid steps **14.6 %** per
point, so the measured peak is the nearest grid value to the prediction. The
curve runs from **20.67 dB** at D = 0.01 up to the peak and down to
**18.13 dB** at D = 2.

## 5. Results: how the optimum behaves

**Figure `e4_robustness`.** The optimum does not move with amplitude: the SR
peak sits at **D = 0.1321**, **0.1321** and **0.1007** for A = 0.1, 0.2 and
0.3. Two-state theory predicts exactly that, since dU/2 has no A in it.

The optimum does move with frequency, rising monotonically: **D = 0.0879**,
**0.1153**, **0.1321** and **0.1986** for f0 = 0.005, 0.01, 0.02 and 0.05. This
is time-scale matching. A faster drive needs a faster hopping rate to stay
synchronised, and hopping is only made faster by more noise. The heatmap panel
shows the same trend across a two-dimensional grid.

One caution on reading that figure. For A = 0.1 and A = 0.2 the largest SNR in
the sweep is not the resonance: it is **31.81 dB** and **37.16 dB** at the
lowest noise, D = 0.01. That is the intra-well response discussed in section 8,
and the figure marks the resonance rather than the global maximum.

**Figure `e5_parameter_tuned`.** When the noise is fixed and the system is
tuned instead, the measured optimum does not follow the prediction. Sweeping
a with b = a at three fixed input noise levels, the interior peaks sit at
**a = 0.394**, **0.281** and **0.251** against predicted 0.20, 0.40 and 0.80,
and at **0.281**, **0.394** and **0.281** at the weaker amplitude against 0.40,
0.80 and 1.20. The peaks do not move right as D_in grows. Section 8 explains
why, and corrects the prediction itself.

## 6. Results: detection and baselines

SNR is a spectral ratio over a long record. A detector sees one short record
and must answer yes or no, so section 6 measures that directly: 500 runs per
hypothesis, statistic the f0 periodogram power of a single 16-period record.

**Figure `e6_detection`.** For the bistable system the answer is layered. At
the baseline amplitude the problem is trivial, AUC **0.990 to 1.000** across
every noise level, so the sweep runs at a weaker A = 0.05 where the metric
carries information. There, AUC starts at **1.000**, collapses to a dip of
**0.665** at D = 0.0533 as hopping begins, recovers to **0.792** at
D = 0.1626, and declines to **0.567** at D = 2. Detection probability at a 5 %
false-alarm rate does the same more sharply, dipping to **0.082** and
recovering to **0.346**, a gain of **+0.264**. So noise does improve detection,
at the predicted place, but only relative to the dip.

For the threshold detector there is no ambiguity. It has no intra-well
response, so below threshold it emits nothing and both hypotheses give an empty
record: AUC is exactly **0.500** and detection probability **0.000**. Adding
noise takes it to AUC **1.000 at sigma = 0.156** and detection probability
**1.000**, and yet more noise returns it to chance, AUC **0.495 at sigma = 300**
and detection probability **0.046**. Chance to perfect to chance, driven by
nothing but noise.

**Figure `e7_baselines`.** Five detectors on the same low-SNR input.

| Method | SNR | AUC | Pd at 5 % |
| --- | --- | --- | --- |
| Linear band-pass on the raw input | **41.60 dB** | 1.000 | 1.000 |
| Threshold detector, no added noise | 30.74 dB | 0.939 | 0.724 |
| Stochastic resonance at the optimal D | 27.99 dB | 0.779 | 0.320 |
| Parameter-tuned stochastic resonance | 25.79 dB | 0.980 | 0.998 |
| Bistable system, no added noise | 19.58 dB | 1.000 | 1.000 |

Adding noise to the bistable system is worth **+8.40 dB**, from 19.58 to 27.99.
That is stochastic resonance doing real work. The linear band-pass still beats
it by **+13.62 dB**. The band-pass is not a weak comparison: for a sine of
known frequency in white Gaussian noise, the f0 periodogram bin is the matched
filter, so it is the strongest linear detector available for this problem.

## 7. Extension: low-light images

Not attempted. `src/srlab/image_sr.py` is a documented stub and no image
experiment was run, so this project makes no claim about images.

## 8. Discussion

**Where the two-state model breaks, and why it matters four times over.** The
model contains only well-to-well hopping. It has no intra-well motion. But a
real particle sitting in one well still responds linearly to the drive, with an
amplitude of roughly A / (2a), and that response is perfectly coherent. As noise
falls, the background falls with it while that coherent tone does not, so the
measured SNR rises as D approaches zero, exactly where the model predicts it
should vanish.

This single omission explains four separate observations. E3's curve turns
upward at its low-noise end instead of falling away. E4 reports its largest SNR
at the lowest noise for the two weaker amplitudes rather than at the resonance.
E5's curves decline monotonically as a falls, because shrinking a both lowers
the barrier and softens the well, and the softer well gives a larger intra-well
response. E6's detection is trivially easy at low noise for the same reason. In
every case the resonance is still present and still near dU/2; it is simply not
the global maximum.

**A correction to the parameter-tuned prediction.** The project plan gives the
tuning optimum as dU = 2 D_in, hence a = 8 D_in. That is the condition for a
sweep over D, where the peak arises from balancing a 1/D^2 prefactor against
exp(-dU/D). Sweeping a at fixed D, x_m = sqrt(a/b) = 1 is constant and the only
a dependence is SNR proportional to a exp(-a / (4 D_in)), whose stationary
point is **a = 4 D_in**, that is dU = D_in. This was confirmed against the
two-state formula itself, whose maximum over a sits at exactly 0.20, 0.40, 0.80
and 1.20 for D_in = 0.05, 0.1, 0.2 and 0.3. The figure marks the corrected
value. Even so the measurement does not follow it, for the intra-well reason
above.

**An unverified prefactor.** The two-state SNR is implemented with the
prefactor pi that the project plan specifies. An independent derivation, taking
the signal as the delta-peak strength and the background as the one-sided noise
PSD, gives pi/2 instead, a factor of two and therefore 3 dB. The difference is a
spectrum convention, not physics, and it does not move the peak. The measured
offset between theory and simulation at the peak is **-2.30 dB**, which
excludes neither: with pi the theory sits 2.30 dB above the measurement, with
pi/2 it would sit 0.70 dB below. **TODO:** settle this against Gammaitoni,
Hanggi, Jung and Marchesoni, Rev. Mod. Phys. 70, 223 (1998), which was not
reachable while this work was done. No equation number is cited here rather
than inventing one.

**Three statements the evidence supports.**

1. Added noise improves a nonlinear, threshold-type detector. The threshold
   detector goes from useless to perfect on noise alone, and the bistable
   system gains **+8.40 dB**. It does not beat an optimal linear filter for a
   sine in white Gaussian noise, which wins by **+13.62 dB**.
2. The optimum is not a universal constant. It is independent of amplitude but
   rises with frequency, from **D = 0.0879** at f0 = 0.005 to **0.1986** at
   f0 = 0.05, so any practical use needs tuning to the signal.
3. The two-state theory is an approximation valid where hopping dominates. It
   predicts the peak location to **7.7 %**, within one grid step, and fails
   qualitatively at low noise where intra-well motion takes over.

## 9. Conclusion

An overdamped bistable system driven by a sub-threshold sine shows a clear
maximum in output SNR at **D = 0.1153**, within **7.7 %** of the predicted
dU/2 = 0.1250 and within the grid resolution, confirmed by two independent
numerical schemes agreeing to better than **0.7 dB**. The optimum is
independent of signal amplitude and rises with signal frequency, as time-scale
matching requires. Noise measurably improves detection, taking a threshold
detector from chance to certainty and back, but for a sine of known frequency
in white Gaussian noise a linear matched filter remains **13.62 dB** better
than the best stochastic-resonance configuration tested.

### References

1. Benzi, Sutera and Vulpiani, J. Phys. A 14, L453 (1981).
2. Gammaitoni, Hanggi, Jung and Marchesoni, Rev. Mod. Phys. 70, 223 (1998).
3. McNamara and Wiesenfeld, Phys. Rev. A 39, 4854 (1989).

**TODO:** references are listed from standard citations of this literature and
were not verified against the papers themselves, which were not reachable. The
prefactor question in section 8 depends on reference 2.
