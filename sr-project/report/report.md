# Stochastic resonance in a bistable system

Every number in this report is reproduced in `report/numbers.md`, which
`report/extract_numbers.py` generates from `results/`. Figures come from
`figures/`, each drawn by the experiment script named in the caption.

## 1. Introduction

A signal too weak to cross a detector's threshold is invisible, and the usual
instinct is to remove noise. Stochastic resonance, first proposed by Benzi,
Sutera and Vulpiani [1] and reviewed by Gammaitoni et al. [3] and by Wiesenfeld
and Moss [12], is the observation that the opposite can help: in a nonlinear
system with a threshold, adding noise can make a sub-threshold signal easier to
see, up to an optimal noise level, beyond which performance falls again.

This report tests that claim on the standard model, an overdamped particle in a
double-well potential driven by a weak sine. The scope is deliberately narrow.
We measure where the optimum sits, how it moves with the signal, whether it
improves detection rather than just a spectral ratio, and how it compares with
doing nothing and with a linear filter. Section 7 repeats the comparison on
low-light images and on sound. We do not claim stochastic resonance is a good
engineering choice; section 8 argues it usually is not for this problem.

## 2. Theory

The model is

    dx/dt = a x - b x^3 + A cos(2 pi f0 t) + sqrt(2 D) xi(t)

with xi unit white Gaussian noise. The potential U(x) = -a x^2 / 2 + b x^4 / 4
has minima at x_m = +-sqrt(a/b) separated by a barrier dU = a^2 / (4 b). At the
baseline a = b = 1 the barrier is **0.2500** and the static threshold amplitude
is **A_c = 0.3849**. The drive amplitude is **A = 0.3**, which is **0.779** of
A_c, so without noise the particle can never leave its well.

Noise lets it escape at the Kramers rate [6] r_K = a / (sqrt(2) pi) exp(-dU /
D). The two-state model of McNamara and Wiesenfeld [4] gives an output SNR
proportional to (A x_m / D)^2 r_K, which peaks where d/dD of its logarithm
vanishes, at **D = dU/2 = 0.1250**. That value is the prediction this project
tests, and it is independent of the prefactor.

Two caveats, both of which matter later. The two-state model describes only
well-to-well hopping and contains no intra-well motion, so it says nothing
about what happens when hopping is rare. And its prefactor is unverified: see
section 8.

## 3. Method

The Langevin equation is integrated with the Euler-Maruyama scheme [13, 14] at
dt = 0.01, with the noise increment scaled as sqrt(2 D dt); scaling by dt
instead of sqrt(dt) is the classic error and is guarded by a unit test. The
solver is vectorised over noise intensity and trial, so a whole sweep steps as
one array.

Output is decimated to fs = 10 and recorded over an exact integer number of
signal periods, which puts f0 on FFT bin 128 with no leakage. The SNR estimator
averages periodograms across trials, then takes the ratio of the f0 bin to the
median of the surrounding background bins. Error bars are bootstrap [15] 16th
and 84th percentiles over trials.

Two numerical checks back the solver. Halving dt moves the E3 peak not at all,
still **D = 0.1153**, with a maximum deviation of **0.69 dB** across the whole
curve. Switching to Heun, a stochastic Runge-Kutta scheme [14], likewise leaves
the peak at **D = 0.1153** with a maximum deviation of **0.66 dB**. Figure
`e3_numerical_checks` shows all three curves.

## 4. Results: stochastic resonance exists

**Figure `e1_regimes`.** Three regimes in the time domain. At D = 0.02 the
particle hops rarely and out of step with the drive, giving **18.88 dB**. At
D = 0.12 it tracks the drive, giving **28.11 dB**. At D = 1.0 hopping is faster
than the signal and coherence is lost, giving **20.68 dB**. The middle panel is
the effect in one picture.

**Figure `e2_threshold`.** The same effect without any dynamics, the threshold
form of stochastic resonance [9]. A hard threshold at theta = 1 with a sub-
threshold drive A = 0.7 produces nothing at all for the **4** lowest noise
levels, rises to **39.50 dB** at **sigma = 0.621**, and falls to **30.58 dB**
by sigma = 3. Noise alone carries the signal across.

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
tuned instead, as in parameter-tuned stochastic resonance [19], the measured
optimum does not follow the prediction. Sweeping a with b = a at three fixed
input noise levels, the interior peaks sit at **a = 0.394**, **0.281** and
**0.251** against predicted 0.20, 0.40 and 0.80, and at **0.281**, **0.394**
and **0.281** at the weaker amplitude against 0.40, 0.80 and 1.20. The peaks do
not move right as D_in grows. Section 8 explains why, and corrects the
prediction itself.

### 5.1 The effect is not specific to a sine

**Figure `e8_pulses`.** A random sub-threshold pulse train replaces the sine,
scored by peak cross-correlation between input and output rather than power at
one frequency, the aperiodic stochastic resonance measure of Collins, Chow and
Imhoff [20]. The curve has the same two regimes as everything else: an intra-
well maximum of **r = 0.582** at the lowest noise, a dip to **r = 0.166 at D =
0.0403** as hopping begins, then a genuine resonance at **r = 0.423 at D =
0.1230**, decaying to **0.071** by D = 2. That resonance sits **1.6 %** from
the optimum the periodic drive gives, which is the result worth having: the
effect does not depend on the input being periodic.

### 5.2 The hopping really is synchronised

**Figure `e9_residence_times`.** Residence times are how long the particle
stays in one well between committed transitions; their distribution is the
synchronisation test of Gammaitoni et al. [21]. If hopping is locked to the
drive rather than merely frequent, those times cluster at odd multiples of half
the drive period. Measuring the fraction of residences falling near odd half
periods against the fraction near even ones gives **813** at D = 0.04, **3** at
D = 0.125 and **0** at D = 0.6.

The structure is sharpest *below* the SNR optimum, not at it. At D = 0.04 the
particle hops rarely enough to skip whole cycles, and the histogram shows five
clean peaks at T/2 through 9T/2 over two decades of probability density, with a
mean residence of **79.3 time units** across **3205** residences. At the
optimum it hops almost every half period, so nearly all the weight collapses
into the first peak, mean **39.5**. At D = 0.6 the timing is set by noise alone
and the structure is gone, mean **4.5**. This is independent confirmation of
the synchronisation that E1 shows only by eye and E3 only as a spectral ratio,
and it uses no SNR estimator at all.

One caution on that figure. The Kramers time printed on each panel is for the
*undriven* system. At D = 0.04 it is 2301 against a measured mean of 79.3,
because with A = 0.3 the drive periodically lowers the barrier and does most of
the work. The gap is the drive's contribution, not an error.

## 6. Results: detection and baselines

SNR is a spectral ratio over a long record. A detector sees one short record
and must answer yes or no, so section 6 measures that directly with ROC
analysis [16]: 500 runs per hypothesis, statistic the f0 periodogram power of a
single 16-period record.

**Figure `e6_detection`.** For the bistable system the answer is layered. At
the baseline amplitude the problem is trivial, AUC **0.990 to 1.000** across
every noise level, so the sweep runs at a weaker A = 0.05 where the metric
carries information. There, AUC starts at **1.000**, collapses to a dip of
**0.665** at D = 0.0533 as hopping begins, recovers to **0.792** at D = 0.1626,
and declines to **0.567** at D = 2. Detection probability at a 5 % false-alarm
rate does the same more sharply, dipping to **0.082** and recovering to
**0.346**, a gain of **+0.264**. So noise does improve detection, at the
predicted place, but only relative to the dip. Kay [18] gives the general
condition under which added noise can help a detector that is not already
optimal.

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
filter [17], so it is the strongest linear detector available for this problem.

## 7. Extensions: images and sound

Both extensions ask section 6's question on signals that are not a sine: does
added noise help a real signal through a threshold, and does it beat the
obvious alternative? The theory of section 2 does not cover either, so these
results are measured, not predicted.

### 7.1 Low-light images (E10)

**Figure `e10_grating`, the two-dimensional foundation.** Before real images,
E10 first repeats E2 with space in place of time: a 128 x 128 grating
A cos(2 pi (kx x + ky y) / N) at A = 0.7, below theta = 1 everywhere, thresholded
pixel by pixel over 32 frames and scored by the spatial SNR at its wavenumber
from the 2-D FFT. The **5** smallest noise levels give a blank frame at the
-30 dB floor; the pattern peaks at **30.51 dB** at **sigma = 0.690** and fades to
**21.46 dB** by sigma = 3. The closed form for independent pixels,
|sum p e^(-ik.x)|^2 / sum p(1 - p), peaks at the same sigma with **30.37 dB**
and stays within **0.31 dB** of the simulation wherever it is above 0 dB. So
two-dimensional threshold SR behaves exactly as predicted; what the image part of E10
tests is whether it is worth using on pictures.

Each image's brightness channel (V in HSV [22]) is scaled to 12 % of the
original, given mild sensor noise and rounded to 8 bits. The threshold theta is
the brightest pixel of the dark image, so every pixel sits at or below it and a
noiseless threshold returns an almost black frame. Stochastic resonance adds
Gaussian noise of standard deviation sigma, thresholds, averages 50 such frames
and stretches the result to full range. This is the threshold-and-average
scheme of Simonotto et al. [5], with the frame average in place of their
viewers' temporal averaging. The baselines are histogram equalisation [23],
CLAHE [24, 8] and gamma 0.45 [23], all deterministic. Quality is SSIM [7] and
PSNR [23] against the original V channel; histogram equalisation, CLAHE and
both metrics come from scikit-image [11]. Six images were used: two
illustrations by @Sento123450 (`MabelleMona1`, `MabelleMona2`) used with the
artist's permission [26], a binary silhouette (`badapple`, a frame from Anira's
shadow-art video of "Bad Apple!!" [27]), a photograph (`black white river`,
Keith Cooper's Hood Canal photograph from his black-and-white test image [28]),
a colour gradient (`gradiant`, the Granger chart from Reichmann's article [29])
and a ring pattern (`images`, Lazur's public-domain test image [30]). They are
third-party test inputs, used here for study only and not redistributed with
the project. Each image's noise seed comes from its file name, so adding or
removing an image leaves the others' results unchanged.

**Figure `e10_sr_curve`.** A resonance on every image. SSIM rises with sigma to
a maximum at sigma between **0.034 and 0.099**, then falls as noise swamps the
frames. Five images start from almost nothing at the smallest sigma. The
silhouette starts at **0.372**, because a hard threshold already suits a
two-level picture, and it peaks earliest, at sigma = 0.034.

**Figure `e10_image_grid`.** SSIM for every method, best per image in bold:

| Image | Dark input | Hist. eq. | CLAHE | Gamma 0.45 | SR at best sigma |
| --- | --- | --- | --- | --- | --- |
| images | 0.067 | 0.800 | **0.906** | 0.530 | 0.653 |
| badapple | **0.470** | 0.112 | 0.372 | 0.379 | 0.463 |
| MabelleMona1 | 0.185 | 0.408 | 0.520 | **0.591** | 0.195 |
| MabelleMona2 | 0.195 | 0.275 | 0.466 | **0.693** | 0.179 |
| black white river | 0.198 | 0.258 | 0.441 | **0.659** | 0.138 |
| gradiant | 0.222 | 0.303 | 0.414 | **0.593** | 0.138 |

The best baseline beats stochastic resonance on **5 of 6** images. The
exception is the silhouette, where SR at **0.463** beats all three baselines
but not the untouched dark input at **0.470**. SR improves on the dark input
for only **2 of 6** images. By PSNR the ranking is simpler: CLAHE is best on
**6 of 6**, from **17.52 to 23.17 dB**, against **7.92 to 16.99 dB** for SR.

The reason is the same as in section 6. Each output pixel is the fraction of
frames in which v + noise exceeds theta, which tends to Phi((v - theta) /
sigma). That is a sigmoid tone curve, limited to 51 levels and carrying
binomial noise from only 50 frames. Gamma and CLAHE remap brightness without
the noise. A small sigma makes the sigmoid too steep, so only pixels close to
theta get through; a large one flattens it and noise dominates. The optimum
falls where the sigmoid's width is comparable to the dark image's 0 to 0.12
brightness range, an amplitude analogue of the time-scale matching in E4.

**Figure `e10_frames`.** More frames remove the grain but not the tone curve.
Each pixel is a binomial fraction, so its noise should be sqrt(p (1 - p) / N)
and fall as 1/sqrt(N). At each image's best sigma, two independent runs of up
to 400 frames were compared. On log axes the grain falls with slope **-0.501
to -0.497**, within **0.8 %** of sqrt(mean p(1 - p) / N) at every N, and the
50-frame point reproduces the table above exactly. SSIM rises with N: on the
river photograph from **0.138** at 50 frames to **0.334** at 400, against
**0.508** in the limit N -> infinity, where the output is the sigmoid alone.
At 400 frames SR beats the dark input on **6 of 6** images, but the best
baseline on only **1 of 6**, the silhouette. Even the grain-free limit beats
the best baseline on only **2 of 6**, the silhouette and the ring pattern
(**0.924** against CLAHE's **0.906**). So 50 frames understate SR, yet on the
photograph, the gradient and the two illustrations the sigmoid is the wrong
tone curve whatever N is. Sigma was held at its 50-frame optimum; with more
frames a smaller sigma becomes affordable, so these limits are lower bounds.

**Figures `e10_denoise` and `e10_denoise_crops`.** Removing the grain
afterwards. The SR output was denoised by total variation (TV) [31, 32] or
non-local means (NLM) [33], each with its strength set from a blind noise
estimate, the median absolute pixel difference scaled as in [34], never from
the original. Denoising keeps the mean, so the image does not go dark again: on
the river photograph the mean stays at **0.412** while SSIM rises from
**0.138** to **0.760**. Denoised SR beats the best undenoised baseline on
**5 of 6** images. The exception is the ring pattern, where denoising lowers
SR's SSIM (**1 of 6**): the estimator takes the fine rings for noise.

Repeating SR then TV on its own output brightens it only if the threshold
moves to the median. At the maximum every pass is darker, from **0.59x** to
**0.31x** the original brightness over five passes, with SSIM falling from
**0.649** to **0.377** (means over the six images): thresholding a full-range
image at its maximum lets only near-white pixels through. At the median the
second pass lifts brightness to **0.92x** and improves SSIM on **4 of 6**
images, river to **0.799**, but the crops show the gain is partly smoothing:
the tree loses its branches. Later passes add nothing.

The fair test gives the baselines the same denoisers, before or after
enhancing. Then CLAHE with NLM or TV reaches **0.926** to **0.945** on the
photograph, the gradient and the two illustrations, and the best denoised
baseline beats every SR variant on **6 of 6** images. CLAHE's low score
without denoising came from the sensor noise it amplifies, not from its tone
curve. Adding noise and taking it out again cannot supply information that a
deterministic method lacks when the threshold is ours to choose. Where the
threshold is physical, as in a sensor that counts single photons as binary
frames [35], SR followed by denoising would be the natural pipeline.

Approach B, the dynamic stochastic resonance of Chouhan, Jha and Biswas [2],
was not implemented. The
paper could not be retrieved, and the related paper that could be read works in
the DCT domain, which is a different method. This report makes no claim about
approach B.

### 7.2 Sound (E11)

No recording was supplied, so a synthetic three-note melody stands in: 220,
277 and 330 Hz with two harmonics each and a plucked envelope. Placing a WAV at
`data/audio/input.wav` replaces it with no code change. Quality is waveform
fidelity, 10 log10(1 / (1 - rho^2)) with rho the correlation to the clean
sound, so silence scores 0 dB. Every case is also written as a WAV under
`results/audio/` for listening.

**Figure `e11_audio`, top left: threshold array.** N comparators with threshold
1, each adding its own noise, averaged; the sound peaks at 0.8 of the
threshold. Without noise the output is silence. Noise lets the sound through,
peaking at **0.20 dB** for one unit, **2.09 dB** for 16 and **4.65 dB at sigma
= 0.67** for 64. This is dither [9] in the array form of Collins, Chow and
Imhoff [10]: each unit's noise averages out across units while the shared
signal does not. A linear sensor array of the same size, computed analytically,
scores **10.28 dB** at that same sigma and beats the threshold array at every
sigma tested, by at least **+0.06 dB** at the largest. Stochastic resonance
pays only when the threshold cannot be removed.

**Top right: bistable model.** The model resonates near one frequency, so time
is stretched, a frequency rescaling in the spirit of Leng et al. [25], until
the melody's dominant **277.3 Hz** lands on f0 = 0.05, at **35** solver steps
per audio sample. With no noise the output follows the sound at **11.75 dB**,
through intra-well motion alone, which is a linear low-pass response. Any noise
destroys it: the best with noise is **0.05 dB at D = 0.5**. Hopping turns the
output into a random two-level signal.

**Bottom left: pure tone.** The same route on a 277 Hz tone separates two
meanings of "works". SNR at 277 Hz has a shallow maximum of **28.15 dB at
D = 0.15**, near dU/2 = 0.125, against **27.50** at D = 0.03 and **23.78** at
D = 0.8. Waveform fidelity stays between **0.04 and 0.30 dB** at every D. The
tone stands far above the background in its own FFT bin, but most of the
output's power is broadband switching. With four trials per point the bump is
small; the contrast between the two metrics is the result.

So stochastic resonance helps a sub-threshold sound cross a threshold and helps
detect a tone. It does not enhance or reproduce sound, and a linear sensor does
better whenever one is available. **TODO:** a check on input that was already
noisy made the output worse, but that run is not in `results/` and is not
quoted here; move it into `e11_audio.py` before relying on it.

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
PSD, gives pi/2 instead, a factor of two and therefore 3 dB. The difference is
a spectrum convention, not physics, and it does not move the peak. The measured
offset between theory and simulation at the peak is **-2.30 dB**, which
excludes neither: with pi the theory sits 2.30 dB above the measurement, with
pi/2 it would sit 0.70 dB below. **TODO:** settle this against Gammaitoni et
al. [3], which has not yet been consulted. No equation number is cited here
rather than inventing one.

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
4. The extensions agree with the sine. On images and on sound, noise helps a
   signal cross a threshold, but a deterministic alternative, a tone curve or
   a linear sensor, wins wherever one is available: on **5 of 6** images, still
   so at 400 frames, on **6 of 6** once both sides are denoised, and at every
   noise level for sound.

## 9. Conclusion

An overdamped bistable system driven by a sub-threshold sine shows a clear
maximum in output SNR at **D = 0.1153**, within **7.7 %** of the predicted
dU/2 = 0.1250 and within the grid resolution, confirmed by two independent
numerical schemes agreeing to better than **0.7 dB**. The optimum is
independent of signal amplitude and rises with signal frequency, as time-scale
matching requires. Noise measurably improves detection, taking a threshold
detector from chance to certainty and back, but for a sine of known frequency
in white Gaussian noise a linear matched filter remains **13.62 dB** better
than the best stochastic-resonance configuration tested. The same pattern holds
on low-light images and on sound: a resonance in noise exists, but a
deterministic method does better.

### References

1. Benzi, Sutera and Vulpiani, "The mechanism of stochastic resonance",
   J. Phys. A 14, L453 (1981). https://doi.org/10.1088/0305-4470/14/11/006
2. Chouhan, Jha and Biswas, "Enhancement of dark and low-contrast images using
   dynamic stochastic resonance", IET Image Processing 7, 174 (2013).
   https://doi.org/10.1049/iet-ipr.2012.0114
3. Gammaitoni, Hanggi, Jung and Marchesoni, "Stochastic resonance",
   Rev. Mod. Phys. 70, 223 (1998). https://doi.org/10.1103/RevModPhys.70.223
4. McNamara and Wiesenfeld, "Theory of stochastic resonance", Phys. Rev. A 39,
   4854 (1989). https://doi.org/10.1103/PhysRevA.39.4854
5. Simonotto, Riani, Seife, Roberts, Twitty and Moss, "Visual perception of
   stochastic resonance", Phys. Rev. Lett. 78, 1186 (1997).
   https://doi.org/10.1103/PhysRevLett.78.1186
6. Kramers, "Brownian motion in a field of force and the diffusion model of
   chemical reactions", Physica 7, 284 (1940).
   https://doi.org/10.1016/S0031-8914(40)90098-2
7. Wang, Bovik, Sheikh and Simoncelli, "Image quality assessment: from error
   visibility to structural similarity", IEEE Trans. Image Process. 13, 600
   (2004). https://doi.org/10.1109/TIP.2003.819861
8. Zuiderveld, "Contrast limited adaptive histogram equalization", in Graphics
   Gems IV, 474 (1994). https://doi.org/10.1016/B978-0-12-336156-1.50061-6
9. Gammaitoni, "Stochastic resonance and the dithering effect in threshold
   physical systems", Phys. Rev. E 52, 4691 (1995).
   https://doi.org/10.1103/PhysRevE.52.4691
10. Collins, Chow and Imhoff, "Stochastic resonance without tuning", Nature
    376, 236 (1995). https://doi.org/10.1038/376236a0
11. van der Walt et al., "scikit-image: image processing in Python", PeerJ 2,
    e453 (2014). https://doi.org/10.7717/peerj.453
12. Wiesenfeld and Moss, "Stochastic resonance and the benefits of noise: from
    ice ages to crayfish and SQUIDs", Nature 373, 33 (1995).
    https://doi.org/10.1038/373033a0
13. Maruyama, "Continuous Markov processes and stochastic equations", Rend.
    Circ. Mat. Palermo 4, 48 (1955). https://doi.org/10.1007/BF02846028
14. Kloeden and Platen, Numerical Solution of Stochastic Differential
    Equations, Springer (1992). https://doi.org/10.1007/978-3-662-12616-5
15. Efron and Tibshirani, An Introduction to the Bootstrap, Springer (1993).
    https://doi.org/10.1007/978-1-4899-4541-9
16. Fawcett, "An introduction to ROC analysis", Pattern Recognit. Lett. 27, 861
    (2006). https://doi.org/10.1016/j.patrec.2005.10.010
17. Kay, Fundamentals of Statistical Signal Processing, Vol. II: Detection
    Theory, Prentice Hall (1998). ISBN 013504135X
18. Kay, "Can detectability be improved by adding noise?", IEEE Signal
    Process. Lett. 7, 8 (2000). https://doi.org/10.1109/97.809511
19. Xu, Duan, Bao and Li, "Stochastic resonance with tuning system parameters:
    the application of bistable systems in signal processing", Chaos, Solitons
    and Fractals 13, 633 (2002). https://doi.org/10.1016/S0960-0779(00)00266-6
20. Collins, Chow and Imhoff, "Aperiodic stochastic resonance in excitable
    systems", Phys. Rev. E 52, R3321 (1995).
    https://doi.org/10.1103/PhysRevE.52.R3321
21. Gammaitoni, Marchesoni, Menichella-Saetta and Santucci, "Stochastic
    resonance in bistable systems", Phys. Rev. Lett. 62, 349 (1989).
    https://doi.org/10.1103/PhysRevLett.62.349
22. Smith, "Color gamut transform pairs", ACM SIGGRAPH Computer Graphics 12,
    12 (1978). https://doi.org/10.1145/965139.807361
23. Gonzalez and Woods, Digital Image Processing, 4th ed., Pearson (2018).
    ISBN 9780133356724
24. Pizer et al., "Adaptive histogram equalization and its variations",
    Comput. Vis. Graph. Image Process. 39, 355 (1987).
    https://doi.org/10.1016/S0734-189X(87)80186-X
25. Leng, Wang, Guo, Xu and Fan, "Engineering signal processing based on
    bistable stochastic resonance", Mech. Syst. Signal Process. 21, 138 (2007).
    https://doi.org/10.1016/j.ymssp.2005.08.002
26. @Sento123450, two illustrations (`MabelleMona1.jpg`, `MabelleMona2.jpg`),
    https://x.com/Sento123450. Used with the artist's permission, given
    directly to the author, for this project only; not for redistribution and
    not for training any AI system.
27. Anira (あにら), "【東方】Bad Apple!! ＰＶ【影絵】", shadow-art music video,
    Niconico sm8628149, 27 October 2009,
    https://www.nicovideo.jp/watch/sm8628149. Song: "Bad Apple!! feat. nomico",
    Alstroemeria Records, arranged by Masayoshi Minoshima, album *Lovelight*
    (2007). Original composition: ZUN (Team Shanghai Alice), stage 3 theme of
    *Touhou Gensokyo ~ Lotus Land Story* (1998).
28. Keith Cooper, "A Black and White test image", Northlight Images (2005),
    https://www.northlight-images.co.uk/a-black-and-white-test-image/ . The
    Hood Canal photograph (Washington State, September 2004), (c) 2012 Keith
    Cooper; the author states the test image is free to use for non-commercial
    purposes.
29. Michael Reichmann, "How to Create Color Test Charts" (the Granger chart),
    Luminous Landscape,
    https://luminous-landscape.com/how-to-create-color-test-charts/
30. Lazur, "black and white test image", Openclipart (6 December 2016), public
    domain, https://openclipart.org/detail/267976/black-and-white-test-image
31. Rudin, Osher and Fatemi, "Nonlinear total variation based noise removal
    algorithms", Physica D 60, 259 (1992).
    https://doi.org/10.1016/0167-2789(92)90242-F
32. Chambolle, "An algorithm for total variation minimization and
    applications", J. Math. Imaging Vis. 20, 89 (2004).
    https://doi.org/10.1023/B:JMIV.0000011325.36760.1e
33. Buades, Coll and Morel, "A non-local algorithm for image denoising", Proc.
    IEEE CVPR 2005, vol. 2, 60. https://doi.org/10.1109/CVPR.2005.38
34. Donoho and Johnstone, "Ideal spatial adaptation by wavelet shrinkage",
    Biometrika 81, 425 (1994). https://doi.org/10.1093/biomet/81.3.425
35. Fossum, Ma, Masoodian, Anzagira and others, "The quanta image sensor: every
    photon counts", Sensors 16, 1260 (2016). https://doi.org/10.3390/s16081260

References 1 to 5 are the proposal's bibliography, in its numbering; 6 to 25
are the sources of the methods used; 26 to 30 are the test images; 31 to 35
were added for the denoising follow-up. Every DOI
was checked against its Crossref record and the two ISBNs against OpenLibrary
on 2026-09-30 (31 to 35 on 2026-10-01), which confirms authors, title, venue,
volume and pages.
Reference 27 was checked against the Niconico video record and Wikipedia;
references 28 to 30 were located from image-search results and checked against
their own pages the same day. That confirms the citations, not what the works
say: **TODO:** references 1 to 25 and 31 to 35 have not been read in full, and the prefactor
question in section 8 depends on reference 3.
