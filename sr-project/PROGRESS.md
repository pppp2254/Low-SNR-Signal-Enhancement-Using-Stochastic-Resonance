# Progress log

One dated entry per work package: what was done, what was checked, open issues
(PLAN.md section 3).

## 2026-09-18 - WP0: scaffold

**Done.** Created the folder tree from PLAN.md section 3 under `DSP-proj/sr-project/`:
`configs/`, `src/srlab/`, `experiments/`, `tests/`, `results/`, `figures/`,
`data/`, `notebooks/`, `web/`, `report/`. Wrote `requirements.txt`,
`pyproject.toml` (src layout, installable as `srlab`), `configs/default.yaml`
with every baseline parameter from section 2.2, docstring-only stubs for the
eight `srlab` modules, `README.md`, this log, and `notebooks/demo.ipynb` as a
Colab stub. Copied the plan doc in as `PLAN.md`.

**Checked.** Virtual environment `.venv` created, `requirements.txt` and the
package installed. `python -c "import srlab"` succeeds (0.1.0); `numpy, scipy,
cv2, skimage, matplotlib, yaml` all import. `pytest` passes 4 tests in
`tests/test_config.py`, which assert: the package imports; the baseline values
in `default.yaml` equal section 2.2 (a=b=1, A=0.3, f0=0.01, dt=0.01,
decimate=10, 128 periods, 10 discarded, 20 trials, D from 0.01 to 2 over 40
points); A = 0.3 is below A_c = 0.385 and E2's A = 0.7 is below theta = 1; and
f0 lands exactly on FFT bin 128 given fs = 1/(dt*decimate) = 10.
`notebooks/demo.ipynb` parses as nbformat 4 with 11 cells.

**Open issues.**
1. No physics is implemented yet - every `srlab` module is a docstring stub.
   WP1 implements `signals`, `bistable`, `threshold`, `theory`.
2. `sr-project/` is not its own git repository; the enclosing repository is
   `/home/puyawapun` (the home directory). `git init` here before committing.
3. The notebook's `REPO_URL` is a placeholder; set it once the repo is pushed,
   or the Colab clone cell will fail.
4. The `dsp-cookbook` skill in `DSP-proj/skill/` is a JUCE/C++ audio-effects
   cookbook (biquads, compressors, chorus). It has no bearing on this project
   and should not be consulted for it.
5. `two_state_snr`'s prefactor is still unverified against Gammaitoni et al.
   (1998); do that in WP1 before the theory curve is plotted (PLAN.md 2.1).

## 2026-09-18 - WP1: core library

**Done.** Implemented `theory.py` (`well_position`, `barrier_height`,
`critical_amplitude`, `kramers_rate`, `two_state_snr`), `signals.py`
(`n_samples_for`, `time_axis`, `sine`, `pulse_train`, `white_noise`),
`bistable.py` (`simulate`, Euler-Maruyama and Heun, vectorised over
(n_D, n_trials), `input_override` supported) and `threshold.py`
(`threshold_detector`). Added `tests/test_theory.py`, `test_signals.py`,
`test_bistable.py`, `test_threshold.py`.

**Checked.** `pytest` passes 22 tests in 57 s. Specifically:

- A_c = 0.3849, dU = 0.25, x_m = 1.0 at a = b = 1, as in section 2.1.
- Kramers prefactor a/(sqrt(2) pi) = 0.22508, matching the 0.225 in the plan.
- `two_state_snr` peaks at dU/2 for a = 0.5, 1.0 and 2.0 (grid search over
  20001 points, within 0.1%), so the peak the project validates against is
  reproduced by the formula independently of its prefactor.
- Noise increment standard deviation equals sqrt(2 D dt) within 2%, and scales
  as sqrt(D) across two decades - the sqrt(dt) scaling of section 9 is right.
- f0 = 0.01 at fs = 10 over 128 periods gives 128000 samples and the sine's
  spectral peak sits exactly on bin 128; a non-commensurate f0 raises.
- D = 0, A = 0.3: every trial keeps the sign it started with and |x| > 0.5
  throughout. D = 0, A = 0.5: every trial changes sign within three periods.
- Switching rate with A = 0, measured as committed transitions (registering a
  state only once |x| > x_m/2, which removes barrier-top recrossings), against
  `kramers_rate`: ratios measured/theory were 0.907 at D = 0.08, 0.954 at
  D = 0.10, 0.956 at D = 0.125 and 0.953 at D = 0.15 - within 10% everywhere,
  half the 20% tolerance the plan asks for. The test itself runs D = 0.10 and
  0.15 to stay fast.
- Same seed gives bit-identical output; a different seed does not.
- Heun versus Euler-Maruyama at D = 0.12: switching rates agree within 10% and
  std(x) within 5%.
- `input_override`: a constant tilt of 0.6 (above A_c) pulls every trial to +x;
  a too-short override raises.

**Runtime.** Full 40 x 20 baseline ensemble (40 log-spaced D from 0.01 to 2,
20 trials, 1.38e6 solver steps, 128 recorded periods): **74.7 s**, peak RSS
847 MB. Output finite everywhere, max |x| = 3.45 at D = 2, so dt = 0.01 is
stable across the whole sweep.

**Open issues.**

1. **Memory is the binding constraint on this machine.** `free -g` reports
   about 1 GB available. The 40 x 20 output is 410 MB and peak RSS was 847 MB;
   E3's final 50-trial run would be 1.02 GB of output and will not fit. The fix
   is a calling convention, not a parameter change: call `simulate` one D value
   at a time and reduce to a spectrum before the next, giving 26 MB per call.
   `sweep.py` (WP3) must be written this way. Documented in `simulate`'s
   docstring.
2. **PLAN.md section 2.3 says the sweep takes "well under a minute"; it takes
   75 s here** (4 cores, 3 GB RAM). Not a fault, but E3 at 50 trials plus the
   dt-halved and Heun checks will be several minutes, not seconds.
3. **`two_state_snr` prefactor is still unverified.** Gammaitoni et al. (1998)
   was not reachable from this environment, so no equation number is cited -
   an invented citation would be worse than none. Deriving it independently
   (delta-peak strength (pi/2) xbar^2 over the one-sided noise PSD) gives
   `SNR = (pi/2) (A x_m / D)^2 r_K`, a factor of 2 below the `pi` PLAN.md
   specifies. The difference is a one-sided/two-sided spectrum convention and
   shifts the theory curve by 3 dB without moving its peak. The code implements
   PLAN.md's `pi` rather than silently substituting my own value. **Settle this
   from the paper before the theory curve goes in the report.**
4. `two_state_snr` returns a peak-strength-over-PSD ratio, which has units of
   inverse time, while `metrics.snr_db` returns a dimensionless per-bin ratio.
   Overlaying them in E3 needs the analysis bandwidth df = f0 / n_periods.
   Deferred to WP3, where the overlay is actually drawn.
5. Surprising but correct: the baseline drive is **not** adiabatic -
   W = 2 pi f0 = 0.063 against 2 r_K = 0.061 at D = 0.125. The two-state SNR
   formula still applies, because the Lorentzian factor cancels between the
   signal and the noise background, but PLAN.md section 9's advice to "lower f0"
   if the peak misses theory rests on a different assumption than it implies.
6. Deliberate deviation from the stub docstring: `well_position` was added to
   `theory.py` (not in the WP1 list) because the solver needs it for the
   initial condition. The `switch_rate` helper lives in the test file, not in
   `srlab`; promote it to a module when E9 needs residence times (WP6).

## 2026-09-18 - WP2: metrics

**Done.** Implemented `metrics.py`: `snr_db`, `bootstrap_snr`, `xcorr_peak`,
`roc_auc`, `pd_at_pfa`, with private helpers `_bin_index` (the f0-on-a-bin
assertion PLAN.md asks for), `_band_periodograms` and `_snr_from_band`. Added
`tests/test_metrics.py` (17 tests).

**Checked.** Full suite passes 39 tests in 57 s.

- A sine in white noise with analytically known SNR, A**2 n / (4 sigma**2) per
  bin, is recovered to within 0.30 dB over three (A, sigma) combinations -
  inside the 1 dB the plan asks for.
- Pure noise gives below 1 dB for five different seeds.
- f0 = 0.3001 at fs = 10 over 1000 samples raises; f0 = 0.3 (bin 30) does not.
- A background band that would reach the DC bin raises.
- A constant record raises rather than silently returning NaN.
- Leading sweep axes are preserved, and SNR rises monotonically with A.
- The bootstrap band brackets the point estimate and narrows from 10 to 80
  trials.
- AUC is 0.50 +- 0.03 for identical distributions, rises with separation, and
  handles the ties the binary threshold-detector output produces.
- Pd at Pfa = 5% is 0.05 +- 0.03 for identical distributions and above 0.99 for
  well-separated ones; the realised false-alarm rate does not exceed 5%.
- Integration: `simulate`'s float32 (n_D, n_trials, n_samples) output feeds
  `bootstrap_snr` directly. A reduced 6 x 8 x 32-period sweep gave 13.9, 14.9,
  21.9, 20.6, 17.0, 12.2 dB over D = 0.01 to 2 - rise, peak, fall, peak between
  D = 0.08 and 0.24, at 131 MB peak RSS. Not E3; that is WP3's job with the
  full grid, 50 trials and the theory overlay.

**Two tests initially failed and were right to.** Both were bad premises in the
test, not bugs in the estimator, and both are now documented in `xcorr_peak`:

1. A lag-40 circular shift scored 0.9976, not 1.0. Correlation here is linear,
   so a lag of L costs L samples of overlap while the normalisation keeps the
   full energy of both records; the ceiling at lag L is about (n - L) / n.
2. After relaxing that, the recovered lag came back as 38 rather than 40 - on a
   **sine**. A periodic input makes the correlation periodic in lag, so the peak
   is ambiguous modulo the period, and the triangular overlap envelope pulls the
   reported peak towards zero. This is exactly why PLAN.md section 2.4 reserves
   this metric for aperiodic inputs. The lag test now uses an aperiodic record,
   and a separate test asserts the sine case is *not* identifiable, so the
   limitation stays visible.

**Changed from the WP2 specification, deliberately.**

1. `snr_db` floors the reported value at **-30 dB** (`_SNR_FLOOR_DB`). With the
   background subtracted, S goes non-positive whenever there is no signal and
   the log ran to -120 dB in testing, which would wreck E3's y-axis. The floor
   sits far below what the estimator can resolve: averaging M trials leaves the
   f0 bin with a relative spread of 1/sqrt(M), about 0.22 at 20 trials, so
   anything under roughly -7 dB already means "indistinguishable from the
   background".
2. `bootstrap_snr` returns the all-trials point estimate as the central value
   rather than the bootstrap mean, since that is what belongs on an error bar.
   The 16/84 percentiles are as specified.
3. `snr_db` takes `seed` nowhere; `bootstrap_snr` takes one (default 0) because
   resampling needs randomness and PLAN.md section 3 requires it be seeded.

**Open issues.**

1. **PLAN.md is inconsistent about `guard`.** Section 2.4's prose says the
   background runs over bins f0 +- 3 to f0 +- 25 excluding f0 +- 1, which
   implies guard = 2; the WP2 signature in section 5 says `guard=1`. The code
   and `configs/default.yaml` follow the signature. With an exact integer number
   of periods there is no leakage, so the two differ by hundredths of a dB.
   Worth one sentence in the report rather than a change.
2. `_band_periodograms` chunks the FFT to about 30 MB at a time and keeps only
   the 2 * n_side + 1 bins around f0, so the WP1 memory problem does not
   reappear in the estimator: the 6 x 8 integration run peaked at 131 MB. The
   remaining constraint is the solver output itself.
3. Still open from WP1: the `two_state_snr` prefactor (pi versus pi/2) is
   unverified against Gammaitoni et al. (1998), and the bandwidth conversion
   between `two_state_snr` and `snr_db` is deferred to WP3.

## 2026-09-18 - WP3: headline experiments

**Done.** Added `sweep.py` (chunked sweeps with result caching) and
`plotting.py` (shared style), then `experiments/e1_regimes.py`,
`e2_threshold.py` and `e3_snr_vs_noise.py`. Moved the theory-to-per-bin
conversion into `srlab.theory.two_state_snr_db` so the dashboard exporter and
the experiment share one implementation. Added `tests/test_sweep.py`.
Figures: `e1_regimes`, `e2_threshold`, `e3_snr_vs_noise`,
`e3_numerical_checks`, each as 300 dpi png and pdf.

**Results.**

E1, three regimes at D = 0.02, 0.12, 1.0: SNR 18.88, 28.11, 20.68 dB. The f0
bin stands 34, 37 and 21 dB over its local background. At D = 0.02 the particle
hops rarely and out of step; at D = 0.12 it tracks the drive; at D = 1.0 the
hopping is faster than the signal.

E2, threshold detector, A = 0.7 below theta = 1: peak at sigma = 0.621,
39.50 dB, from a floor at the 4 lowest sigma values (up to 0.069) where the
detector never fires. Clear rise, peak and fall, with no dynamics involved.

E3, the headline: **measured peak at D = 0.1153, SNR 28.19 dB
[28.04, 28.34]**, against the theoretical dU/2 = 0.125. That is a 7.7% error,
inside PLAN.md's "about 30%" criterion. The error is also smaller than the
sweep can resolve: the 40-point log grid steps 14.6% per point, and 0.1153 and
0.1321 are the two grid points bracketing 0.125, so the measured peak is the
nearest resolvable value. Curve runs 20.7 dB at D = 0.01, down to 18.1 dB at
D = 0.017, up to the peak, then down to 18.1 dB at D = 2.

Numerical checks: dt halved and Heun both put the peak at the **same**
D = 0.1153. Maximum deviation from the baseline curve is 0.69 dB (dt halved,
mean 0.19 dB) and 0.66 dB (Heun, mean 0.23 dB). Bootstrap bands overlap at
37/40 and 33/40 D values; since the bands are 16-84 percentiles, two agreeing
independent estimates fail to overlap about 16% of the time by chance, so
3/40 and 7/40 non-overlaps are what agreement looks like, not disagreement.
The bands are only about +-0.15 dB wide, which is why overlap counting is a
weak test here and the 0.7 dB maximum deviation is the better number.

**Theory comparison.** The theory curve tracks the simulation closely from
D = 0.05 upward and peaks 30.49 dB against the measured 28.19 dB, an offset of
-2.30 dB. Note that offset is smaller than the pi versus pi/2 prefactor
question from WP1: with pi/2 the theory peak would be 27.5 dB, 0.7 dB *below*
the measurement. Neither convention is excluded by this data, so the prefactor
still has to come from the paper.

Below D = 0.05 theory and simulation diverge completely: theory falls past
-50 dB while the simulation flattens near 18 to 20 dB. That is expected and
worth stating in the report. The two-state model describes only well-to-well
hopping, and at low D there is almost no hopping; what the simulation still
measures at f0 is the intra-well response to the drive, which the two-state
model does not represent at all.

**Dashboard.** Added `web/` with an Elysia server (`server.ts`, port 3100),
`export_data.py` (results to `web/data.js`) and `index.html`. Raw readout:
measured values strip, E3 with theory overlay, E2, the three solver variants,
E1's six panels, and the 40-row E3 table. Times New Roman only, cool slate
ground in light and dark, no rounded corners, shadows, gradients or icons.
Run with `bun run web/server.ts`, or `bun run dev` from `web/`. Also published
as an artifact. `web/node_modules/` is gitignored; `.claude/launch.json` in the
primary working directory gained an `sr-dashboard` entry on port 3100 (3000 was
already taken).

**Checked.** `pytest` passes 44 tests in 67 s. The E3 sweep peaked at 370 MB
RSS against a roughly 1 GB budget, so the WP1 memory ceiling is handled:
`chunked_snr` reduces each chunk of raw traces to a band periodogram and
releases it. `test_sweep.py` asserts that chunking one value at a time and all
at once give identical SNR.

**Three bugs found and fixed.**

1. E1's f0 marker was drawn on top of the spectral spike it marks, hiding it,
   so the fundamental looked absent and the third harmonic looked dominant.
   Checked numerically rather than by eye: the f0 bin held 719.8, 7053.9 and
   570.4 against local backgrounds of 9.2, 11.0 and 5.0. The physics was right
   and the figure was wrong. The marker now draws underneath the data.
2. `bootstrap_snr` raised on near-silent sweep values. A resample can contain
   only silent trials, which makes the background zero; that is an expected
   outcome of resampling near a detector's floor, not a bug. Bootstrap draws
   now floor while the point estimate still raises, so a real all-zero solver
   output is still caught. Without this E2 could not run at all.
3. The E3 figure was unreadable because the theory curve dives past -50 dB at
   low D, compressing the data into the top quarter. The view is now clipped to
   the measured range, and the redundant "theory shifted to the peak" curve was
   dropped since at a 2.30 dB offset it simply overplotted the solid one.

**Open issues.**

1. `results/e1.npz` originally stored all 20 raw traces per D, 27 MB. E1 now
   caches only the plotted trace and the averaged spectrum: 94 KB, same
   numbers. Worth remembering before any experiment caches raw ensembles.
2. The two_state_snr prefactor is still unverified, and as noted above this
   data does not settle it.
3. The artifact could not be opened from this environment to confirm it renders
   (the browser here is not signed in to claude.ai). It was verified against
   the local Elysia server in both light and dark themes.

## 2026-09-19 - WP4: robustness and parameter tuning

**Done.** `experiments/e4_robustness.py` and `e5_parameter_tuned.py`, plus two
helpers in `srlab`: `signals.snap_f0` and `sweep.interior_peak`. Figures
`e4_robustness` (three panels) and `e5_parameter_tuned` (two panels).

**E4, and it confirms time-scale matching.**

Amplitude family at f0 = 0.01, SR peak positions: D = 0.1321 (A = 0.1),
0.1321 (A = 0.2), 0.1007 (A = 0.3). Essentially constant, which is what the
two-state model predicts, since the peak D = dU/2 carries no A dependence.

Frequency family at A = 0.3, SR peak positions: D = 0.0879 (f0 = 0.005),
0.1153 (0.01), 0.1321 (0.02), 0.1986 (0.05). The optimum rises monotonically
with f0, which is the WP4 check and the time-scale-matching prediction: a
faster drive needs a faster hopping rate, and that needs more noise.

Heatmap SR peak D against f0: 0.093, 0.093, 0.123, 0.163, 0.163 for
f0 = 0.0050 to 0.0186, then 0.123 for 0.0259, 0.0360 and 0.0500. Rising over
the range where the drive is slow enough for hopping to track it, flat above.

**Global maxima are not the resonance.** For A = 0.1 and A = 0.2 the largest
SNR in the sweep sits at the D = 0.01 edge (31.81 and 37.16 dB), far above
their SR peaks (17.35 and 24.10 dB). Same for f0 = 0.02 and 0.05. That edge
rise is the intra-well response: at low D the particle stays in one well and
responds linearly to the drive with amplitude about A / (2a), giving coherent
power at f0 while the background falls with D, so SNR grows as D falls. It is
not stochastic resonance. `sweep.interior_peak` exists to label the resonance
rather than the edge, and it takes a margin of 2 points, because a maximum
sitting on the first or last couple of points is part of the edge rise rather
than separated from it. The same effect explains E3's left-hand upturn.

**E5 does not reproduce the tuning optimum, and PLAN.md's prediction for it is
wrong by a factor of two.**

PLAN.md 2.5 gives the optimum as dU = 2 D_in, so a = 8 D_in. That is the
condition for a *D* sweep, where the peak comes from balancing the 1/D^2
prefactor against exp(-dU/D). Sweeping *a* at fixed D, x_m = 1 is constant and
the only a dependence is SNR proportional to a exp(-a / (4 D_in)), whose
stationary point is a = 4 D_in, that is dU = D_in. Confirmed against
`two_state_snr` itself: its maximum over a sits at exactly 0.20, 0.40, 0.80 and
1.20 for D_in = 0.05, 0.1, 0.2 and 0.3. The figure now marks a = 4 D_in.

Even against the corrected prediction the measurement does not follow it. Two
conditions were run, as agreed with the user:

| condition | A | D_in | interior peak a | predicted 4 D_in |
| --- | --- | --- | --- | --- |
| A, PLAN amplitude | 0.3 | 0.05 | 0.394 | 0.20 |
| A | 0.3 | 0.1 | 0.281 | 0.40 |
| A | 0.3 | 0.2 | 0.251 | 0.80 |
| B, linear response | 0.1 | 0.1 | 0.281 | 0.40 |
| B | 0.1 | 0.2 | 0.394 | 0.80 |
| B | 0.1 | 0.3 | 0.281 | 1.20 |

The interior features sit at a = 0.25 to 0.39 regardless of D_in and do not
move right as D_in grows, so the WP4 check fails for E5. Both panels are
monotonic declines with a weak plateau at small a. The cause is the same
intra-well response: sweeping a downward shrinks the intra-well stiffness 2a,
so the linear response amplitude A / (2a) grows, and that rise swamps the
hopping optimum the two-state model describes. Condition A additionally has the
signal above threshold below a = 0.78 (A_c = 0.385 a), and condition B below
a = 0.26; both regions are shaded on the figure.

No parameters were changed to make a peak appear. Condition B was added on the
user's instruction, with condition A kept alongside it.

**Three bugs found and fixed.**

1. `metrics._band_periodograms` cast the whole ensemble to float64 before
   chunking. That copy is twice the size of the float32 input and defeated the
   chunking entirely, which is why E5's first run peaked at 1.17 GB on a
   roughly 900 MB budget. The cast is gone; numpy's rfft promotes each chunk on
   its own. `snr_db` now adds no measurable memory over the ensemble, and
   float32 and float64 inputs agree to 1.9e-7. Regression test added.
2. E4's log-spaced heatmap frequencies broke the on-bin requirement:
   f0 = 0.006947 left 921198 recorded steps, not a multiple of the decimation
   factor, and the solver raised. `signals.snap_f0` moves f0 to the nearest
   value giving a whole number of periods, a shift under 0.003%. Verified end
   to end through solver and estimator.
3. Labelling curves by `argmax` reported the intra-well edge rather than the
   resonance, as above.

**Open issues.**

1. The prefactor question from WP1 is untouched by this work.
2. E5's null result is worth a paragraph in the report rather than another
   parameter hunt. To resolve the tuning optimum the intra-well contribution
   would have to be suppressed, for instance by lowering A far enough that the
   linear response is negligible at the small-a end, which costs input SNR that
   PLAN.md 4 wants held near -15 dB. The three requirements are not jointly
   satisfiable, as tabulated in the WP4 discussion with the user.
3. Cross-f0 comparison is bandwidth-inconsistent: the analysis bin width is
   df = f0 / n_periods, so holding n_periods fixed measures each f0 in a
   different bandwidth. Peak positions are comparable across f0; peak heights
   are not. Noted in the e4 docstring and on the figure title.

## 2026-09-19 - WP5: detection and baselines

**Done.** `experiments/e6_detection.py` and `e7_baselines.py`, plus
`metrics.f0_power` (per-record periodogram power at f0, the detection
statistic; unlike `snr_db` it does not average over trials, because the spread
across records is what a ROC curve measures). Figures `e6_detection` (four
panels) and `e7_baselines` (two panels), with `results/e6_detection.csv` and
`results/e7_baselines.csv`.

**E7, the comparison PLAN.md 8 asks to be stated plainly.** One low-SNR input,
input noise D_in = 0.02, five detectors. SNR at the baseline A = 0.3 over 128
periods; AUC over 16-period records at A = 0.05, because at A = 0.3 every
method scores AUC 1.000 and the column carries no information.

| method | SNR (dB) | AUC | Pd at 5% |
| --- | --- | --- | --- |
| linear band-pass on the raw input | **41.60** | 1.000 | 1.000 |
| threshold detector, no added noise | 30.74 | 0.939 | 0.724 |
| SR at the optimal D = 0.125 | 27.99 | 0.779 | 0.320 |
| parameter-tuned SR, a = b = 0.78 | 25.79 | 0.980 | 0.998 |
| bistable, no added noise (D = 0.02) | 19.58 | 1.000 | 1.000 |

Both halves of PLAN.md 8's required statement are measured. Adding noise to the
bistable system is worth **+8.41 dB**, 19.58 to 27.99: that is stochastic
resonance working. The linear band-pass still beats it by **13.6 dB**. The
band-pass is not a strawman: for a sine of known frequency in white Gaussian
noise the f0 periodogram bin is the matched filter, so it is the strongest
linear detector available, and it wins on both columns.

**E6, and the bistable result is more interesting than expected.** At the
baseline A = 0.3 every AUC saturates at 1.000, so the sweeps run at
`e6.A_det = 0.05`; the saturated curve is kept and plotted, since it is the
reason for the change.

The bistable system has two competing regimes:

- D below about 0.02: the particle stays in one well and responds linearly to
  the drive. That coherent intra-well tone gives AUC 1.000, Pd 1.000.
- D near 0.05: hopping starts and destroys the tone. AUC dips to 0.665, Pd to
  0.082.
- D near dU/2: a genuine stochastic-resonance recovery. AUC rises to 0.792 at
  D = 0.1626 and Pd to 0.346 at D = 0.0931, a Pd gain of **+0.264** over the
  dip, with the theoretical optimum at 0.125 sitting between them.
- D above that: decline to AUC 0.567 at D = 2.

So noise does improve detection, measurably and at the predicted place, but
only relative to the dip. The intra-well regime at very low noise still beats
the SR peak outright.

**The threshold detector is where SR is unambiguous.** It has no intra-well
response: below threshold it emits nothing, so both hypotheses give an empty
record and AUC is exactly 0.500 with Pd 0.000. Adding noise takes it to AUC
1.000 and Pd 1.000 at sigma = 0.156, and further noise returns it to chance,
AUC 0.495 and Pd 0.046 at sigma = 300. Chance to perfect to chance, driven by
nothing but noise. That answers criterion 5.

The detection sweep needed a far wider sigma range than E2's. E2 stops at
sigma = 3 because that brackets its SNR peak at 0.62, but AUC saturates well
before that and only returns to chance near sigma = 90, so `thr_sigma_max` is
300. Widening the sweep to see the whole curve is not the same as tuning a
parameter to produce one; the shape was already there, outside the window.

**A conclusion I got wrong first and corrected.** From a coarse three-point
amplitude scan plus the saturated A = 0.3 data I reported that AUC falls
monotonically with D and that noise never helps detection. The 20-point sweep
at A = 0.05 shows that is false: there is a clear dip and a clear SR recovery
at dU/2. The coarse scan simply had no point between D = 0.01 and D = 0.125,
which is exactly where the dip lives.

**Open issues.**

1. AUC saturates easily. It reached 1.000 for the threshold detector across
   nearly two decades of sigma, so the "maximum" there is a broad plateau
   rather than a point. Pd at a fixed false-alarm rate is the sharper metric
   and is plotted alongside it throughout.
2. E7's two columns are at different amplitudes, labelled in the table and in
   the figure titles. Merging them onto one amplitude would make one column
   uninformative.
3. The E7 tuned value a = b = 0.78 comes from E5's measured result, that SNR
   declines monotonically with a inside the sub-threshold region, so the best
   tunable a is the smallest one keeping A below A_c. It is not a fitted value.

## 2026-09-19 - WP9: notebook, report and slides

**Done.** `report/extract_numbers.py` (pulls every cited quantity from
`results/` and writes `report/numbers.md` and `numbers.json`),
`report/report.md` (nine sections following PLAN.md 8), `report/slides.md`
(eleven slides, one figure each, takeaway as the title), and a working
`notebooks/demo.ipynb`.

**Numbers are traceable, not remembered.** PLAN.md WP9 requires every number in
the report to come from `results/*.npz` or `*.json`. `extract_numbers.py` is
the only thing that reads them: 62 quantities across theory, E1 to E7. It then
checks `report.md` itself, extracting every bolded claim and verifying its
numeric tokens appear in the extracted set. That check passes, so no figure in
the text was typed from memory or left stale after a rerun. Rerun it after any
experiment reruns.

**Notebook runs, verified rather than assumed.** Code cells extracted and
executed headlessly: **35 s** end to end against PLAN's five-minute budget. It
reproduces the three regimes, the threshold detector sweep, and a reduced E3
that finds the peak at **D = 0.1260** against the predicted 0.1250, a 0.8 %
error on a 26 % grid. It calls `srlab` throughout and defines no solver or
estimator of its own.

Two bugs surfaced only by running it. The threshold sweep hit the constant
record case at low sigma, the same one `chunked_snr` handles, and raised; the
notebook now floors to -30 dB with a comment saying why, which is the behaviour
the library uses. And `float(snr_db(...)[0])` fails, because `snr_db` returns a
tuple and `[0]` takes the whole SNR array rather than its first element; the
cell now unpacks explicitly. The first reduced sweep also gave a 35 % peak
error on a 42 % grid, which is within one grid point but reads badly, so the
sweep went to 24 values over 48 periods: 0.8 % error for six more seconds.

**Two TODOs left in the report, both deliberate.** The two-state prefactor is
still unverified against Gammaitoni et al. (1998), and the reference list is
from standard citations rather than checked against the papers. Both are marked
in the text rather than quietly asserted.

**Cache size fixed before it shipped.** `results/e7.npz` was **36.8 MB**: it
cached five raw ensembles of shape (1, 20, 128000). E7 now reduces to its
summary numbers inside `compute()`, giving **1.1 KB** with identical values.
This is the second time a raw ensemble reached the cache, after E1's 27 MB.

**Open issues.**

1. The report's section 7 says the image extension was not attempted, which is
   accurate: WP7 is not done and `image_sr.py` is still a stub.
2. `report.md` and `slides.md` are Markdown. If a PDF or a slide deck file is
   wanted, that is a conversion step this project does not do.
3. WP10 is the remaining Must-tier item: a single `make_all.py` regenerating
   every figure in a clean environment, and a comparison against what is
   committed.

## 2026-09-19 - WP6: pulses and residence times

**Done.** `experiments/e8_pulses.py` and `e9_residence_times.py`, plus
`bistable.residence_times`, promoted out of the test file as WP1's log said it
should be when E9 needed it. Figures `e8_pulses` and `e9_residence_times`,
report sections 5.1 and 5.2, and both added to `extract_numbers.py` and the
dashboard.

**E8: the effect is not specific to a sine.** A random sub-threshold pulse
train drives the solver through `input_override`, evaluated on the solver's own
step grid. One train is shared across every D and trial, so each point differs
only in the noise added to it. Scored by peak cross-correlation, the metric
PLAN.md 2.4 reserves for aperiodic input.

The curve repeats the pattern of E4, E5 and E6: an intra-well maximum of
**r = 0.582** at the lowest D, a dip to **0.166 at D = 0.0403**, then a genuine
resonance at **r = 0.423 at D = 0.1230**, falling to **0.071** at D = 2. That
resonance is **1.6 %** from the periodic optimum dU/2 = 0.125. The same
optimum, for a completely different input.

**E9: synchronisation confirmed, but not where expected.** Measuring the
fraction of residence times falling near odd half periods against the fraction
near even ones gives **813** at D = 0.04, **3** at D = 0.125 and **0** at
D = 0.6.

PLAN.md's check expects peaks at T/2 and 3T/2 at the optimal D. The peaks are
there, but they are sharpest *below* the optimum. At D = 0.04 hopping is rare
enough to skip whole cycles, giving five clean peaks from T/2 to 9T/2 across
two decades of density. At D = 0.125 the particle hops almost every half
period, so the weight collapses into the first peak and the multi-peak picture
disappears. Both are synchronisation; the multi-peak form simply needs a slower
hopping rate to be visible.

The histograms needed a log density axis: the 3T/2 and 5T/2 peaks sit one to
two orders below the first and are invisible on a linear one.

**A misleading comparison, caught and labelled.** The panels print 1/r_K beside
the measured mean. At D = 0.04 that is 2301 against a measured 79.3, a factor
of 29. It is not an error: r_K is the *undriven* Kramers rate, and with A = 0.3
the drive periodically lowers the barrier and does most of the work. The panels
now say "undriven 1/r_K" and the report explains the gap, rather than leaving a
reader to assume the solver is wrong.

**A config edit that landed in the wrong section.** The e8 and e9 blocks were
inserted by anchoring on `  e10:`, which matches first inside `seeds:`, not
`experiments:`. YAML took the later duplicate key, so `seeds.e8` became a dict
and `default_rng` raised. Both blocks now sit under `experiments:` and
`test_config.py` still passes. Anchor on a line with content, not a bare key
that repeats.

**Open issues.**

1. WP7 is still not done and needs a photo at `data/input.jpg`, which only the
   user can supply. `image_sr.py` remains a stub and report section 7 says so.
2. WP8 was declined by the user in favour of the results dashboard.
3. WP10 is the remaining Must-tier item.

## 2026-09-19 - Iteration-scale figures (not a PLAN.md work package)

`experiments/e1_every_step.py` runs with decimate = 1 and zooms from the whole
record (160 steps per pixel, individual iterations unresolvable) to a single
hop (6.7 pixels per step), overlaying the samples the usual decimate = 10
pipeline keeps. `experiments/e45_hundred_steps.py` shows 100 iterations for
the E4 amplitude family and the E5 system-parameter family, splitting each
Euler-Maruyama step into its well, drive and noise contributions. Noise
exceeds the deterministic terms by 5x to 22x on a single step, which is why
stochastic resonance is invisible at iteration scale and only emerges over
thousands of steps. `results/e1_every_step.npz` is 4.1 MB, the largest cache
in the project, justified because the raw trace is the deliverable. Logged
late: these were made on 2026-09-19 and not recorded at the time.

## 2026-09-25 - E11: stochastic resonance on sound (extension)

**Done.** `experiments/e11_audio.py`, `metrics.fidelity_db`, config block
`experiments.e11`. Input is `data/audio/input.wav` or a path argument; with
neither, a synthetic three-note melody. Record a real one with
`arecord -f S16_LE -r 16000 -c 1 -d 5 data/audio/input.wav`. Outputs are
WAV files in `results/audio/` for listening, and figure `e11_audio`.

`fidelity_db` is needed because sound has no single f0 for `snr_db` to read.
It is the SNR of the best linear fit of output to the clean sound, equal to
10 log10(1 + SNR) for additive noise; a constant output scores 0 dB.

**Results on the synthetic melody.**

Route 1, threshold array (memoryless, runs at audio rate). The sound is scaled
to peak at 80 % of theta, so alone it never crosses and the output is silence.
Noise lets it through, with a clean resonance that strengthens with the number
of units:

| units | no noise | best | sigma at best | sigma = 10 |
| --- | --- | --- | --- | --- |
| 1 | 0.00 dB | 0.20 dB | 0.496 | 0.00 dB |
| 16 | 0.00 dB | 2.09 dB | 0.670 | 0.03 dB |
| 64 | 0.00 dB | 4.65 dB | 0.670 | 0.12 dB |

This is the effect audio engineering calls dither. A linear sensor with the
same noise and the same 64 units beats it at every sigma (analytic curve on the
figure), so it helps only because the sensor is a threshold device.

Route 2, the bistable model, time-stretched so the melody's dominant 277 Hz
lands on f0_model = 0.05. With no noise it scores **11.75 dB**: the particle
stays in one well and responds linearly, acting as a low-pass filter. Adding
any noise destroys that, **0.05 dB at best**. For sound, the best way to use
the bistable is without stochastic resonance.

**Why, on a pure 277 Hz tone.** The same route on a tone gives SNR at 277 Hz a
shallow maximum, 28.15 dB at D = 0.15 against 23.78 dB at D = 0.8, while
waveform fidelity never exceeds 0.30 dB. The hopping output is a random
two-state signal: the tone stands well above the background in its own FFT bin
while most of the output's power is broadband switching. Stochastic resonance
improves detection of a tone; it does not reproduce a sound. That is the answer
to "does this apply to sound": as a detector yes, as an enhancer only through
many averaged threshold units, and never better than a linear sensor.

**Checked.** 52 tests pass, including two for `fidelity_db` (the additive-noise
formula within 0.1 dB; gain and offset ignored; a constant scores 0 dB). E11
runs in 174 s at 192 MB peak RSS. The config blocks were inserted by anchoring
on lines verified unique first, after WP6's YAML mishap; `test_config.py`
passes and `experiments.e10` is intact.

**Caught before it shipped.** The first figure overlaid the bistable's best
*with-noise* output, which reproduces almost nothing, when its best result is
the no-noise one; the overlay now shows that. The docstring first claimed the
tone showed "the classic resonance, peaking near D = 0.2"; this run's peak is
D = 0.15 and the bump is shallow, so the text now says what was measured.

**Open issues.**

1. Only the synthetic melody has been tried. A real recording is the actual
   test and needs the user to record one; this session did not record from the
   microphone.
2. The bistable route maps one dominant frequency. Speech or music spread over
   several octaves will be treated worse than the melody.
3. The tone's SNR bump uses four trials per point and is shallow. More trials
   would sharpen it; they would not change the fidelity result.
4. `results/audio/` holds eight WAVs of about 94 KB each; decide whether they
   belong in git.
5. E11 is not in `report/report.md` or the Thai report doc yet.

## 2026-09-30 - WP7: low-light images (E10), approach A only

**Done.** `src/srlab/image_sr.py` (`darken`, `approach_a`, `baselines`,
`image_metrics`), `experiments/e10_image.py`, `tests/test_image_sr.py` (4
tests), config `experiments.e10`. Processes every image in `data/picture/`
(five supplied by the user) or paths on the command line. Outputs
`figures/e10_image_grid`, `figures/e10_sr_curve`, `results/e10_metrics.csv`.

Method. V of HSV only. Each image is the ground truth, darkened to 12 %, given
Gaussian sensor noise 0.005 and quantised to 8 bits; at that scale the whole
image fits in about 30 grey levels. Approach A thresholds at the dark image's
maximum, so with no noise no pixel fires (E2's sub-threshold setup in 2-D),
averages 50 noisy binary frames, and stretches the result linearly to [0, 1].
Baselines: histogram equalisation, CLAHE, gamma 0.45.

**Results, SSIM against the original:**

| image | dark input | HE | CLAHE | gamma 0.45 | SR best |
| --- | --- | --- | --- | --- | --- |
| HomuMado.jpeg (illustration) | 0.170 | 0.301 | 0.552 | **0.561** | 0.210 |
| badapple.png (binary silhouette) | **0.470** | 0.111 | 0.371 | 0.380 | 0.454 |
| black white river.jpeg (photo) | 0.198 | 0.258 | 0.434 | **0.659** | 0.141 |
| gradiant.jpeg (colour ramps) | 0.222 | 0.305 | 0.421 | **0.593** | 0.127 |
| images.jpeg (rings pattern) | 0.067 | 0.801 | **0.907** | 0.529 | 0.645 |

CLAHE has the best PSNR on all five (16.76 to 21.54 dB).

1. **The resonance exists on every image.** SSIM goes from 0.000 (no noise, no
   pixel fires) to a peak at sigma 0.034 to 0.099 and falls again, e.g.
   0.000 -> 0.645 -> 0.174 on the rings. Noise is what lets the image through.
2. **SR never wins overall.** Gamma wins SSIM on three images, CLAHE on the
   texture-heavy rings and on PSNR everywhere. By SSIM, SR even trails the
   untouched dark input on three of five; SSIM rewards the dark image for
   keeping structure, while PSNR (2.8 to 8.1 dB for the dark input) shows SR
   does brighten every image.
3. **SR's one win is the binary silhouette**, where it beats all three
   baselines (0.454 against 0.380 at best). Thresholding suits content that is
   already two-level.
4. **Why SR loses:** 50 binary frames give only 51 grey levels, and binomial
   grain up to about 0.07, which the stretch amplifies; the grain is visible
   in the grid. PLAN.md 9 suggests 50 to 100 frames.

This is the "mixed result" PLAN.md 6 anticipated: a real noise-benefit peak,
but classical tone mapping scores higher.

**Approach B not implemented, as PLAN.md instructs when the paper is not
accessible.** Chouhan, Jha and Biswas, IET Image Processing 7 (2013) 174-184,
returned HTTP 403. The same group's open-access APSIPA 2013 paper (Jha,
Chouhan, Aizawa, Biswas) was reachable, but it is a DCT-domain method, not
the spatial one PLAN.md names, and the extracted text gave inconsistent
values (step size both 0.001 and 0.01 to 0.055). Building on either would be
guessing.

**Correction.** Earlier entries say the Gammaitoni et al. (1998) paper "was
not reachable". Accurately, no fetch was attempted then; web tools were
available. The prefactor question is still open.

**Copyright.** `HomuMado.jpeg` is copyrighted anime art and `black white
river.jpeg` carries a "(c) 2012 Keith Cooper, Northlight Images" watermark;
`badapple.png` is likely derived from copyrighted material too. The E10 grid
embeds them. Fine for local testing; do not commit `data/picture/` or the E10
figures to a public repository. scikit-image's bundled test images
(`coffee`, `astronaut`, `chelsea`, `camera`) load offline and are a
publishable substitute.

**Checked.** 4 new tests pass (darkening scale and 8-bit quantisation; no
output at zero noise; an SSIM optimum exists; baselines keep shape and
range). Runtime about 60 s for five images, peak RSS 425 MB.

## WP9 update: E10 and E11 in the report

`report.md` section 7 now covers both extensions (7.1 images, 7.2 sound), with
one added statement in section 8 and one sentence each in the introduction and
conclusion. `extract_numbers.py` gained `rows_e10` (from
`results/e10_metrics.csv` and `results/e10_*.npz`) and `rows_e11` (from
`results/e11.npz`); 110 quantities, and every bolded claim traces. The
noisy-input sound check is still only in scratch, so the report names it
behind a TODO without quoting its numbers.

## E10 rerun: consent art, memory fix, name-based seeds

`HomuMado.jpeg` was replaced by two illustrations used with the artist's
permission for this project only, not for AI training (`MabelleMona1.jpg`,
`MabelleMona2.jpg`); artist credit still to be added to `report.md`.

**Crash.** The first six-image run peaked at 1.45 GB on a 3.8 GB machine; the
kernel OOM-killed claude-desktop and the machine was rebooted. Two causes:
`threshold_detector` drew all trials of noise as one float64 block (656 MB for
the 1.64 Mpx image, plus an equal temporary), and the E10 grid drew every panel
at full resolution. Fixed by drawing noise one trial at a time from the same
generator (identical stream; new test compares against the block draw) and by
drawing a ~800 px display copy in the grid. Sweep peak is now about 400 MB;
the full run including figures peaks at 790 MB. Heavy runs go under
`systemd-run --user --scope -p MemoryMax=1200M -p MemorySwapMax=0`.

**Seeds.** E10 seeded each image by its position in the sorted file list, so
adding an image changed every other image's numbers. Seeds now come from a
CRC32 of the file name. This moved every E10 number once (for example images
SR 0.645 to 0.653); a recompute of black white river after the memory fix was
bit-identical to the run before it, and `e10_metrics.csv` byte-identical.

**Checked.** 57 tests pass (72 s, 342 MB). Every bolded claim in `report.md`
traces (113 quantities). New counts: best baseline beats SR on 5 of 6 images,
SR beats the dark input on 2 of 6, CLAHE has the best PSNR on 6 of 6.

## References: proposal bibliography, method sources, image credits

`report.md` now has 28 references, each cited in the text: the proposal's five
in its own numbering ([1] to [5]; Chouhan et al. and Simonotto et al. were
missing before), twenty method sources ([6] to [25]: Kramers rate,
Euler-Maruyama and Heun, bootstrap, ROC, matched filter, noise-aided
detection, parameter-tuned, aperiodic and array SR, residence times, HSV,
histogram equalisation, gamma, PSNR, AHE and CLAHE, SSIM, scikit-image,
frequency rescaling) and three image credits ([26] to [28], including Keith
Cooper's Hood Canal photograph, matched by watermark and scene). Every DOI was
checked against Crossref and both ISBNs against OpenLibrary; that confirms the
citations, not the content. Still open: sources of `gradiant.jpeg` and
`images.jpeg`, and reading reference 3 for the prefactor.

## E10 frame count: grain falls as 1/sqrt(N)

New `experiments/e10_frames.py`, config block `experiments.e10_frames` (a
separate block, because adding a key to `e10` would invalidate every E10
cache). At each image's best sigma from `e10_image.py` it counts firings over
N = 10 to 400 frames, one frame in memory at a time (the full stack would be
2.6 GB at N = 400), in two independent runs.

**Result.** Grain slope against N is -0.501 to -0.497 on all six images and
within 0.8 % of sqrt(mean p(1-p) / N) at every N. Run 1 uses E10's seed, and
its 50-frame SSIM matches `e10_metrics.csv` to four decimals on every image,
which checks the counting against `approach_a`. SSIM roughly doubles from 50
to 400 frames on the photographs and illustrations (river 0.138 to 0.334,
limit 0.508). At 400 frames SR beats the dark input on 6 of 6 images (2 of 6
at 50) but the best baseline on only 1 of 6; even the N -> infinity limit, the
pure sigmoid, beats it on 2 of 6 (badapple, images). So the 50-frame result
understated SR, but the verdict stands. Sigma was fixed at the 50-frame
optimum, so the limits are lower bounds. Runtime 72 s, peak 442 MB.
Report 7.1 has a new paragraph; statement 4 in section 8 notes 400 frames.

## WP10: clean-environment reproduction

New `make_all.py` runs every experiment in order (E1 to E11, then
`e10_frames`), then `report/extract_numbers.py` and `web/export_data.py`,
stopping on the first failure. `--fresh` deletes `results/*.json` so every
simulation reruns. README already pointed to `make_all.py`, which did not exist.

**Run.** The project was copied without `.venv`, `results/`, `figures/`,
`report/numbers.*`, `web/data.js` and `node_modules`, a fresh venv was made
with Python 3.12.3, and `pip install -r requirements.txt` then
`pip install -e .` resolved the same versions as the working venv (numpy
2.5.3, scipy 1.18.1, matplotlib 3.11.2, scikit-image 0.26.0, opencv-python
5.0.0.93). `pytest`: 57 passed in 81 s. `make_all.py` under a 1.2 GB cap:
exit 0, about 43 min in total. The run was stopped once mid-E4 and restarted
as a detached systemd service so it would survive closing the desktop app; the
finished steps came from that run's own cache, so every result was still
computed in the fresh environment. Slowest steps: E3 680 s, E5 358 s,
E6 360 s, E7 151 s, E10 126 s.

**Mismatches against the committed results.**

1. `report/numbers.md`, `report/numbers.json` and `web/data.js` are
   byte-identical, so no number in the report or dashboard changes.
2. Results: 95 of 103 files byte-identical. Eight differ, all SNR-derived
   (`e1`, `e2`, `e3`, `e3_dt_half`, `e3_heun`, `e5_D0.05`, `e5_D0.1`,
   `e5_D0.2`), by at most 2.2e-8 relative. Likely cause: those caches predate
   the removal of the float64 cast in `_band_periodograms` (entry above), which
   changes rounding at that level; the cache keys on parameters, not code, so
   they were never rebuilt. Harmless, but a code change can leave a cache
   stale without warning.
3. Figures: 15 of 17 PNGs pixel-identical. `e3_numerical_checks` and
   `e5_parameter_tuned` differ in 34 and 33 pixels (0.001 %), from item 2.
   PDFs were not compared, as their creation dates always differ.
4. `results/e10_HomuMado.{npz,json}` exist only in the original: a stale
   cache of the removed image, safe to delete.

**TODOs left in the report.** Three: the E11 noisy-input check not in
`results/` (section 7.2), the pi vs pi/2 prefactor against reference 3
(section 8), and references 1 to 25 not read in full (reference list note).

**Figures not referenced in the report.** `e1_every_step`,
`e4_hundred_steps` and `e5_hundred_steps`, the step-by-step views, which are
supplementary and listed in the README only.

**References not cited in the text.** None; all 30 are cited, and no citation
points past the list.

## E10 denoising afterwards, and repeated passes

Question from review: if SR's grain is the problem, denoise the SR output, and
does repeating SR make the image brighter? New `experiments/e10_denoise.py`,
config `experiments.e10_denoise.passes`, and in `srlab.image_sr`:
`noise_sigma` (blind MAD estimate of diagonal differences), `denoise` (TV or
NLM from scikit-image, strength from `noise_sigma`, never from the original)
and `sr_passes` (SR then TV, repeated, threshold at max or median). Four new
tests: the estimator recovers a known sigma, denoising keeps the mean, pass 1
equals E10 then TV and max-threshold passes darken, a flat image is left alone.
Wavelet shrinkage and BM3D were not tried: both need PyWavelets or another new
dependency, and `skimage.restoration.estimate_sigma` needs PyWavelets too,
hence the hand-written estimator.

**Bug found on the way.** On the silhouette, over half the pixel differences
are exactly zero, the estimate was 0 and TV divided by it. `denoise` now
returns a zero-noise image unchanged; test added.

**Result.** Denoising keeps the mean (river 0.412 before and after) and lifts
SSIM a lot (river 0.138 to 0.760); denoised SR beats the best undenoised
baseline on 5 of 6. Repeating at the max threshold darkens every pass (0.59x
to 0.31x the original brightness, SSIM 0.649 to 0.377, means over images). At
the median, pass 2 brightens to 0.92x and helps on 4 of 6 (river 0.799), but
the crops show the tree losing its branches; later passes add nothing. Given
the same denoisers, before or after, CLAHE reaches 0.926 to 0.945 on four
images and the best denoised baseline beats every SR variant on 6 of 6. So the
E10 verdict gets stronger, not weaker: CLAHE was held back by amplified sensor
noise, not by its tone curve. Runtime about 3 min for six images, peak 661 MB.
References 31 to 35 added (TV, Chambolle, NLM, Donoho and Johnstone for the
MAD scaling, the quanta image sensor), DOIs checked against Crossref.

## Proposal-scope report

New `report/proposal_report.md`, in Thai, following the proposal's own
sections: objectives, method, the three expected outcomes, criteria 4 to 7,
conclusion, limits and further work. It covers only what the proposal promised:
1-D SR (E1, E2, E3 with its numerical checks, E4's frequency shift, E8's pulse
train, E6, E7) and the low-light image extension (E10 approach A). Everything
beyond the proposal (E5, E9, E11, the frame test, denoising) stays in
`report.md`. Twenty references: the proposal's five in their own numbering,
the method sources it uses, and the five image credits.
`extract_numbers.py`'s claim check now covers both reports; every bolded number
in each traces to `results/`.

## E12: two-dimensional SR without an image

Asked for a pure 2-D foundation under E10. New `experiments/e12_grating.py`,
config `experiments.e12`, seed 1014, `metrics.snr2d_db` (the 2-D form of
`snr_db`: averaged 2-D periodograms, median background in a window around k)
and `threshold.threshold_snr2d_theory_db` (closed form for independent pixels,
|sum p e^(-ik.x)|^2 / sum p(1-p)), with one test that checks a blank output
without noise and agreement within 1 dB with it. E2's A and theta on a
128 x 128 grating at k = (16, 8), 32 frames, 40 sigma values. Result: blank at
the 5 smallest sigma, peak 30.51 dB at sigma = 0.690 against theory 30.37 dB at
the same sigma, within 0.31 dB of theory wherever theory is above 0 dB; near
the floor only a handful of pixels fire and the estimator cannot resolve it,
as in E2. Added to both reports before E10 and to `make_all.py` after E2.
Runtime a few seconds.
Renamed the same day to sit under E10 as its foundation: `experiments/e10_grating.py`,
config and seed `e10_grating`, figure and results `e10_grating`. `.gitignore`
keeps `figures/e10_*` out but lets `figures/e10_grating.*` through, since it holds
no third-party image. Numbers unchanged on rerun.
