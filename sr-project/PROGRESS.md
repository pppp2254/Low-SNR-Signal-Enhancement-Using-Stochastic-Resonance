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
