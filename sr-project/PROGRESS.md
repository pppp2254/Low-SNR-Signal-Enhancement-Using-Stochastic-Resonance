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
