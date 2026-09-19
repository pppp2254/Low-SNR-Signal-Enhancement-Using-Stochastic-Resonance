# Stochastic resonance (sr-project)

Does adding noise to a weak sub-threshold signal make it *easier* to detect?
This repo simulates an overdamped particle in a double well driven by a weak
sine plus white noise, and measures output SNR as noise intensity D is swept.
The headline result (E3) is an SNR curve that rises, peaks near D = 0.125, then
falls.

The full specification is in [PLAN.md](PLAN.md): the model in section 2, the
repository rules in section 3, the ten experiments in section 4, and the work
packages in section 5.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
pytest
```

## Reproduce a figure

Each experiment is a thin script that imports `srlab`; none of them define
their own solver or SNR estimator.

```bash
python experiments/e3_snr_vs_noise.py
```

Every script writes raw arrays to `results/<name>.npz`, the exact parameters to
`results/<name>.json`, and the figure to `figures/<name>.png` and `.pdf`.
A figure can be rebuilt from `results/` without rerunning the simulation.

To regenerate everything from scratch:

```bash
python make_all.py
```

## Where the figures come from

Every figure is drawn by matplotlib inside an experiment script. Nothing is
drawn by hand, and the dashboard does not redraw anything: it displays the same
PNG files.

The chain is one way:

```
configs/default.yaml      parameters, the single source for every number
   |
src/srlab/                solver, metrics, theory  (no plotting logic)
   |
experiments/eN_*.py       run the sweep, then draw the figure
   |-> results/<name>.npz raw arrays,  results/<name>.json the exact parameters
   |-> figures/<name>.png 300 dpi,  figures/<name>.pdf vector
   |
web/export_data.py        reads results/, writes web/data.js (numbers only, 3 KB)
   |
web/index.html            shows figures/*.png plus the tables
web/server.ts             Elysia on Bun serves the page and the figures, port 3100
```

| Figure | Drawn by |
| --- | --- |
| `e1_regimes` | `experiments/e1_regimes.py` |
| `e2_threshold` | `experiments/e2_threshold.py` |
| `e3_snr_vs_noise` | `experiments/e3_snr_vs_noise.py` |
| `e3_numerical_checks` | `experiments/e3_snr_vs_noise.py` |
| `e4_robustness` | `experiments/e4_robustness.py` |
| `e5_parameter_tuned` | `experiments/e5_parameter_tuned.py` |
| `e6_detection` | `experiments/e6_detection.py` |
| `e7_baselines` | `experiments/e7_baselines.py` |

Style lives in `src/srlab/plotting.py`, so the figures match each other:

- `use_style()` sets matplotlib rcParams once: 300 dpi for saved files, tight
  bounding box, 9 pt base font, grid at 0.25 alpha, no top or right spine.
- `COLORS` fixes the series palette: `sim` #1f4e79, `theory` #c0392b,
  `clean` #7f8c8d, `optimal` #e67e22, `alt1` #2e8b57, `alt2` #8e44ad. The E4
  heatmap uses matplotlib's `viridis` through `pcolormesh`.
- `save(fig, name)` writes both `figures/<name>.png` and `.pdf`.
- The backend is `Agg`, so scripts run headless and never open a window.

Two consequences worth knowing:

1. A figure can be redrawn from `results/` without rerunning the solver. The
   scripts cache on the exact parameter set, so rerunning one is close to
   instant unless a parameter changed. Delete `results/<name>.json` to force a
   rebuild.
2. The dashboard cannot drift from the figures, because it serves the same
   files. An earlier version redrew the charts in JavaScript and did drift: its
   colours resolved to the dark-theme palette instead of the matplotlib one.
   That code was deleted rather than patched.

## Layout

| Path | Contents |
| --- | --- |
| `src/srlab/` | all reusable code: solver, metrics, theory, plotting |
| `experiments/` | one script per experiment, `e1_regimes.py` to `e10_image.py` |
| `configs/default.yaml` | baseline parameters (PLAN.md section 2.2) |
| `tests/` | pytest unit tests |
| `results/` | `.npz` data and `.json` parameters |
| `figures/` | 300 dpi `.png` and `.pdf` |
| `notebooks/demo.ipynb` | Colab walkthrough that calls `srlab` |
| `web/` | self-contained interactive demo |
| `report/` | report and slides |

## Rules

- All randomness goes through `numpy.random.default_rng(seed)`; seeds live in
  the config.
- Missing function? Add it to `srlab` with a test, not to a script.
- Never tune parameters silently to make a peak appear. Any change from
  `configs/default.yaml` is logged in `PROGRESS.md` with the reason.

## Status

WP0 to WP5 done. `signals`, `bistable`, `threshold`, `theory`, `metrics`,
`sweep` and `plotting` are implemented and tested (44 tests). E1, E2 and E3
have run; E3 puts the SNR peak at D = 0.1153 against a theoretical 0.125, and
the halved-step and Heun checks agree to within 0.7 dB. E4 confirms that the
optimal D rises with f0. E5 does not reproduce its tuning optimum, and found
that PLAN.md 2.5's a = 8 D_in is a D-sweep condition misapplied to an a sweep;
the correct stationary point is a = 4 D_in. E6 and E7 answer the detection
criterion: a linear band-pass beats every nonlinear method by 13.6 dB, while
adding noise is still worth +8.41 dB to the bistable system, and the threshold
detector goes from chance to perfect to chance as noise alone is increased.

Not started: E8 to E10, the report and the notebook (its experiment cells are
still stubs). `image_sr` is a docstring stub. See `PROGRESS.md` for open
issues, including an unverified theory prefactor.

## Dashboard

```bash
cd web && bun install && bun run dev
```

Serves a raw results readout on port 3100: the measured values, the six
figures, and the full E3 table.

The page has a rerun control per experiment. Clicking one runs that experiment
script, then re-runs `web/export_data.py`, then reloads the page so the new
figures and numbers appear. Because each script caches on its exact parameter
set, a rerun with nothing changed only redraws and takes seconds; change a
value in `configs/default.yaml` and the same button recomputes it.

The controls need the Elysia server. The published artifact is static, so
`/run/status` is unreachable there and the controls stay hidden.

`server.ts` runs only the scripts in a fixed whitelist and never builds a
command from request input. It is local development tooling; do not expose the
port to a network.

Rebuild the dashboard numbers by hand if you ran a script from the shell:

```bash
.venv/bin/python web/export_data.py
```
