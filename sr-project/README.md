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

WP0 to WP3 done. `signals`, `bistable`, `threshold`, `theory`, `metrics`,
`sweep` and `plotting` are implemented and tested (44 tests). E1, E2 and E3
have run; E3 puts the SNR peak at D = 0.1153 against a theoretical 0.125, and
the halved-step and Heun checks agree to within 0.7 dB.

Not started: E4 to E10, the report and the notebook (its experiment cells are
still stubs). `image_sr` is a docstring stub. See `PROGRESS.md` for open
issues, including an unverified theory prefactor.

## Dashboard

```bash
cd web && bun install && bun run dev
```

Serves a raw results readout on port 3100. Rebuild its data after any
experiment reruns:

```bash
.venv/bin/python web/export_data.py
```
