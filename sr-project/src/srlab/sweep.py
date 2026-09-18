"""Sweep runner with result caching (PLAN.md 3).

Ensembles are built in chunks so peak memory stays bounded: the raw traces of
one chunk are reduced to a band periodogram and thrown away before the next.
The 40 x 50 E3 ensemble would be 1.0 GB held whole.
"""

import json
import pathlib

import numpy as np

from .metrics import _SNR_FLOOR_DB, bootstrap_snr

__all__ = ["ROOT", "load_config", "d_grid", "chunked_snr", "cached"]

# Works for an editable install, where __file__ stays in the repo.
ROOT = pathlib.Path(__file__).resolve().parents[2]

# Per-chunk budget for the raw ensemble. About 300 MB leaves room on a 1 GB box.
_MAX_CHUNK_BYTES = 3.0e8


def load_config(path=None):
    """Load configs/default.yaml, or another config given by path."""
    import yaml

    with open(path or ROOT / "configs" / "default.yaml") as f:
        return yaml.safe_load(f)


def d_grid(cfg):
    """Log-spaced noise intensities from the config."""
    n = cfg["noise"]
    return np.logspace(np.log10(n["D_min"]), np.log10(n["D_max"]), n["n_D"])


def chunked_snr(
    make_ensemble,
    n_values,
    bytes_per_value,
    fs,
    f0,
    n_side=25,
    guard=1,
    n_boot=200,
    seed=0,
    max_bytes=_MAX_CHUNK_BYTES,
):
    """SNR with bootstrap band over a swept parameter, a chunk at a time.

    make_ensemble(lo, hi) returns the raw ensemble for sweep values [lo, hi)
    with shape (hi - lo, n_trials, n_samples). bytes_per_value is that array's
    size per swept value, used to pick the chunk width.

    A sweep value whose record is constant (the threshold detector never fires
    at low sigma) gets the SNR floor instead of raising, and says so.

    Returns (snr, lo, hi) as 1-D arrays of length n_values, in dB.
    """
    per_chunk = max(1, int(max_bytes // bytes_per_value))
    snr = np.empty(n_values)
    band_lo = np.empty(n_values)
    band_hi = np.empty(n_values)

    for start in range(0, n_values, per_chunk):
        stop = min(start + per_chunk, n_values)
        x = make_ensemble(start, stop)
        dead = x.std(axis=(-2, -1)) == 0
        for arr in (snr, band_lo, band_hi):
            arr[start:stop][dead] = _SNR_FLOOR_DB
        if dead.any():
            print(f"  [floor] {dead.sum()} sweep value(s) gave a constant record")
        if not dead.all():
            s, l, h = bootstrap_snr(
                x[~dead], fs, f0, n_side=n_side, guard=guard, n_boot=n_boot,
                seed=seed + start,
            )
            snr[start:stop][~dead] = s
            band_lo[start:stop][~dead] = l
            band_hi[start:stop][~dead] = h
        del x

    return snr, band_lo, band_hi


def cached(name, params, compute):
    """Return results/<name>.npz if it was built from these params, else build it.

    A figure can then be redrawn from results/ without rerunning the solver
    (PLAN.md 3). Delete the .json to force a rebuild.
    """
    results = ROOT / "results"
    results.mkdir(exist_ok=True)
    npz, meta = results / f"{name}.npz", results / f"{name}.json"

    if npz.exists() and meta.exists():
        with open(meta) as f:
            if json.load(f) == json.loads(json.dumps(params, default=_jsonable)):
                print(f"  [cached] {npz.relative_to(ROOT)}")
                return dict(np.load(npz))

    arrays = compute()
    np.savez_compressed(npz, **arrays)
    with open(meta, "w") as f:
        json.dump(params, f, indent=2, default=_jsonable, sort_keys=True)
    print(f"  [saved] {npz.relative_to(ROOT)}")
    return arrays


def _jsonable(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer, np.floating)):
        return o.item()
    raise TypeError(f"cannot serialise {type(o)}")
