"""WP3: chunking and caching (PLAN.md 3)."""
import numpy as np
import pytest

from srlab import metrics
from srlab.sweep import cached, chunked_snr, d_grid, load_config

FS, F0, N = 10.0, 0.01, 8000  # 8 periods, f0 on bin 8


def ensemble(n_values, n_trials=6, amps=None, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(N) / FS
    amps = np.zeros(n_values) if amps is None else amps
    x = np.empty((n_values, n_trials, N))
    for i, A in enumerate(amps):
        x[i] = A * np.cos(2 * np.pi * F0 * t) + rng.standard_normal((n_trials, N))
    return x


def test_config_grid_matches_plan():
    cfg = load_config()
    D = d_grid(cfg)
    assert D.size == 40
    assert D[0] == pytest.approx(0.01) and D[-1] == pytest.approx(2.0)


def test_chunking_does_not_change_the_answer():
    """Same ensemble, one chunk or several: identical SNR."""
    amps = np.array([0.05, 0.1, 0.2, 0.4])
    x = ensemble(4, amps=amps)
    per_value = x[0].nbytes

    def make(lo, hi):
        return x[lo:hi]

    whole = chunked_snr(make, 4, per_value, FS, F0, n_side=5, n_boot=20,
                        max_bytes=per_value * 4)
    split = chunked_snr(make, 4, per_value, FS, F0, n_side=5, n_boot=20,
                        max_bytes=per_value)
    assert np.allclose(whole[0], split[0])
    assert np.all(np.diff(whole[0]) > 0)


def test_constant_sweep_value_is_floored_not_raised():
    """The threshold detector emits nothing at low sigma; the sweep must survive."""
    x = ensemble(3, amps=[0.0, 0.2, 0.4])
    x[0] = 0.0  # detector never fired

    snr, lo, hi = chunked_snr(lambda a, b: x[a:b], 3, x[0].nbytes, FS, F0,
                              n_side=5, n_boot=20)
    assert snr[0] == metrics._SNR_FLOOR_DB
    assert np.isfinite(snr).all() and snr[2] > snr[1]


def test_degenerate_bootstrap_draw_is_floored_but_point_estimate_raises():
    """Resampling only silent trials is expected near a detector's floor."""
    x = np.zeros((1, 8, N))
    x[0, 0] = np.random.default_rng(0).standard_normal(N)  # one trial fires

    snr, lo, hi = metrics.bootstrap_snr(x, FS, F0, n_side=5, n_boot=50)
    assert np.isfinite(snr).all() and lo <= snr <= hi

    with pytest.raises(ValueError, match="background power is zero"):
        metrics.snr_db(np.zeros((1, 8, N)), FS, F0, n_side=5)


def test_cache_round_trips_and_rebuilds_on_param_change(tmp_path, monkeypatch):
    import srlab.sweep as sweep

    monkeypatch.setattr(sweep, "ROOT", tmp_path)
    calls = []

    def compute():
        calls.append(1)
        return {"y": np.arange(5.0)}

    first = sweep.cached("demo", {"n": 1}, compute)
    second = sweep.cached("demo", {"n": 1}, compute)
    assert np.array_equal(first["y"], second["y"])
    assert len(calls) == 1, "second call should have hit the cache"

    sweep.cached("demo", {"n": 2}, compute)
    assert len(calls) == 2, "changed params must rebuild"
