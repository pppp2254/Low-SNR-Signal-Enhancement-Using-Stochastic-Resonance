"""WP1: the threshold detector (PLAN.md section 2.5)."""
import numpy as np
import pytest

from srlab import signals
from srlab.threshold import threshold_detector


def signal():
    t = signals.time_axis(f0=0.01, fs=10, n_periods=8)
    return signals.sine(0.7, 0.01, t)


def test_output_is_binary_with_the_right_shape():
    y = threshold_detector(signal(), sigma=[0.2, 1.0], theta=1.0, n_trials=3, seed=1)
    assert y.shape == (2, 3, 8000)
    assert y.dtype == np.float32
    assert set(np.unique(y)) <= {0.0, 1.0}


def test_sub_threshold_signal_never_fires_without_noise():
    """A = 0.7 < theta = 1: no noise, no output at all. That is the point."""
    y = threshold_detector(signal(), sigma=0.0, theta=1.0, n_trials=2, seed=1)
    assert y.sum() == 0


def test_firing_rate_grows_with_noise():
    sigmas = [0.05, 0.3, 1.0, 3.0]
    y = threshold_detector(signal(), sigma=sigmas, theta=1.0, n_trials=4, seed=2)
    rates = y.mean(axis=(1, 2))
    assert np.all(np.diff(rates) > 0)
    assert rates[0] < 1e-3 and rates[-1] > 0.2


def test_same_seed_reproduces():
    kw = dict(signal=signal(), sigma=0.5, theta=1.0, n_trials=3)
    assert np.array_equal(
        threshold_detector(**kw, seed=9), threshold_detector(**kw, seed=9)
    )
    assert not np.array_equal(
        threshold_detector(**kw, seed=9), threshold_detector(**kw, seed=10)
    )
