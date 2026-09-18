"""WP1: generators (PLAN.md section 2)."""
import numpy as np
import pytest

from srlab import signals


def test_f0_on_a_bin():
    assert signals.n_samples_for(f0=0.01, fs=10, n_periods=128) == 128_000
    t = signals.time_axis(f0=0.01, fs=10, n_periods=128)
    x = signals.sine(1.0, 0.01, t)
    spectrum = np.abs(np.fft.rfft(x))
    assert spectrum.argmax() == 128  # bin index equals the number of periods
    with pytest.raises(ValueError):
        signals.n_samples_for(f0=0.0033, fs=10, n_periods=128)


def test_noise_scales_with_sqrt_dt():
    """sqrt(2 D dt), not 2 D dt (PLAN.md section 9)."""
    rng = np.random.default_rng(0)
    D, dt = 0.3, 0.01
    w = signals.white_noise((200_000,), D, dt, rng)
    assert w.std() == pytest.approx(np.sqrt(2 * D * dt), rel=0.02)


def test_noise_broadcasts_over_D():
    rng = np.random.default_rng(0)
    D = np.array([[0.01], [1.0]])
    w = signals.white_noise((2, 100_000), D, 0.01, rng)
    assert w.std(axis=1)[1] / w.std(axis=1)[0] == pytest.approx(10.0, rel=0.02)


def test_pulse_train_is_sub_threshold_and_rectangular():
    rng = np.random.default_rng(1)
    t = np.arange(0, 2000, 0.1)
    y = signals.pulse_train(t, rate=0.01, width=20.0, amplitude=0.3, rng=rng)
    assert set(np.unique(y)) <= {0.0, 0.3}
    assert 0.0 < y.mean() < 0.3
