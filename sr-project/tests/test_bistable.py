"""WP1: the Langevin solver (PLAN.md sections 2.1, 2.3).

The checks are the ones PLAN.md's WP1 prompt names: a sub-threshold signal
cannot escape its well without noise, a super-threshold one can, the noise
strength is right (switching rate against Kramers), the seed reproduces, and
Heun agrees with Euler-Maruyama.
"""
import numpy as np
import pytest

from srlab.bistable import simulate
from srlab.theory import critical_amplitude, kramers_rate

BASE = dict(a=1.0, b=1.0, f0=0.01, dt=0.01, decimate=10)


def switch_rate(x, dt_sample, x_m=1.0):
    """Transitions per unit time, counting only committed well-to-well moves.

    A bare sign change overcounts: near the barrier top the trajectory
    recrosses many times per real transition. Requiring |x| > x_m / 2 before a
    state is registered removes those.
    """
    committed = np.sign(x[np.abs(x) > x_m / 2])
    n = np.count_nonzero(np.diff(committed) != 0)
    return n / (x.size * dt_sample)


def test_sub_threshold_signal_cannot_escape_without_noise():
    """A = 0.3 < A_c = 0.385: the noiseless particle stays where it started."""
    assert 0.3 < critical_amplitude(1, 1)
    _, x = simulate(**BASE, A=0.3, D=0.0, n_periods=8, discard_periods=0,
                    n_trials=6, seed=11)
    signs = np.sign(x[0])
    assert np.all(signs == signs[:, :1]), "trajectory left its starting well at D = 0"
    assert np.all(np.abs(x) > 0.5)


def test_super_threshold_signal_crosses_without_noise():
    """A = 0.5 > A_c: the well vanishes each half cycle, so the particle must cross."""
    assert 0.5 > critical_amplitude(1, 1)
    _, x = simulate(**BASE, A=0.5, D=0.0, n_periods=3, discard_periods=0,
                    n_trials=6, seed=11)
    for trial in x[0]:
        assert np.any(np.sign(trial) != np.sign(trial[0])), "no crossing at A = 0.5"


def test_switching_rate_matches_kramers():
    """With no drive, the measured rate is within 20% of r_K for D = 0.08 to 0.15."""
    Ds = [0.10, 0.15]
    t, x = simulate(**BASE, A=0.0, D=Ds, n_periods=48, discard_periods=2,
                    n_trials=6, seed=7)
    dt_sample = t[1] - t[0]
    for i, D in enumerate(Ds):
        measured = np.mean([switch_rate(trial, dt_sample) for trial in x[i]])
        assert measured == pytest.approx(kramers_rate(1, 1, D), rel=0.20), f"D={D}"


def test_same_seed_reproduces():
    kw = dict(**BASE, A=0.3, D=[0.05, 0.2], n_periods=4, discard_periods=1, n_trials=3)
    _, first = simulate(**kw, seed=42)
    _, again = simulate(**kw, seed=42)
    _, other = simulate(**kw, seed=43)
    assert np.array_equal(first, again)
    assert not np.array_equal(first, other)


def test_heun_agrees_with_euler_maruyama():
    """The numerical-accuracy check of PLAN.md section 2.3, at the statistics level."""
    kw = dict(**BASE, A=0.0, D=[0.12], n_periods=48, discard_periods=2, n_trials=6, seed=3)
    t, x_em = simulate(**kw, method="em")
    _, x_heun = simulate(**kw, method="heun")
    dt_sample = t[1] - t[0]
    r_em = np.mean([switch_rate(trial, dt_sample) for trial in x_em[0]])
    r_heun = np.mean([switch_rate(trial, dt_sample) for trial in x_heun[0]])
    assert r_em == pytest.approx(r_heun, rel=0.10)
    assert x_em.std() == pytest.approx(x_heun.std(), rel=0.05)


def test_output_shape_and_dtype():
    t, x = simulate(**BASE, A=0.3, D=[0.05, 0.1, 0.2], n_periods=2, discard_periods=0,
                    n_trials=4, seed=1)
    assert x.shape == (3, 4, 2 * 100 * 10)  # n_D, n_trials, n_periods * (1/f0) * fs
    assert x.dtype == np.float32
    assert t[0] == 0.0 and t[1] - t[0] == pytest.approx(0.1)


def test_input_override_drives_the_system():
    """E5 pushes one pre-built noisy input through systems with different a."""
    n_total = int(round(3 / BASE["f0"] / BASE["dt"])) + 1
    drive = np.full(n_total, 0.6)  # constant tilt, above A_c
    _, x = simulate(**BASE, A=0.0, D=0.0, n_periods=3, discard_periods=0,
                    n_trials=4, seed=5, input_override=drive)
    assert np.all(x[0, :, -1] > 0), "a constant super-critical tilt should pull every trial to +x"
    with pytest.raises(ValueError):
        simulate(**BASE, A=0.0, D=0.0, n_periods=3, discard_periods=0, n_trials=1,
                 seed=5, input_override=np.zeros(10))
