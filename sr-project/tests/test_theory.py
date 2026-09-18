"""WP1: closed-form quantities (PLAN.md section 2.1)."""
import numpy as np
import pytest

from srlab import theory


def test_values_at_baseline():
    assert theory.well_position(1, 1) == pytest.approx(1.0)
    assert theory.barrier_height(1, 1) == pytest.approx(0.25)
    assert theory.critical_amplitude(1, 1) == pytest.approx(0.385, abs=5e-4)


def test_kramers_prefactor():
    """r_K = a / (sqrt(2) pi) * exp(-dU / D); prefactor 0.225 at a = b = 1."""
    assert theory.kramers_rate(1, 1, 1e9) == pytest.approx(1 / (np.sqrt(2) * np.pi))
    assert theory.kramers_rate(1, 1, 0.25) == pytest.approx(0.2251 * np.exp(-1), rel=1e-3)


def test_two_state_snr_peaks_at_half_the_barrier():
    """The peak sits at D = dU / 2 whatever the prefactor is."""
    for a in (0.5, 1.0, 2.0):
        barrier = theory.barrier_height(a, a)
        D = np.logspace(-3, 1, 20001)
        peak = D[theory.two_state_snr(a, a, 0.3, D).argmax()]
        assert peak == pytest.approx(barrier / 2, rel=1e-3)
