"""Closed-form quantities for the double well (PLAN.md 2.1).

U(x) = -a x^2/2 + b x^4/4, driven by A cos(2 pi f0 t) plus noise of intensity D.
Dimensionless model units throughout, so rates are per unit model time.
"""

import numpy as np

__all__ = [
    "well_position",
    "barrier_height",
    "critical_amplitude",
    "kramers_rate",
    "two_state_snr",
    "two_state_snr_db",
]


def well_position(a, b):
    """Minima at x_m = +-sqrt(a/b). 1.0 at a = b = 1."""
    return np.sqrt(a / b)


def barrier_height(a, b):
    """dU = U(0) - U(x_m) = a^2/(4b). 0.25 at a = b = 1."""
    return a**2 / (4 * b)


def critical_amplitude(a, b):
    """A_c = sqrt(4a^3/(27b)). 0.385 at a = b = 1.

    Below A_c the tilted potential keeps two minima, so a noiseless particle
    never escapes. That is what "sub-threshold" means here.
    """
    return np.sqrt(4 * a**3 / (27 * b))


def kramers_rate(a, b, D):
    """Escape rate from one well, r_K = a/(sqrt(2) pi) exp(-dU/D).

    From r = sqrt(|U''(0)| U''(x_m))/(2 pi) exp(-dU/D) with U''(0) = -a and
    U''(x_m) = 2a. Prefactor 0.2251 at a = b = 1. Mean residence time is 1/r_K.
    Asymptotic for dU >> D, so only a rough guide at the optimum (dU/D = 2).
    """
    D = np.asarray(D, dtype=float)
    return (a / (np.sqrt(2) * np.pi)) * np.exp(-barrier_height(a, b) / D)


def two_state_snr(a, b, A, D):
    """Two-state SNR, SNR = pi (A x_m / D)^2 r_K(D) (PLAN.md 2.1).

    Units of inverse time, not dB. Converting to the per-bin ratio that
    metrics.snr_db reports needs the analysis bandwidth df = f0/n_periods;
    see theory_snr_db in this module's caller (experiments/e3_snr_vs_noise.py).

    Peak is at D = dU/2 whatever the prefactor, since
    d/dD log SNR = -2/D + dU/D^2.

    Frequency drops out: the response amplitude and the noise background share
    the same Lorentzian in 4 r_K^2 + W^2, so it cancels. That matters because
    the baseline is not adiabatic (W = 0.063 against 2 r_K = 0.061 at D = 0.125).

    Prefactor unverified. PLAN.md asks for a check against Gammaitoni et al.
    (1998); the paper was unreachable here, so no equation number is cited.
    An independent derivation gives pi/2, a factor of 2 below this, which is a
    one-sided vs two-sided spectrum convention: 3 dB vertical, peak unmoved.
    PLAN.md's pi is used so the project matches its own spec. See PROGRESS.md.
    """
    D = np.asarray(D, dtype=float)
    return np.pi * (A * well_position(a, b) / D) ** 2 * kramers_rate(a, b, D)


def two_state_snr_db(a, b, A, D, f0, n_periods):
    """two_state_snr converted to the per-bin ratio metrics.snr_db reports.

    snr_db gives xbar^2 T_rec / (2 S_N) and two_state_snr gives pi xbar^2 / S_N,
    so the per-bin value is R T_rec / (2 pi) with T_rec = n_periods / f0. Assumes
    the S_N in the two-state formula is a one-sided PSD in ordinary frequency.
    That assumption and the pi prefactor are both unverified; neither moves the
    peak. See PROGRESS.md.
    """
    R = two_state_snr(a, b, A, D)
    return 10 * np.log10(R * (n_periods / f0) / (2 * np.pi))
