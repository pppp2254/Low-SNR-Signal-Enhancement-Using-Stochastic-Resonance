"""Threshold detector, y = 1 if s(t) + noise > theta else 0 (PLAN.md 2.5).

Simplest system that shows SR: with A < theta the clean signal never crosses,
so noise is what carries it over. Reused by the image extension, approach A.
"""

import numpy as np
from scipy.special import erfc

__all__ = ["threshold_detector", "threshold_snr2d_theory_db"]


def threshold_detector(signal, sigma, theta, n_trials, seed):
    """Push a 1-D signal through a noisy hard threshold.

    sigma is a standard deviation, not the intensity D the solver takes. It may
    be an array, which becomes axis 0. Returns float32 zeros and ones of shape
    (n_sigma, n_trials, n_samples).
    """
    rng = np.random.default_rng(seed)
    signal = np.asarray(signal, dtype=float)
    if signal.ndim != 1:
        raise ValueError(f"signal must be 1-D, got shape {signal.shape}")
    sigma = np.atleast_1d(np.asarray(sigma, dtype=float))

    y = np.empty((sigma.size, n_trials, signal.size), dtype=np.float32)
    # One trial at a time: same random stream as a single (n_trials, n) draw,
    # without an n_trials x n float64 temporary (656 MB for a 1.6 Mpx image).
    for i, s in enumerate(sigma):
        for t in range(n_trials):
            y[i, t] = signal + s * rng.standard_normal(signal.size) > theta
    return y


def threshold_snr2d_theory_db(signal2d, theta, sigma, k):
    """Expected metrics.snr2d_db of threshold_detector on a 2-D pattern (E10 grating).

    Pixels fire independently with p = P(s + sigma xi > theta). The averaged
    periodogram is |sum p e^(-i k.x)|^2 + sum p(1 - p) at k and sum p(1 - p)
    elsewhere, so SNR = |sum p e^(-i k.x)|^2 / sum p(1 - p), floored at -30 dB.
    sigma may be an array.
    """
    s = np.asarray(signal2d, dtype=float)
    ny, nx = s.shape
    y, x = np.mgrid[0:ny, 0:nx]
    phase = np.exp(-2j * np.pi * (k[0] * x / nx + k[1] * y / ny))
    out = []
    for sg in np.atleast_1d(sigma):
        p = 0.5 * erfc((theta - s) / (sg * np.sqrt(2)))
        var = np.sum(p * (1 - p))
        out.append(-30.0 if var <= 0 else
                   max(10 * np.log10(np.abs(np.sum(p * phase)) ** 2 / var), -30.0))
    return np.array(out)
