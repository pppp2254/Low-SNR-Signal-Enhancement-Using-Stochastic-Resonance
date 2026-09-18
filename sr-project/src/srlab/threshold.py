"""Threshold detector, y = 1 if s(t) + noise > theta else 0 (PLAN.md 2.5).

Simplest system that shows SR: with A < theta the clean signal never crosses,
so noise is what carries it over. Reused by the image extension, approach A.
"""

import numpy as np

__all__ = ["threshold_detector"]


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
    for i, s in enumerate(sigma):
        y[i] = signal + s * rng.standard_normal((n_trials, signal.size)) > theta
    return y
