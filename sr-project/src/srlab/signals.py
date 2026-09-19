"""Signal and noise generators (PLAN.md 2).

Callers supply a numpy.random.Generator built from a seed in the config.
"""

import numpy as np

__all__ = ["snap_f0", "n_samples_for", "time_axis", "sine", "pulse_train",
           "white_noise"]


def snap_f0(f0, fs, n_periods):
    """Nearest frequency to f0 giving a whole number of samples in n_periods.

    An arbitrary f0, from a log-spaced grid say, leaves n_periods * fs / f0
    non-integer, so f0 misses its FFT bin and the record is not a whole number
    of solver steps. Snapping keeps the leakage-free property. The shift is
    tiny: on a 40-point grid it is well under one grid step.
    """
    n = round(n_periods * fs / f0)
    if n < 1:
        raise ValueError(f"f0={f0} is too high for {n_periods} periods at fs={fs}")
    return n_periods * fs / n


def n_samples_for(f0, fs, n_periods):
    """Samples in an exact integer number of periods, so f0 lands on a bin.

    Raises if f0, fs and n_periods do not give a whole number of samples.
    """
    n = n_periods * fs / f0
    if abs(n - round(n)) > 1e-9:
        raise ValueError(
            f"f0={f0}, fs={fs}, n_periods={n_periods} give {n} samples, not an "
            "integer; f0 would not land on an FFT bin"
        )
    return int(round(n))


def time_axis(f0, fs, n_periods):
    """Sample times [0, n_periods/f0), spacing 1/fs."""
    return np.arange(n_samples_for(f0, fs, n_periods)) / fs


def sine(A, f0, t):
    """Weak periodic drive A cos(2 pi f0 t)."""
    return A * np.cos(2 * np.pi * f0 * np.asarray(t, dtype=float))


def pulse_train(t, rate, width, amplitude, rng):
    """Poisson arrivals of rectangular pulses, for E8.

    Keep amplitude below critical_amplitude to stay sub-threshold.
    """
    t = np.asarray(t, dtype=float)
    starts = np.sort(rng.uniform(t[0], t[-1], size=rng.poisson(rate * (t[-1] - t[0]))))
    y = np.zeros_like(t)
    for s in starts:
        y[(t >= s) & (t < s + width)] = amplitude
    return y


def white_noise(shape, D, dt, rng):
    """Euler-Maruyama increment sqrt(2 D dt) randn. D broadcasts against shape.

    Scaling is sqrt(dt), not dt (PLAN.md 9).
    """
    return np.sqrt(2 * np.asarray(D, dtype=float) * dt) * rng.standard_normal(shape)
