"""Vectorised Langevin solver for the double well (PLAN.md 2.1, 2.3).

dx/dt = a x - b x^3 + A cos(2 pi f0 t) + sqrt(2 D) xi(t)

One array of shape (n_D, n_trials) holds the whole grid and is stepped together.
"""

import numpy as np

from .theory import well_position

__all__ = ["simulate", "residence_times"]


def simulate(
    a,
    b,
    A,
    f0,
    D,
    dt,
    n_periods,
    discard_periods,
    decimate,
    n_trials,
    seed,
    method="em",
    input_override=None,
):
    """Integrate the Langevin equation over a grid of noise intensities.

    D is a scalar or 1-D array and becomes axis 0 of the output. dt must be
    0.01 or smaller or the cubic term blows up at large D. discard_periods is
    a whole number of periods, so the drive is back at phase zero when
    recording starts and the output compares directly with cos(2 pi f0 t).
    decimate sets fs = 1/(dt*decimate). method is "em" or "heun"; the noise is
    additive, so Heun reuses the same increment in both stages.

    input_override replaces the drive with a pre-built one, shaped
    (..., n_discard_steps + n_record_steps + 1) and broadcastable against
    (n_D, n_trials). Used by E5. It holds one value per solver step, not per
    recorded sample, so it is large.

    Returns (t, x) with t starting at 0 and x float32 of shape
    (n_D, n_trials, n_samples).

    Memory: n_D * n_trials * n_samples * 4 bytes, 410 MB for the 40 x 20
    baseline and 1.0 GB at 50 trials. Call in chunks of D on a small machine.
    """
    if method not in ("em", "heun"):
        raise ValueError(f"method must be 'em' or 'heun', got {method!r}")

    rng = np.random.default_rng(seed)
    D = np.atleast_1d(np.asarray(D, dtype=float))
    if D.ndim != 1:
        raise ValueError(f"D must be scalar or 1-D, got shape {D.shape}")
    n_D = D.size

    n_discard = int(round(discard_periods / f0 / dt))
    n_record = int(round(n_periods / f0 / dt))
    if n_record % decimate:
        raise ValueError(f"{n_record} recorded steps is not a multiple of {decimate}")
    n_samples = n_record // decimate
    n_total = n_discard + n_record

    if input_override is None:
        drive = A * np.cos(2 * np.pi * f0 * (np.arange(n_total + 1) * dt))
    else:
        drive = np.asarray(input_override, dtype=float)
        if drive.shape[-1] < n_total + 1:
            raise ValueError(
                f"input_override has {drive.shape[-1]} steps, need {n_total + 1}"
            )

    # Random +-x_m per trial (PLAN.md 2.3).
    x_m = well_position(a, b)
    x = x_m * (2.0 * rng.integers(0, 2, size=(n_D, n_trials)) - 1.0)
    noise_scale = np.sqrt(2 * D * dt).reshape(n_D, 1)
    shape = (n_D, n_trials)
    out = np.empty((n_D, n_trials, n_samples), dtype=np.float32)

    def step(k):
        nonlocal x
        f = a * x - b * x * x * x + drive[..., k]
        dW = noise_scale * rng.standard_normal(shape)
        if method == "em":
            x = x + f * dt + dW
        else:
            xp = x + f * dt + dW
            fp = a * xp - b * xp * xp * xp + drive[..., k + 1]
            x = x + 0.5 * (f + fp) * dt + dW

    for k in range(n_discard):
        step(k)

    for j in range(n_samples):
        out[:, :, j] = x
        for i in range(decimate):
            step(n_discard + j * decimate + i)

    return np.arange(n_samples) * (dt * decimate), out


def residence_times(x, dt_sample, x_m=1.0):
    """Durations spent in one well between committed transitions.

    A state counts only once |x| exceeds x_m / 2. Bare sign changes overcount
    badly: near the barrier top a trajectory recrosses many times per real
    transition. Returns an empty array when the trace never crosses.

    At the optimal noise the histogram of these times peaks at odd multiples of
    half the drive period, which is what synchronisation means: the particle
    waits for the potential to tilt its way, hops, then waits a half period.
    """
    x = np.asarray(x)
    idx = np.flatnonzero(np.abs(x) > x_m / 2)
    if idx.size < 2:
        return np.empty(0)
    jumps = np.flatnonzero(np.diff(np.sign(x[idx])) != 0)
    if jumps.size < 2:
        return np.empty(0)
    return np.diff(idx[jumps + 1]) * dt_sample
