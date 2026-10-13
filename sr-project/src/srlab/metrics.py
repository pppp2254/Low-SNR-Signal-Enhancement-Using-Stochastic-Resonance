"""SNR, correlation and detection metrics (PLAN.md 2.4).

snr_db is a per-bin ratio: power in the f0 bin over background power per bin.
It depends on the analysis bandwidth df = fs/n_samples, so comparisons must
use the same record length.
"""

import numpy as np
from scipy.signal import correlate, correlation_lags
from scipy.stats import rankdata

__all__ = ["snr_db", "snr2d_db", "bootstrap_snr", "f0_power", "fidelity_db", "xcorr_peak",
           "roc_auc", "pd_at_pfa"]

# Keeps peak FFT memory near 30 MB whatever the ensemble size.
_FFT_CHUNK_ELEMENTS = 2_000_000

# With the background subtracted, S goes non-positive when there is no signal
# and the log runs to -inf. Floor it well below what the estimator can resolve:
# averaging M trials leaves the f0 bin with a spread of 1/sqrt(M), about 0.22 at
# 20 trials, so anything under roughly -7 dB already means "not there".
_SNR_FLOOR_DB = -30.0


def _bin_index(n_samples, fs, f0):
    """Index of the f0 bin. Raises unless f0 lands on one exactly (PLAN.md 9)."""
    k = f0 * n_samples / fs
    if abs(k - round(k)) > 1e-9:
        raise ValueError(
            f"f0={f0} at fs={fs} over {n_samples} samples falls on bin {k}, not an "
            "integer; use an exact integer number of signal periods"
        )
    return int(round(k))


def _band_periodograms(x, fs, f0, n_side):
    """Per-trial one-sided periodogram over bins f0-n_side .. f0+n_side.

    Shape x.shape[:-1] + (2*n_side+1,), with f0 at index n_side. Only the band
    is kept, so a 410 MB ensemble reduces to kilobytes without ever holding a
    full spectrum.
    """
    # Do not cast the whole ensemble to float64 here: that copy is twice the
    # size of the float32 input and defeats the chunking below. numpy's rfft
    # promotes internally, so each chunk is upcast on its own.
    x = np.asarray(x)
    if x.ndim < 2:
        raise ValueError(
            f"expected (..., n_trials, n_samples), got {x.shape}; wrap one trace as x[None]"
        )
    n = x.shape[-1]
    k0 = _bin_index(n, fs, f0)
    lo, hi = k0 - n_side, k0 + n_side + 1
    if lo < 1:
        raise ValueError(
            f"f0 is on bin {k0}; n_side={n_side} would reach the DC bin. "
            "Record more periods or lower n_side."
        )

    flat_x = x.reshape(-1, n)
    out = np.empty((flat_x.shape[0], 2 * n_side + 1))
    chunk = max(1, _FFT_CHUNK_ELEMENTS // n)
    scale = 2.0 / (fs * n)
    for i in range(0, flat_x.shape[0], chunk):
        spec = np.fft.rfft(flat_x[i : i + chunk], axis=-1)[:, lo:hi]
        out[i : i + chunk] = scale * (spec.real**2 + spec.imag**2)
    return out.reshape(x.shape[:-1] + (2 * n_side + 1,))


def _snr_from_band(band, n_side, guard, on_degenerate="raise"):
    """Average the periodograms, then take the ratio.

    Averaging first is what stabilises the background (PLAN.md 9). A single
    trial's periodogram is exponential, whose median sits at 0.693 of the mean
    and biases N low by 1.6 dB; by 20 trials that is under 0.1 dB.

    on_degenerate says what to do when the background is zero. "raise" is for
    the point estimate, where a constant record means something is wrong.
    "floor" is for bootstrap draws, where resampling only silent trials is an
    expected outcome near the bottom of a detector sweep, not a bug.
    """
    p = band.mean(axis=-2)
    offsets = np.arange(-n_side, n_side + 1)
    side = np.abs(offsets) > guard
    if not side.any():
        raise ValueError(f"guard={guard} leaves no background bins at n_side={n_side}")

    N = np.median(p[..., side], axis=-1)
    degenerate = N <= 0
    if np.any(degenerate):
        if on_degenerate == "raise":
            raise ValueError(
                "background power is zero, so SNR is undefined; the record is constant "
                "over the analysis band (an all-zero detector output does this)"
            )
        N = np.where(degenerate, 1.0, N)

    S = np.maximum(p[..., n_side] - N, N * 10 ** (_SNR_FLOOR_DB / 10))
    snr = 10 * np.log10(S / N)
    return np.where(degenerate, _SNR_FLOOR_DB, snr), S, N


def snr_db(x, fs, f0, n_side=25, guard=1):
    """Output SNR at f0 in dB, floored at -30 dB (PLAN.md 2.4).

    x is (..., n_trials, n_samples); periodograms are averaged over axis -2 and
    leading axes are preserved. guard excludes bins within that distance of f0,
    so the background is the median of bins f0 +-(guard+1) .. f0 +-n_side.
    Default 1 follows the WP2 signature; PLAN.md 2.4's prose implies 2. With an
    exact integer number of periods there is no leakage and the two agree to
    hundredths of a dB.

    Returns (snr, S, N): dB, signal power at f0 less background, background per bin.
    """
    return _snr_from_band(_band_periodograms(x, fs, f0, n_side), n_side, guard)


def snr2d_db(frames, k, n_side=6, guard=1):
    """Spatial SNR of a 2-D pattern at wavenumber k = (kx, ky), in dB (E10 grating).

    snr_db in two dimensions: frames is (n_frames, ny, nx), periodograms are
    averaged over frames, and the background is the median of a
    (2 n_side + 1)^2 window around k with a (2 guard + 1)^2 hole at k. Keep the
    window clear of DC and of 2k. Floored at -30 dB like snr_db.
    """
    p = (np.abs(np.fft.fft2(np.asarray(frames, dtype=float))) ** 2).mean(axis=0)
    kx, ky = k
    off = np.arange(-n_side, n_side + 1)
    win = p[np.ix_((ky + off) % p.shape[0], (kx + off) % p.shape[1])]
    hole = (np.abs(off)[:, None] <= guard) & (np.abs(off)[None, :] <= guard)
    noise = np.median(win[~hole])
    if noise <= 0:
        return _SNR_FLOOR_DB
    sig = max(p[ky % p.shape[0], kx % p.shape[1]] - noise, noise * 10 ** (_SNR_FLOOR_DB / 10))
    return float(10 * np.log10(sig / noise))


def bootstrap_snr(x, fs, f0, n_side=25, guard=1, n_boot=200, seed=0):
    """SNR with a bootstrap band over trials (PLAN.md 2.4).

    Periodograms are computed once and reused, so the resamples are near free.
    Returns (snr, lo, hi): the all-trials point estimate to plot as the marker,
    and the 16th and 84th bootstrap percentiles for the error bar.
    """
    band = _band_periodograms(x, fs, f0, n_side)
    snr, _, _ = _snr_from_band(band, n_side, guard)

    rng = np.random.default_rng(seed)
    n_trials = band.shape[-2]
    draws = np.empty((n_boot,) + np.shape(snr))
    for i in range(n_boot):
        idx = rng.integers(0, n_trials, size=n_trials)
        draws[i] = _snr_from_band(band[..., idx, :], n_side, guard, "floor")[0]
    lo, hi = np.percentile(draws, [16, 84], axis=0)
    return snr, lo, hi


def f0_power(x, fs, f0):
    """Periodogram power at f0 for each trial separately.

    The detection statistic of E6. Unlike snr_db this does not average over
    trials: each record is one observation, and the spread across records is
    what the ROC curve measures.
    """
    return _band_periodograms(x, fs, f0, 1)[..., 1]


def fidelity_db(reference, output):
    """How well output reproduces a broadband reference waveform, in dB.

    The SNR of the best linear fit a + b * output to the reference:
    10 log10(var(reference) / var(residual)). Equivalent to
    10 log10(1 / (1 - rho^2)) with rho the correlation, and to 10 log10(1 + SNR)
    when output is the reference plus independent noise. A constant output
    carries nothing and scores 0 dB. For sound, where there is no single f0
    for snr_db to look at.
    """
    r = np.asarray(reference, dtype=float).ravel()
    o = np.asarray(output, dtype=float).ravel()
    r = r - r.mean()
    o = o - o.mean()
    if not (o @ o) or not (r @ r):
        return 0.0
    rho2 = (r @ o) ** 2 / ((r @ r) * (o @ o))
    return float(10 * np.log10(1 / max(1 - rho2, 1e-12)))


def xcorr_peak(clean, output):
    """Peak normalised cross-correlation and its lag in samples.

    For aperiodic inputs such as E8's pulse train (PLAN.md 2.4). Both records
    are mean-removed, so the result is a Pearson correlation per lag. Divide the
    lag by fs for model time.

    On a periodic input the correlation is periodic in lag, so the peak is
    ambiguous modulo the period and the triangular overlap envelope pulls the
    reported lag towards zero. Even aperiodic, the coefficient cannot reach 1 at
    a non-zero lag: overlap costs |lag| samples while the normalisation keeps
    the full energy of both records.
    """
    a = np.asarray(clean, dtype=float).ravel()
    b = np.asarray(output, dtype=float).ravel()
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt((a @ a) * (b @ b))
    if denom == 0:
        raise ValueError("cross-correlation is undefined for a constant input")
    r = correlate(b, a, mode="full") / denom
    lags = correlation_lags(b.size, a.size, mode="full")
    i = np.argmax(np.abs(r))
    return float(r[i]), int(lags[i])


def roc_auc(stat_h0, stat_h1):
    """Area under the ROC curve, from the rank-sum identity.

    0.5 means the statistic cannot separate the hypotheses. Mid-ranks handle
    ties, which matters for the binary threshold-detector output.
    """
    h0 = np.asarray(stat_h0, dtype=float).ravel()
    h1 = np.asarray(stat_h1, dtype=float).ravel()
    ranks = rankdata(np.concatenate([h0, h1]))
    return float(
        (ranks[h0.size :].sum() - h1.size * (h1.size + 1) / 2) / (h0.size * h1.size)
    )


def pd_at_pfa(stat_h0, stat_h1, pfa=0.05):
    """Detection probability at a fixed false-alarm rate.

    Threshold is the (1-pfa) quantile of the signal-absent statistic, taken with
    the "higher" rule so the realised false-alarm rate never exceeds pfa.
    """
    h0 = np.asarray(stat_h0, dtype=float).ravel()
    h1 = np.asarray(stat_h1, dtype=float).ravel()
    return float(np.mean(h1 > np.quantile(h0, 1 - pfa, method="higher")))
