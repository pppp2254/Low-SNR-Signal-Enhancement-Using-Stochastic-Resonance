"""Low-light image extension, E10 (PLAN.md 6).

Works on the V channel of HSV only; hue and saturation are left untouched.

Approach A, threshold SR: add Gaussian noise to the dark V channel, threshold
every pixel at theta, and average n_frames binary frames. theta is set at the
dark image's maximum, so with no noise not a single pixel fires; that is E2's
sub-threshold setup in two dimensions. The average of the frames is the firing
probability per pixel, a smooth tone curve whose shape the noise controls.

Denoising afterwards (e10_denoise.py): TV or non-local means removes the
binomial grain but keeps the mean, so the image does not go dark again.
sr_passes repeats SR plus TV on its own output.

Approach B, dynamic SR after Chouhan, Jha and Biswas (2013), is not here. That
paper returned 403 when fetched, and PLAN.md says to report that rather than
guess its parameters and stopping rule.
"""

import numpy as np
from skimage import exposure, measure, metrics, restoration

from .threshold import threshold_detector

__all__ = ["darken", "approach_a", "baselines", "image_metrics", "noise_sigma", "denoise",
           "sr_passes"]


def darken(v, scale, noise_sigma, rng):
    """A dark test image from a well-lit V channel, so a ground truth exists.

    Scaled, given mild Gaussian sensor noise, clipped and quantised to 8 bits.
    Quantisation matters: at scale 0.12 the whole image fits in about 30 grey
    levels, which is part of what makes a dark image hard.
    """
    dark = scale * np.asarray(v, dtype=float) + noise_sigma * rng.standard_normal(np.shape(v))
    return np.round(np.clip(dark, 0.0, 1.0) * 255) / 255


def approach_a(v_dark, sigma, n_frames, seed, theta=None):
    """Threshold SR: average of n_frames noisy binary frames, stretched to [0, 1].

    The raw average is a firing probability that never exceeds about 0.5 when
    theta sits at the image maximum, so it is stretched linearly to the full
    range. That stretch is the only post-processing and is applied the same
    way at every sigma.
    """
    v_dark = np.asarray(v_dark, dtype=float)
    theta = float(v_dark.max()) if theta is None else theta
    frames = threshold_detector(v_dark.ravel(), sigma, theta, n_frames, seed)[0]
    y = frames.mean(axis=0).reshape(v_dark.shape).astype(float)
    top = y.max()
    return y / top if top > 0 else y


def baselines(v_dark, gamma):
    """The classical methods PLAN.md 6 compares against."""
    v_dark = np.asarray(v_dark, dtype=float)
    return {
        "histogram eq": exposure.equalize_hist(v_dark),
        "CLAHE": exposure.equalize_adapthist(v_dark),
        f"gamma {gamma}": v_dark ** gamma,
    }


def image_metrics(reference, output):
    """PSNR and SSIM against the original; entropy, brightness, RMS contrast alone."""
    ref = np.asarray(reference, dtype=float)
    out = np.clip(np.asarray(output, dtype=float), 0.0, 1.0)
    return {
        "psnr_db": float(metrics.peak_signal_noise_ratio(ref, out, data_range=1.0)),
        "ssim": float(metrics.structural_similarity(ref, out, data_range=1.0)),
        "entropy_bits": float(measure.shannon_entropy(np.round(out * 255))),
        "mean": float(out.mean()),
        "rms_contrast": float(out.std()),
    }


def noise_sigma(x):
    """Blind noise level from the image alone: MAD of diagonal pixel differences.

    A difference of two independent noisy pixels has std sqrt(2) sigma, and
    0.6745 turns a median absolute value into a std. Edges and fine texture
    count as noise too, so a busy image gets an overestimate.
    """
    x = np.asarray(x, dtype=float)
    return float(np.median(np.abs(x[1:, 1:] - x[:-1, :-1])) / 0.6745 / np.sqrt(2))


def denoise(x, method):
    """TV (Chambolle) or non-local means, strength from noise_sigma, no ground truth.

    Returned unchanged when the estimate is zero, which happens when over half
    the image is flat (a silhouette): no measurable noise, nothing to remove.
    """
    s = noise_sigma(x)
    if s == 0:
        return np.asarray(x, dtype=float)
    if method == "TV":
        return restoration.denoise_tv_chambolle(x, weight=s)
    if method == "NLM":
        # h = 0.8 sigma is scikit-image's suggestion for fast mode with sigma given.
        return restoration.denoise_nl_means(x, h=0.8 * s, sigma=s, patch_size=5,
                                            patch_distance=6, fast_mode=True)
    raise ValueError(f"unknown method {method!r}")


def sr_passes(v_dark, sigma, n_frames, seed, passes, theta="max"):
    """Threshold SR then TV, repeated on its own output; yields each pass.

    Pass 1 is approach_a at theta = max, exactly as E10. Later passes threshold
    the previous output at its "max" or its "median", with sigma scaled by the
    output's range so the noise keeps the same size relative to the image.
    """
    x, span = np.asarray(v_dark, dtype=float), np.ptp(v_dark)
    for i in range(passes):
        th = None if theta == "max" or i == 0 else float(np.median(x))
        x = denoise(approach_a(x, sigma * np.ptp(x) / span, n_frames, seed + 1000 * i,
                               theta=th), "TV")
        yield x
