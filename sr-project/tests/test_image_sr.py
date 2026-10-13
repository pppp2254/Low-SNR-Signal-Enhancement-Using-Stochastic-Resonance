"""E10: threshold SR on images (PLAN.md 6)."""
import numpy as np
import pytest

from srlab.image_sr import (approach_a, baselines, darken, denoise, image_metrics,
                            noise_sigma, sr_passes)


def ramp_image():
    """A smooth horizontal ramp with a bright square: tones plus an edge."""
    v = np.tile(np.linspace(0.2, 1.0, 96), (64, 1))
    v[20:44, 30:60] = 1.0
    return v


def test_darken_scales_and_quantises():
    v = ramp_image()
    dark = darken(v, 0.12, 0.0, np.random.default_rng(0))
    assert dark.max() == pytest.approx(0.12, abs=1 / 255)
    assert np.allclose(dark * 255, np.round(dark * 255))


def test_no_noise_means_no_output():
    """theta at the image maximum: sub-threshold everywhere, like E2."""
    dark = darken(ramp_image(), 0.12, 0.0, np.random.default_rng(0))
    assert approach_a(dark, 0.0, 10, seed=1).max() == 0.0


def test_noise_has_an_optimum():
    """SSIM rises from nothing, peaks, and falls: stochastic resonance in 2-D."""
    v = ramp_image()
    dark = darken(v, 0.12, 0.002, np.random.default_rng(0))
    ssim = [image_metrics(v, approach_a(dark, s, 30, seed=2))["ssim"]
            for s in (0.001, 0.03, 3.0)]
    assert ssim[1] > ssim[0] and ssim[1] > ssim[2]


def test_baselines_keep_shape_and_range():
    dark = darken(ramp_image(), 0.12, 0.002, np.random.default_rng(0))
    for name, out in baselines(dark, 0.45).items():
        assert out.shape == dark.shape, name
        assert 0.0 <= out.min() and out.max() <= 1.0 + 1e-9, name


def test_noise_sigma_recovers_known_noise():
    x = 0.5 + 0.05 * np.random.default_rng(0).standard_normal((200, 200))
    assert noise_sigma(x) == pytest.approx(0.05, rel=0.05)


def test_denoise_removes_grain_but_keeps_brightness():
    """The point of denoising after SR: less grain, same mean, so not dark again."""
    dark = darken(ramp_image(), 0.12, 0.0, np.random.default_rng(0))
    sr = approach_a(dark, 0.05, 20, seed=2)
    for method in ("TV", "NLM"):
        out = denoise(sr, method)
        assert noise_sigma(out) < 0.5 * noise_sigma(sr)
        assert out.mean() == pytest.approx(sr.mean(), rel=0.02)


def test_first_pass_is_e10_then_tv_and_max_darkens():
    dark = darken(ramp_image(), 0.12, 0.0, np.random.default_rng(0))
    by_max = list(sr_passes(dark, 0.05, 20, seed=3, passes=3, theta="max"))
    by_median = list(sr_passes(dark, 0.05, 20, seed=3, passes=3, theta="median"))
    first = denoise(approach_a(dark, 0.05, 20, seed=3), "TV")
    assert np.array_equal(by_max[0], first) and np.array_equal(by_median[0], first)
    assert by_max[-1].mean() < by_max[0].mean() < by_median[1].mean()


def test_denoise_leaves_a_flat_image_alone():
    """Zero estimated noise used to divide by zero inside TV."""
    x = np.zeros((32, 32))
    x[:4] = 1.0
    for method in ("TV", "NLM"):
        assert np.array_equal(denoise(x, method), x)
