"""E10 frame count: does threshold-SR grain fall as 1/sqrt(N)?

At each image's best sigma from e10_image.py, count firings over N frames
and score two things per N. Grain: two independent runs differ only by
sampling noise, so rms(f1 - f2)/sqrt(2) is the per-pixel noise of one run,
against theory sqrt(mean p(1-p) / N) with p = P(v + sigma xi > theta).
Quality: SSIM of run 1 against the original, stretched as in approach_a,
against the N -> infinity limit, the expected image p / max p.

The N checkpoints are nested (one run of max N frames, read off on the way),
so points along N share frames. Run 1 uses E10's seed, so at N = 50 it is the
same output as the SR column of e10_metrics.csv.
"""

import hashlib
import pathlib
import zlib

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc
from skimage import color

from e10_image import EXTS, load_rgb
from srlab import plotting
from srlab.image_sr import darken, image_metrics
from srlab.sweep import ROOT, cached, load_config


def process(cfg, path):
    e = cfg["experiments"]["e10"]
    seed = cfg["seeds"]["e10"] + zlib.crc32(path.name.encode()) % 100_000
    best = np.load(ROOT / "results" / f"e10_{path.stem}.npz")
    k = int(np.argmax(best["ssim"]))
    sigma = float(best["sigma"][k])
    checkpoints = cfg["experiments"]["e10_frames"]["frames"]

    def compute():
        v = color.rgb2hsv(load_rgb(path))[..., 2]
        dark = darken(v, e["darken"], e["sensor_noise"], np.random.default_rng(seed))
        theta = dark.max()
        p = 0.5 * erfc((theta - dark) / (sigma * np.sqrt(2)))
        rngs = [np.random.default_rng(seed + k), np.random.default_rng(seed + 10_000)]
        counts = [np.zeros(dark.shape, dtype=np.uint16) for _ in rngs]
        grain, ssim = [], []
        # One frame in memory at a time; the full stack is 2.6 GB at N = 400.
        for n in range(1, max(checkpoints) + 1):
            for c, rng in zip(counts, rngs):
                c += dark + sigma * rng.standard_normal(dark.shape) > theta
            if n in checkpoints:
                f1, f2 = counts[0] / n, counts[1] / n
                grain.append(np.sqrt(np.mean((f1 - f2) ** 2) / 2))
                ssim.append(image_metrics(v, f1 / f1.max())["ssim"])
        return {"n": np.array(checkpoints), "sigma": sigma, "grain": np.array(grain),
                "theory": np.sqrt(np.mean(p * (1 - p)) / np.array(checkpoints)),
                "ssim": np.array(ssim), "ssim_limit": image_metrics(v, p / p.max())["ssim"]}

    params = dict(e, frames=checkpoints, image=hashlib.sha1(path.read_bytes()).hexdigest(),
                  seed=seed, sigma=sigma)
    return cached(f"e10_frames_{path.stem}", params, compute)


def main():
    cfg = load_config()
    e = cfg["experiments"]["e10"]
    paths = sorted(p for p in (ROOT / e["images"]).iterdir() if p.suffix.lower() in EXTS)
    results = [process(cfg, p) for p in paths]

    plotting.use_style()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    dashes, marks = ["-", "--", ":", "-."], ["o", "s", "^", "D", "v", "P", "X"]
    for i, (path, r) in enumerate(zip(paths, results)):
        slope = np.polyfit(np.log(r["n"]), np.log(r["grain"]), 1)[0]
        style = dict(ls="none", marker=marks[i % 7], ms=4, color=plotting.COLORS["sim"])
        axes[0].plot(r["n"], r["grain"], label=f"{path.name}, slope {slope:.3f}", **style)
        axes[0].plot(r["n"], r["theory"], ls=dashes[i % 4], lw=0.8,
                     color=plotting.COLORS["theory"])
        axes[1].plot(r["n"], r["ssim"], label=f"{path.name}, limit {r['ssim_limit']:.3f}",
                     **dict(style, ls=dashes[i % 4]))
        axes[1].axhline(r["ssim_limit"], ls=dashes[i % 4], lw=0.8, color=plotting.COLORS["clean"])
        print(f"  {path.name}: sigma {r['sigma']:.4f}, grain slope {slope:.3f}, "
              f"measured/theory {np.round(r['grain'] / r['theory'], 3)}, "
              f"SSIM {np.round(r['ssim'], 3)} -> limit {r['ssim_limit']:.3f}")
    axes[0].plot([], [], color=plotting.COLORS["theory"], lw=0.8,
                 label="theory sqrt(mean p(1-p) / N)")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("grain, rms of one run's firing fraction")
    axes[0].set_title("Grain against frame count", loc="left")
    axes[1].set_xscale("log")
    axes[1].set_ylabel("SSIM against the original")
    axes[1].set_title("SSIM against frame count, grey lines N -> infinity", loc="left")
    for ax in axes:
        ax.set_xlabel("frames averaged, N")
        ax.legend(fontsize=7)
    fig.tight_layout()
    plotting.save(fig, "e10_frames")


if __name__ == "__main__":
    main()
