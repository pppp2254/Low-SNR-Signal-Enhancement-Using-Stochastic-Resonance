"""E10 denoising: does removing SR's grain afterwards change the verdict?

At each image's best sigma from e10_image.py: SR as in E10, then TV or
non-local means; SR plus TV repeated on its own output, thresholding at the
max or at the median; and, for a fair test, every baseline denoised the same
two ways, before and after enhancing. Denoiser strength comes from a blind
noise estimate, never from the original. Only metrics are cached, no pixels.
"""

import hashlib
import zlib

import matplotlib.pyplot as plt
import numpy as np
from skimage import color

from e10_image import EXTS, load_rgb
from srlab import plotting
from srlab.image_sr import (approach_a, baselines, darken, denoise, image_metrics,
                            sr_passes)
from srlab.sweep import ROOT, cached, load_config

CROP_IMAGE = "black white river.jpeg"


def setup(cfg, path):
    e = cfg["experiments"]["e10"]
    seed = cfg["seeds"]["e10"] + zlib.crc32(path.name.encode()) % 100_000
    best = np.load(ROOT / "results" / f"e10_{path.stem}.npz")
    k = int(np.argmax(best["ssim"]))
    return seed, k, float(best["sigma"][k])


def outputs(cfg, path):
    """The original V, and every method's output, one at a time to save memory."""
    e = cfg["experiments"]["e10"]
    seed, k, sigma = setup(cfg, path)
    v = color.rgb2hsv(load_rgb(path))[..., 2]
    dark = darken(v, e["darken"], e["sensor_noise"], np.random.default_rng(seed))

    def gen():
        # seed + k: the same SR output as E10's table.
        sr = approach_a(dark, sigma, e["n_frames"], seed + k)
        yield "SR", sr
        yield "SR then TV", denoise(sr, "TV")
        yield "SR then NLM", denoise(sr, "NLM")
        for mode in ("median", "max"):
            for i, x in enumerate(sr_passes(dark, sigma, e["n_frames"], seed + k,
                                            cfg["experiments"]["e10_denoise"]["passes"],
                                            mode), 1):
                yield f"pass {i}, theta {mode}", x
        yield from baselines(dark, e["gamma"]).items()
        for how in ("TV", "NLM"):
            for name, x in baselines(dark, e["gamma"]).items():
                yield f"{name} then {how}", denoise(x, how)
            for name, x in baselines(denoise(dark, how), e["gamma"]).items():
                yield f"{how} then {name}", x

    return v, gen()


def process(cfg, path):
    seed, k, sigma = setup(cfg, path)

    def compute():
        v, gen = outputs(cfg, path)
        names, ssim, psnr, mean = [], [], [], []
        for name, x in gen:
            m = image_metrics(v, x)
            names.append(name)
            ssim.append(m["ssim"])
            psnr.append(m["psnr_db"])
            mean.append(m["mean"])
        return {"names": np.array(names), "ssim": np.array(ssim), "psnr": np.array(psnr),
                "mean": np.array(mean), "v_mean": float(v.mean())}

    params = dict(cfg["experiments"]["e10"], denoise=cfg["experiments"]["e10_denoise"],
                  image=hashlib.sha1(path.read_bytes()).hexdigest(), seed=seed, sigma=sigma)
    r = cached(f"e10_denoise_{path.stem}", params, compute)
    return {str(n): (s, p, m) for n, s, p, m in zip(r["names"], r["ssim"], r["psnr"], r["mean"])}, \
        float(r["v_mean"])


def summary(scores, gamma):
    """Best of each family by SSIM: (name, ssim)."""
    raw = ("histogram eq", "CLAHE", f"gamma {gamma}")
    fam = {
        "SR": ["SR"],
        "SR denoised": ["SR then TV", "SR then NLM"],
        "SR passes, theta median": [n for n in scores if n.endswith("theta median")],
        "baseline": list(raw),
        "baseline denoised": [n for n in scores if " then " in n and not n.startswith(("SR", "pass"))],
    }
    return {f: max(((n, scores[n][0]) for n in names), key=lambda t: t[1])
            for f, names in fam.items()}


def main():
    cfg = load_config()
    e = cfg["experiments"]["e10"]
    passes = cfg["experiments"]["e10_denoise"]["passes"]
    paths = sorted(p for p in (ROOT / e["images"]).iterdir() if p.suffix.lower() in EXTS)
    results = [process(cfg, p) for p in paths]
    sums = [summary(s, e["gamma"]) for s, _ in results]

    for path, (scores, v_mean), sm in zip(paths, results, sums):
        print(f"  {path.name} (original mean {v_mean:.3f})")
        for fam, (name, s) in sm.items():
            print(f"      {fam:24s} {name:24s} SSIM {s:.3f}  mean {scores[name][2]:.3f}")

    plotting.use_style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.4, 1, 1]})
    c = plotting.COLORS
    # Marker shape and fill carry the method, so nothing rests on colour alone.
    styles = {
        "SR": dict(marker="o", mfc="white", color=c["sim"]),
        "SR denoised": dict(marker="o", color=c["sim"]),
        "SR passes, theta median": dict(marker="D", color=c["sim"]),
        "baseline": dict(marker="s", mfc="white", color=c["theory"]),
        "baseline denoised": dict(marker="s", color=c["theory"]),
    }
    y = np.arange(len(paths))
    ax = axes[0]
    for fam, st in styles.items():
        ax.plot([sm[fam][1] for sm in sums], y, ls="none", ms=7, label=fam, **st)
    ax.set_yticks(y, [p.name for p in paths])
    ax.invert_yaxis()
    ax.set_xlabel("SSIM against the original")
    ax.set_xlim(0, 1)
    ax.set_title("Best of each method family", loc="left")

    n = np.arange(1, passes + 1)
    for ax, key, lab in ((axes[1], 0, "SSIM against the original"),
                         (axes[2], 2, "mean brightness / original mean")):
        for mode, col, ls in (("median", c["sim"], "-"), ("max", c["clean"], "--")):
            rows = np.array([[s[f"pass {i}, theta {mode}"][key] / (vm if key == 2 else 1)
                              for i in n] for s, vm in results])
            for row in rows:
                ax.plot(n, row, color=col, ls=ls, lw=0.7, alpha=0.6)
            ax.plot(n, rows.mean(axis=0), color=col, ls=ls, lw=2.2, marker="o", ms=5,
                    label=f"theta at {mode}: mean of images (thin: each image)")
        ax.set_xticks(n)
        ax.set_xlabel("pass of SR then TV")
        ax.set_ylabel(lab)
    axes[2].axhline(1.0, color=c["clean"], lw=0.8, ls=":")
    axes[1].set_title("Repeating SR then TV: quality", loc="left")
    axes[2].set_title("Repeating SR then TV: brightness", loc="left")
    # Legends below the axes, where no data can sit under them.
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower left",
               bbox_to_anchor=(0.04, 0.0), ncol=3, fontsize=8)
    fig.legend(*axes[1].get_legend_handles_labels(), loc="lower right",
               bbox_to_anchor=(0.99, 0.0), ncol=1, fontsize=8)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    plotting.save(fig, "e10_denoise")

    crop = ROOT / e["images"] / CROP_IMAGE
    if crop in paths:
        i = paths.index(crop)
        sm = sums[i]
        want = ["SR", "SR then TV", sm["SR passes, theta median"][0], sm["baseline"][0],
                sm["baseline denoised"][0]]
        v, gen = outputs(cfg, crop)
        imgs = {name: x for name, x in gen if name in want}
        h, w = v.shape
        sl = (slice(max(h // 2 - 110, 0), h // 2 + 110), slice(max(w // 2 - 160, 0), w // 2 + 160))
        fig, axes = plt.subplots(2, 3, figsize=(12, 7.4), layout="constrained")
        for ax, name in zip(axes.ravel(), ["original"] + want):
            x = v if name == "original" else imgs[name]
            ax.imshow(np.clip(x[sl], 0, 1), cmap="gray", vmin=0, vmax=1, interpolation="nearest")
            ax.set_title(name if name == "original" else
                         f"{name}\nSSIM {results[i][0][name][0]:.3f}, mean {results[i][0][name][2]:.2f}",
                         fontsize=9)
            ax.axis("off")
        fig.suptitle(f"{CROP_IMAGE}, centre crop, V channel")
        plotting.save(fig, "e10_denoise_crops")


if __name__ == "__main__":
    main()
