"""E10: does stochastic resonance help low-light images (PLAN.md 6)?

Every image in data/picture (or the paths given on the command line) is taken
as a well-lit ground truth, darkened to 12 % brightness with mild sensor noise,
then recovered by threshold SR (approach A) and by three classical methods.
Scored against the original with PSNR and SSIM, and on its own with entropy,
mean brightness and RMS contrast. Only V of HSV is processed.

Approach B, dynamic SR after Chouhan, Jha and Biswas (2013), is not run: that
paper was not accessible (HTTP 403), and PLAN.md says to report rather than
guess its parameters. See srlab/image_sr.py.
"""

import csv
import hashlib
import pathlib
import sys
import zlib

import matplotlib.pyplot as plt
import numpy as np
from skimage import color, io, transform

from srlab import plotting
from srlab.image_sr import approach_a, baselines, darken, image_metrics
from srlab.sweep import ROOT, cached, load_config

EXTS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"}


def load_rgb(path):
    im = io.imread(path)
    if im.ndim == 2:
        im = np.stack([im] * 3, axis=-1)
    im = im[..., :3].astype(float)
    return im / 255.0 if im.max() > 1 else im


def thumb(img, side=800):
    """Display copy for the grid; a panel is ~780 px at 300 dpi. Metrics use full size."""
    k = max(img.shape[:2]) / side
    if k <= 1:
        return img
    return transform.resize(img, (round(img.shape[0] / k), round(img.shape[1] / k)),
                            anti_aliasing=True)


def with_v(rgb, v):
    """Same hue and saturation, new brightness."""
    hsv = color.rgb2hsv(rgb)
    hsv[..., 2] = np.clip(v, 0.0, 1.0)
    return color.hsv2rgb(hsv)


def process(cfg, path):
    e = cfg["experiments"]["e10"]
    # Seeded by file name, so adding or removing an image leaves the others unchanged.
    seed = cfg["seeds"]["e10"] + zlib.crc32(pathlib.Path(path).name.encode()) % 100_000
    v = color.rgb2hsv(load_rgb(path))[..., 2]
    dark = darken(v, e["darken"], e["sensor_noise"], np.random.default_rng(seed))
    sigma = np.logspace(np.log10(e["sigma_min"]), np.log10(e["sigma_max"]), e["n_sigma"])

    def compute():
        ssim, psnr = [], []
        for i, s in enumerate(sigma):
            m = image_metrics(v, approach_a(dark, s, e["n_frames"], seed + i))
            ssim.append(m["ssim"])
            psnr.append(m["psnr_db"])
        return {"sigma": sigma, "ssim": np.array(ssim), "psnr": np.array(psnr)}

    params = dict(e, image=hashlib.sha1(pathlib.Path(path).read_bytes()).hexdigest(),
                  seed=seed)
    r = cached(f"e10_{pathlib.Path(path).stem}", params, compute)
    k = int(np.argmax(r["ssim"]))
    outputs = {"dark input": dark}
    outputs.update(baselines(dark, e["gamma"]))
    outputs[f"SR, sigma {r['sigma'][k]:.3f}"] = approach_a(dark, r["sigma"][k],
                                                           e["n_frames"], seed + k)
    rows = [dict(image=pathlib.Path(path).name, method=name, **image_metrics(v, out))
            for name, out in outputs.items()]
    return v, dark, outputs, r, rows


def main():
    cfg = load_config()
    e = cfg["experiments"]["e10"]
    paths = [pathlib.Path(p) for p in sys.argv[1:]] or sorted(
        p for p in (ROOT / e["images"]).iterdir() if p.suffix.lower() in EXTS)
    if not paths:
        raise SystemExit(f"no images in {ROOT / e['images']}")

    results = [process(cfg, p) for p in paths]
    all_rows = [row for *_, rows in results for row in rows]

    out_csv = ROOT / "results" / "e10_metrics.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0]))
        w.writeheader()
        w.writerows(all_rows)
    print(f"  [csv] results/{out_csv.name}")

    plotting.use_style()
    cols = ["original"] + list(results[0][2])
    fig, axes = plt.subplots(len(paths), len(cols), figsize=(2.6 * len(cols), 2.2 * len(paths)),
                             squeeze=False)
    for row, (path, (v, dark, outputs, r, rows)) in enumerate(zip(paths, results)):
        rgb = thumb(load_rgb(path))
        panels = [("original", rgb)] + [(n, with_v(rgb, thumb(o))) for n, o in outputs.items()]
        for col, (name, img) in enumerate(panels):
            ax = axes[row, col]
            ax.imshow(np.clip(img, 0, 1))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            ssim = "" if name == "original" else \
                f"\nSSIM {next(x['ssim'] for x in rows if x['method'] == name):.3f}"
            ax.set_title(name.split(",")[0] + ssim, fontsize=8)
        axes[row, 0].set_ylabel(path.name, fontsize=8)
    fig.suptitle("E10: low-light enhancement, V channel only, darkened to "
                 f"{e['darken']:.0%}", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e10_image_grid")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    # Line style and marker both change per image, so the curves stay apart in
    # black and white and never repeat for up to 25 images.
    dashes, marks = ["-", "--", ":", "-."], ["o", "s", "^", "D", "v"]
    for i, (path, (_, _, _, r, _)) in enumerate(zip(paths, results)):
        k = int(np.argmax(r["ssim"]))
        style = dict(ls=dashes[i % 4], marker=marks[i % 5], ms=3.5,
                     color=plotting.COLORS["sim"])
        axes[0].plot(r["sigma"], r["ssim"], label=f"{path.name}, best {r['ssim'][k]:.3f}",
                     **style)
        axes[1].plot(r["sigma"], r["psnr"], label=path.name, **style)
    for ax, lab in zip(axes, ("SSIM against the original", "PSNR against the original (dB)")):
        ax.set_xscale("log")
        ax.set_xlabel("added noise sigma (V units)")
        ax.set_ylabel(lab)
        ax.legend(fontsize=7.5)
    axes[0].set_title("Threshold SR: SSIM against noise", loc="left")
    axes[1].set_title("Threshold SR: PSNR against noise", loc="left")
    fig.tight_layout()
    plotting.save(fig, "e10_sr_curve")

    print()
    for path, (_, _, _, r, rows) in zip(paths, results):
        k = int(np.argmax(r["ssim"]))
        print(f"  {path.name}: SR SSIM {r['ssim'][0]:.3f} at sigma {r['sigma'][0]} -> "
              f"{r['ssim'][k]:.3f} at {r['sigma'][k]:.4f} -> {r['ssim'][-1]:.3f} at "
              f"{r['sigma'][-1]}")
        for x in rows:
            print(f"      {x['method']:22s} SSIM {x['ssim']:.3f}  PSNR {x['psnr_db']:6.2f} dB  "
                  f"entropy {x['entropy_bits']:.2f}  mean {x['mean']:.3f}  "
                  f"contrast {x['rms_contrast']:.3f}")


if __name__ == "__main__":
    main()
