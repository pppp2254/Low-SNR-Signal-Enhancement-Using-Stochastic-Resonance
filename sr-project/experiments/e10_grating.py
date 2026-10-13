"""E10 grating: stochastic resonance in two dimensions, without an image.

The foundation under E10's real images.

E2 with space in place of time: a grating s(x, y) = A cos(2 pi (kx x + ky y) / N)
lies below theta everywhere, so a noiseless threshold returns a blank frame.
Noise is added, every pixel is thresholded, and n_frames binary frames are
averaged. Quality is the spatial SNR at k from the 2-D FFT, against the
closed form for independent pixels. Synthetic, so the right answer is known.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting
from srlab.metrics import snr2d_db
from srlab.sweep import cached, load_config
from srlab.threshold import threshold_detector, threshold_snr2d_theory_db


def grating(e):
    y, x = np.mgrid[0:e["N"], 0:e["N"]]
    return e["A"] * np.cos(2 * np.pi * (e["k"][0] * x + e["k"][1] * y) / e["N"])


def main():
    cfg = load_config()
    e, seed = cfg["experiments"]["e10_grating"], cfg["seeds"]["e10_grating"]
    sigma = np.logspace(np.log10(e["sigma_min"]), np.log10(e["sigma_max"]), e["n_sigma"])
    s = grating(e)
    k, n = tuple(e["k"]), e["N"]

    def compute():
        snr, means = [], []
        for i, sg in enumerate(sigma):
            frames = threshold_detector(s.ravel(), sg, e["theta"], e["n_frames"], seed + i)[0]
            frames = frames.reshape(e["n_frames"], n, n)
            snr.append(snr2d_db(frames, k, e["n_side"], e["guard"]))
            means.append(frames.mean(axis=0))
        snr = np.array(snr)
        show = [int(np.argmin(np.abs(sigma - 0.15))), int(np.argmax(snr)), sigma.size - 1]
        return {"sigma": sigma, "snr": snr, "theory": threshold_snr2d_theory_db(s, e["theta"], sigma, k),
                "show": np.array(show), "means": np.array([means[j] for j in show], dtype=np.float32)}

    r = cached("e10_grating", dict(e, seed=seed), compute)
    sigma, snr, theory, show = r["sigma"], r["snr"], r["theory"], r["show"]
    kp, kt = int(np.argmax(snr)), int(np.argmax(theory))

    plotting.use_style()
    fig = plt.figure(figsize=(11, 6.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.15])
    panels = [("input grating, all below theta", s, -1.0, 1.0)] + [
        (f"output, sigma = {sigma[j]:.2f}", m, 0.0, None) for j, m in zip(show, r["means"])]
    for col, (title, img, lo, hi) in enumerate(panels):
        ax = fig.add_subplot(gs[0, col])
        ax.imshow(img[:64, :64], cmap="gray", vmin=lo, vmax=hi if hi is not None else max(img.max(), 1e-9),
                  interpolation="nearest")
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    ax = fig.add_subplot(gs[1, :])
    ax.plot(sigma, snr, "o", ms=4, color=plotting.COLORS["sim"],
            label=f"simulated, {e['n_frames']} frames of {n} x {n}")
    ax.plot(sigma, theory, "-", lw=1.2, color=plotting.COLORS["theory"],
            label="theory, |sum p e^(-ik.x)|^2 / sum p(1-p)")
    ax.axvline(sigma[kp], color=plotting.COLORS["optimal"], ls="--", lw=1.0,
               label=f"peak at sigma = {sigma[kp]:.2f}, {snr[kp]:.1f} dB")
    for j in show:
        ax.plot(sigma[j], snr[j], "s", ms=8, mfc="none", color=plotting.COLORS["clean"])
    ax.set_xscale("log")
    ax.set_xlabel("noise standard deviation sigma (dimensionless)")
    ax.set_ylabel(f"spatial SNR at k = {k} (dB)")
    ax.set_title(f"E10 foundation: 2-D grating, A = {e['A']} below theta = {e['theta']}; squares mark the "
                 "panels above (top-left 64 x 64 shown)", loc="left")
    ax.legend(fontsize=8)
    fig.tight_layout()
    plotting.save(fig, "e10_grating")

    ok = snr > -30
    print(f"  silent below sigma = {sigma[np.argmax(ok)]:.3f} ({int((~ok).sum())} sigma values at the floor)")
    print(f"  peak: sigma = {sigma[kp]:.3f}, {snr[kp]:.2f} dB; theory peak sigma = {sigma[kt]:.3f}, "
          f"{theory[kt]:.2f} dB")
    print(f"  ends: {snr[-1]:.2f} dB at sigma = {sigma[-1]:.2f}; "
          f"max |sim - theory| above the floor {np.abs(snr - theory)[ok].max():.2f} dB")


if __name__ == "__main__":
    main()
