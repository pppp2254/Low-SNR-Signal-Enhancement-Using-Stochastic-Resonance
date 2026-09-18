"""E2: a plain threshold detector shows SR too (PLAN.md 4).

y = 1 if A cos(2 pi f0 t) + noise > theta, with A = 0.7 below theta = 1. No
double well, no dynamics: noise alone carries the signal over the threshold.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.sweep import cached, chunked_snr, load_config
from srlab.threshold import threshold_detector


def main():
    cfg = load_config()
    e2, s, m = cfg["experiments"]["e2"], cfg["solver"], cfg["metrics"]
    f0, fs = cfg["signal"]["f0"], s["fs"]
    n_trials, seed = cfg["ensemble"]["n_trials"], cfg["seeds"]["e2"]

    sigma = np.logspace(np.log10(e2["sigma_min"]), np.log10(e2["sigma_max"]), e2["n_sigma"])
    t = signals.time_axis(f0, fs, s["n_periods"])
    sig = signals.sine(e2["A"], f0, t)

    params = dict(A=e2["A"], theta=e2["theta"], f0=f0, fs=fs,
                  n_periods=s["n_periods"], sigma=sigma, n_trials=n_trials,
                  seed=seed, n_side=m["n_side"], guard=m["guard"], n_boot=m["n_boot"])

    def compute():
        def make(lo, hi):
            return threshold_detector(sig, sigma[lo:hi], e2["theta"], n_trials, seed + lo)

        snr, band_lo, band_hi = chunked_snr(
            make, sigma.size, n_trials * t.size * 4, fs, f0,
            n_side=m["n_side"], guard=m["guard"], n_boot=m["n_boot"], seed=seed)
        return {"sigma": sigma, "snr": snr, "lo": band_lo, "hi": band_hi}

    r = cached("e2", params, compute)
    sigma, snr, lo, hi = r["sigma"], r["snr"], r["lo"], r["hi"]
    k = int(np.argmax(snr))

    plotting.use_style()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.fill_between(sigma, lo, hi, color=plotting.COLORS["sim"], alpha=0.2,
                    label="bootstrap 16-84%")
    ax.plot(sigma, snr, "o-", ms=3, color=plotting.COLORS["sim"], label="threshold detector")
    ax.axvline(sigma[k], color=plotting.COLORS["optimal"], ls="--", lw=1.0,
               label=f"peak at sigma = {sigma[k]:.2f}, {snr[k]:.1f} dB")
    ax.set_xscale("log")
    ax.set_xlabel("noise standard deviation sigma (dimensionless)")
    ax.set_ylabel("output SNR at f0 (dB)")
    ax.set_title(f"E2: threshold detector, A = {e2['A']} below theta = {e2['theta']}, "
                 f"{n_trials} trials", loc="left")
    ax.legend()
    fig.tight_layout()
    plotting.save(fig, "e2_threshold")

    print(f"  peak: sigma = {sigma[k]:.3f}, SNR = {snr[k]:.2f} dB "
          f"[{lo[k]:.2f}, {hi[k]:.2f}]")
    print(f"  ends: {snr[0]:.2f} dB at sigma = {sigma[0]:.2f}, "
          f"{snr[-1]:.2f} dB at sigma = {sigma[-1]:.2f}")


if __name__ == "__main__":
    main()
