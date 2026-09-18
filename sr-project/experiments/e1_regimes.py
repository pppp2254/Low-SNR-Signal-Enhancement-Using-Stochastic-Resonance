"""E1: what SR looks like in time (PLAN.md 4).

Three noise regimes, D = 0.02, 0.12, 1.0: too little to hop, matched to the
signal, too much to follow it. Time series with the clean sine overlaid, and
the spectrum beside each.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.bistable import simulate
from srlab.metrics import snr_db
from srlab.sweep import cached, load_config

LABELS = ["too little noise", "near optimal", "too much noise"]
N_TRIALS = 20
SHOW_PERIODS = 6


def main():
    cfg = load_config()
    a, b = cfg["system"]["a"], cfg["system"]["b"]
    A, f0 = cfg["signal"]["A"], cfg["signal"]["f0"]
    fs, s = cfg["solver"]["fs"], cfg["solver"]
    D = cfg["experiments"]["e1"]["D"]

    params = dict(a=a, b=b, A=A, f0=f0, D=D, dt=s["dt"], n_periods=s["n_periods"],
                  discard_periods=s["discard_periods"], decimate=s["decimate"],
                  n_trials=N_TRIALS, seed=cfg["seeds"]["e1"], method=s["method"],
                  show_periods=SHOW_PERIODS)

    def compute():
        # Store only what the figure needs. Keeping all 20 raw traces per D made
        # a 27 MB npz; the panels use one trace and the averaged spectrum.
        t, x = simulate(**{k: v for k, v in params.items() if k != "show_periods"})
        spec = ((2.0 / (fs * t.size)) * np.abs(np.fft.rfft(x, axis=-1)) ** 2).mean(axis=1)
        freq = np.fft.rfftfreq(t.size, 1 / fs)
        band = freq < 10 * f0
        snr, _, _ = snr_db(x, fs, f0, cfg["metrics"]["n_side"], cfg["metrics"]["guard"])
        n_show = int(SHOW_PERIODS / f0 * fs)
        return {"t": t[:n_show], "trace": x[:, 0, :n_show], "freq": freq[band],
                "spec": spec[:, band], "snr": snr,
                "k0": np.array(round(f0 * t.size / fs))}

    r = cached("e1", params, compute)
    t, trace, freq, spec, snr = r["t"], r["trace"], r["freq"], r["spec"], r["snr"]
    k0 = int(r["k0"])

    plotting.use_style()
    fig, axes = plt.subplots(3, 2, figsize=(9, 7),
                             gridspec_kw={"width_ratios": [1.9, 1]})
    clean = signals.sine(1.0, f0, t)

    for i, d in enumerate(D):
        ax = axes[i, 0]
        ax.plot(t, trace[i], color=plotting.COLORS["sim"], lw=0.9,
                label="output x(t)")
        ax.plot(t, clean, color=plotting.COLORS["clean"], lw=1.6, ls="--",
                label="clean signal (scaled)")
        ax.set_ylim(-2.2, 2.2)
        ax.set_ylabel("x (dimensionless)")
        ax.set_title(f"D = {d}, {LABELS[i]}: SNR = {snr[i]:.1f} dB", loc="left")
        if i == 0:
            ax.legend(loc="upper right", ncol=2)

        ax = axes[i, 1]
        band = freq > 0
        # Marker under the data: drawn on top it hides the very spike it marks.
        ax.axvline(f0, color=plotting.COLORS["optimal"], lw=1.0, ls=":", zorder=0)
        ax.semilogy(freq[band], spec[i][band], color=plotting.COLORS["sim"], lw=0.9,
                    zorder=2)
        ax.plot(f0, spec[i][k0], "o", ms=4, color=plotting.COLORS["optimal"], zorder=3,
                label=f"f0 bin: {10 * np.log10(spec[i][k0] / np.median(spec[i][band])):.0f} dB over background")
        ax.set_ylabel("PSD (1/Hz)")
        ax.set_title("spectrum", loc="left")
        ax.legend(loc="upper right")

    axes[-1, 0].set_xlabel("time (model units)")
    axes[-1, 1].set_xlabel("frequency (model units)")
    fig.suptitle(
        f"E1: three noise regimes, a = b = 1, A = {A} (sub-threshold), f0 = {f0}",
        y=0.98)
    fig.tight_layout()
    plotting.save(fig, "e1_regimes")

    for d, s_ in zip(D, snr):
        print(f"  D = {d:<5} SNR = {s_:6.2f} dB")


if __name__ == "__main__":
    main()
