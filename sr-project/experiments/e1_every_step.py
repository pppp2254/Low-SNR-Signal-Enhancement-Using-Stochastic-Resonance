"""Every solver iteration, undecimated (companion to E1).

The other figures show the recorded output, which keeps every 10th solver step
(configs/default.yaml, solver.decimate). This one runs with decimate = 1 so
nothing is thrown away, and zooms until individual iterations are larger than a
pixel.

That zoom is not a stylistic choice. A 320000-step record drawn 2000 pixels
wide is 160 steps per pixel, so the leftmost panel cannot show iterations no
matter how it is drawn. Only the rightmost one, a few hundred steps across,
actually resolves them. The right column also overlays what the decimated
pipeline would have kept, which is what every other figure in this project is
built from.

The record is shorter than E1's, 32 periods against 128, because storing every
step costs ten times the samples.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.bistable import simulate
from srlab.sweep import cached, load_config

N_PERIODS = 32
LABELS = ["too little noise", "near optimal", "too much noise"]


def find_hop(x, x_m=1.0):
    """Index of a committed well-to-well transition, for the zoom panel."""
    idx = np.flatnonzero(np.abs(x) > x_m / 2)
    if idx.size < 2:
        return x.size // 2
    jumps = np.flatnonzero(np.diff(np.sign(x[idx])) != 0)
    return int(idx[jumps[len(jumps) // 2] + 1]) if jumps.size else x.size // 2


def main():
    cfg = load_config()
    s = cfg["solver"]
    a, b, A = cfg["system"]["a"], cfg["system"]["b"], cfg["signal"]["A"]
    f0, dt, decimate = cfg["signal"]["f0"], s["dt"], s["decimate"]
    D = cfg["experiments"]["e1"]["D"]

    params = dict(a=a, b=b, A=A, f0=f0, D=D, dt=dt, n_periods=N_PERIODS,
                  discard_periods=s["discard_periods"], decimate=1, n_trials=1,
                  seed=cfg["seeds"]["e1"])

    def compute():
        # decimate = 1: every Euler-Maruyama step is kept. One trial only,
        # because full rate costs ten times the samples of the usual output.
        t, x = simulate(**params)
        print(f"  {x.shape[-1]:,} steps per trial, {x.nbytes / 1e6:.1f} MB")
        return {"t": t, "x": x[:, 0, :].astype(np.float32)}

    r = cached("e1_every_step", params, compute)
    t, x = r["t"], r["x"]
    n = t.size
    period_steps = int(round(1 / f0 / dt))

    plotting.use_style()
    fig, axes = plt.subplots(len(D), 3, figsize=(14, 8.4))
    widths = ["whole record", "one drive period", "a single hop"]

    for i, d in enumerate(D):
        hop = find_hop(x[i])
        half = period_steps // 2
        spans = [(0, n), (max(0, hop - half), min(n, hop + half)),
                 (max(0, hop - 150), min(n, hop + 150))]

        for j, (lo, hi) in enumerate(spans):
            ax = axes[i, j]
            sl = slice(lo, hi)
            steps = hi - lo
            if j == 2:
                ax.plot(t[sl], x[i][sl], "-o", ms=2.4, lw=0.8,
                        color=plotting.COLORS["sim"], label="every step")
                dec = np.arange(lo, hi)
                dec = dec[(dec % decimate) == 0]
                ax.plot(t[dec], x[i][dec], "s", ms=4.5, mfc="none",
                        color=plotting.COLORS["theory"],
                        label=f"kept by decimate = {decimate}")
                ax.legend(loc="lower right", fontsize=7)
            else:
                ax.plot(t[sl], x[i][sl], lw=0.6 if j == 0 else 0.9,
                        color=plotting.COLORS["sim"])
            if j == 0:
                ax.plot(t[sl], signals.sine(1.0, f0, t[sl]),
                        color=plotting.COLORS["clean"], ls="--", lw=1.0)
            ax.set_ylim(-2.4, 2.4)
            ax.set_title(f"D = {d}, {widths[j]}: {steps:,} steps"
                         + (f", {steps / 2000:.0f} per pixel" if j == 0 else ""),
                         loc="left", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"x, {LABELS[i]}")
            if i == len(D) - 1:
                ax.set_xlabel("time (model units)")

    fig.suptitle(f"Every solver iteration, dt = {dt}, {N_PERIODS} periods, "
                 f"{n:,} steps per trial (no decimation)", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e1_every_step")

    print(f"\n  steps recorded per trial : {n:,}")
    print(f"  the usual pipeline keeps : {n // decimate:,} (every {decimate}th)")
    print(f"  whole-record panel       : {n / 2000:.0f} steps per pixel, "
          f"individual iterations not resolvable")
    print(f"  single-hop panel         : 300 steps across, "
          f"{2000 / 300:.1f} pixels per step, every iteration drawn as a marker")


if __name__ == "__main__":
    main()
