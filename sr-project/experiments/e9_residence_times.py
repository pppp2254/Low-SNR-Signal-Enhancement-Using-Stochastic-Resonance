"""E9: is the switching really synchronised with the drive (PLAN.md 4)?

Residence times are how long the particle stays in one well between committed
transitions. If the hopping is locked to the drive rather than merely frequent,
those times cluster at odd multiples of half the drive period: the particle
waits for the potential to tilt its way, hops, and waits again.

This is independent evidence for what E1 shows only by eye and E3 shows only as
a spectral ratio, and it does not use the SNR estimator at all.

The structure is sharpest below the SNR optimum, not at it. At D = 0.04 the
particle hops rarely enough to skip whole cycles, so the histogram shows clean
separate peaks at T/2, 3T/2 and 5T/2. At the optimum it hops almost every half
period, so nearly all the weight collapses into the first peak. Both are
synchronisation; the multi-peak picture just needs a slower hopping rate to be
visible.

The 1/r_K printed on each panel is the Kramers time for the *undriven* system.
At low noise it is wildly longer than the measured mean, because with A = 0.3
the drive periodically lowers the barrier and does most of the work. The gap is
the drive's contribution, not an error.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting
from srlab.bistable import residence_times, simulate
from srlab.sweep import cached, load_config
from srlab.theory import barrier_height, kramers_rate

LABELS = ["too little noise", "near optimal", "too much noise"]


def main():
    cfg = load_config()
    s, e9 = cfg["solver"], cfg["experiments"]["e9"]
    a, b, A = cfg["system"]["a"], cfg["system"]["b"], cfg["signal"]["A"]
    f0, dt = cfg["signal"]["f0"], s["dt"]
    period = 1.0 / f0
    D = e9["D"]
    d_opt = barrier_height(a, b) / 2

    params = dict(a=a, b=b, A=A, f0=f0, D=D, dt=dt, n_periods=e9["n_periods"],
                  discard_periods=s["discard_periods"], decimate=s["decimate"],
                  n_trials=e9["n_trials"], seed=cfg["seeds"]["e9"])

    def compute():
        print(f"  {len(D)} noise levels, {e9['n_trials']} trials, "
              f"{e9['n_periods']} periods each")
        t, x = simulate(**params)
        dt_sample = t[1] - t[0]
        out = {"dt_sample": np.array(dt_sample)}
        for i, d in enumerate(D):
            times = np.concatenate([residence_times(trial, dt_sample)
                                    for trial in x[i]])
            out[f"times{i}"] = times.astype(np.float32)
            print(f"    D = {d}: {times.size} residences, "
                  f"mean {times.mean():.1f} time units")
        return out

    res = cached("e9", params, compute)

    plotting.use_style()
    fig, axes = plt.subplots(1, len(D), figsize=(4.6 * len(D), 4.2), sharey=True)
    summary = []
    for i, (d, ax) in enumerate(zip(D, axes)):
        times = res[f"times{i}"]
        upper = 5.5 * period
        # How concentrated the times are near odd half periods against even
        # ones, which is what synchronisation should favour.
        phase = (times % period) / period
        near_odd = float(np.mean(np.abs(phase - 0.5) < 0.125))
        near_even = float(np.mean(np.minimum(phase, 1 - phase) < 0.125))
        ax.hist(times[times < upper], bins=e9["n_bins"], range=(0, upper),
                density=True, color=plotting.COLORS["sim"], alpha=0.85)
        # Log density: the 3T/2 and 5T/2 peaks are one to two orders below the
        # first and are invisible on a linear axis.
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 0.2)
        for n in (1, 3, 5, 7, 9):
            x_line = n * period / 2
            if x_line < upper:
                ax.axvline(x_line, color=plotting.COLORS["optimal"], ls="--", lw=1.0)
        summary.append((d, times.size, float(times.mean()), near_odd, near_even))

        ax.set_xlabel("residence time (model units)")
        ax.set_title(f"D = {d}, {LABELS[i]}\n"
                     f"mean {times.mean():.0f}, undriven 1/r_K = "
                     f"{1 / kramers_rate(a, b, d):.0f}, "
                     f"odd/even {near_odd / near_even if near_even else float('inf'):.0f}",
                     loc="left")
        if i == 0:
            ax.set_ylabel("probability density (1/time)")
    axes[0].plot([], [], color=plotting.COLORS["optimal"], ls="--",
                 label="odd multiples of T/2")
    axes[0].legend(loc="upper right")

    fig.suptitle(f"E9: residence times, drive period T = {period:.0f}, "
                 f"A = {A}, optimum at D = {d_opt}", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e9_residence_times")

    print()
    for d, n, mean, odd, even in summary:
        print(f"  D = {d:<6} {n:6d} residences, mean {mean:7.1f}, "
              f"near odd T/2 {odd:.3f}, near even {even:.3f}, "
              f"ratio {odd / even if even else float('inf'):.2f}")


if __name__ == "__main__":
    main()
