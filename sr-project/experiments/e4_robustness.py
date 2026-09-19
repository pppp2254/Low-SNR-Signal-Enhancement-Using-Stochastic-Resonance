"""E4: how the optimum moves with the signal (PLAN.md 4).

SNR against D for three amplitudes and four frequencies, plus a heatmap over
(D, f0). The claim being tested is that the optimal D rises with f0, which is
what time-scale matching predicts: the hopping rate has to keep up with a
faster drive, and that needs more noise.

Bandwidth caveat. snr_db is a per-bin ratio and the bin width is
df = fs / n_samples = f0 / n_periods, so holding n_periods fixed means each f0
is measured in a different bandwidth. Curves for different f0 are therefore
offset from each other by 10 log10 of the bandwidth ratio. Compare peak
positions across f0, not peak heights.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.bistable import simulate
from srlab.sweep import cached, chunked_snr, d_grid, interior_peak, load_config
from srlab.theory import barrier_height


def snr_vs_D(cfg, name, D, A, f0, n_trials, n_periods, discard, seed):
    s, m = cfg["solver"], cfg["metrics"]
    fs = s["fs"]
    n_samples = int(round(n_periods / f0 * fs))
    params = dict(a=cfg["system"]["a"], b=cfg["system"]["b"], A=A, f0=f0, D=D,
                  dt=s["dt"], n_periods=n_periods, discard_periods=discard,
                  decimate=s["decimate"], n_trials=n_trials, seed=seed,
                  n_side=m["n_side"], guard=m["guard"], n_boot=m["n_boot"])

    def compute():
        def make(lo, hi):
            _, x = simulate(a=params["a"], b=params["b"], A=A, f0=f0, D=D[lo:hi],
                            dt=s["dt"], n_periods=n_periods, discard_periods=discard,
                            decimate=s["decimate"], n_trials=n_trials, seed=seed + lo)
            return x

        snr, lo, hi = chunked_snr(make, D.size, n_trials * n_samples * 4, fs, f0,
                                  n_side=m["n_side"], guard=m["guard"],
                                  n_boot=m["n_boot"], seed=seed)
        return {"D": D, "snr": snr, "lo": lo, "hi": hi}

    print(f"  {name}: A={A}, f0={f0}, {n_trials} trials, {n_periods} periods")
    return cached(name, params, compute)


def main():
    cfg = load_config()
    s, e4 = cfg["solver"], cfg["experiments"]["e4"]
    D = d_grid(cfg)
    seed = cfg["seeds"]["e4"]
    n_trials, f0_base = cfg["ensemble"]["n_trials"], cfg["signal"]["f0"]
    d_opt = barrier_height(cfg["system"]["a"], cfg["system"]["b"]) / 2

    amps = cfg["sweeps"]["A"]
    a_runs = [snr_vs_D(cfg, f"e4_A{A}", D, A, f0_base, n_trials, s["n_periods"],
                       s["discard_periods"], seed + 100 * i)
              for i, A in enumerate(amps)]

    freqs = cfg["sweeps"]["f0"]
    f_runs = [snr_vs_D(cfg, f"e4_f{f}", D, cfg["signal"]["A"], f, n_trials,
                       s["n_periods"], s["discard_periods"], seed + 200 + 100 * i)
              for i, f in enumerate(freqs)]

    # Snap each log-spaced f0 so the record stays a whole number of periods.
    hm_f0 = np.array([signals.snap_f0(f, s["fs"], e4["heatmap_n_periods"])
                      for f in np.logspace(np.log10(e4["heatmap_f0_min"]),
                                           np.log10(e4["heatmap_f0_max"]),
                                           e4["heatmap_n_f0"])])
    hm_D = np.logspace(np.log10(cfg["noise"]["D_min"]),
                       np.log10(cfg["noise"]["D_max"]), e4["heatmap_n_D"])
    hm = [snr_vs_D(cfg, f"e4_hm{i}", hm_D, cfg["signal"]["A"], f,
                   e4["heatmap_n_trials"], e4["heatmap_n_periods"],
                   e4["heatmap_discard"], seed + 400 + 10 * i)
          for i, f in enumerate(hm_f0)]
    grid = np.array([r["snr"] for r in hm])

    plotting.use_style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    ax = axes[0]
    for A, r, c in zip(amps, a_runs, ("alt1", "sim", "alt2")):
        ax.fill_between(r["D"], r["lo"], r["hi"], color=plotting.COLORS[c], alpha=0.18)
        k = interior_peak(r["snr"])
        ax.plot(r["D"], r["snr"], "o-", ms=3, color=plotting.COLORS[c],
                label=f"A = {A}, SR peak D = {r['D'][k]:.3f}")
        ax.plot(r["D"][k], r["snr"][k], "*", ms=11, color=plotting.COLORS[c], zorder=5)
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("output SNR at f0 (dB)")
    ax.set_title(f"SNR against D for three amplitudes, f0 = {f0_base}", loc="left")
    ax.legend()

    ax = axes[1]
    for f, r, c in zip(freqs, f_runs, ("alt1", "sim", "alt2", "theory")):
        k = interior_peak(r["snr"])
        ax.plot(r["D"], r["snr"], "o-", ms=3, color=plotting.COLORS[c],
                label=f"f0 = {f}, SR peak D = {r['D'][k]:.3f}")
        ax.plot(r["D"][k], r["snr"][k], "*", ms=11, color=plotting.COLORS[c], zorder=5)
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("output SNR at f0 (dB)")
    ax.set_title("SNR against D for four frequencies (bandwidth differs, see docstring)",
                 loc="left")
    ax.legend()

    ax = axes[2]
    mesh = ax.pcolormesh(hm_D, hm_f0, grid, shading="nearest", cmap="viridis")
    peak_D = np.array([hm_D[interior_peak(row) or int(np.argmax(row))] for row in grid])
    ax.plot(peak_D, hm_f0, "o-", color="white", ms=4, lw=1.2, label="peak D per f0")
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("signal frequency f0 (dimensionless)")
    ax.set_title(f"SNR over (D, f0), {e4['heatmap_n_periods']} periods, "
                 f"{e4['heatmap_n_trials']} trials", loc="left")
    ax.legend(loc="upper left")
    fig.colorbar(mesh, ax=ax, label="output SNR at f0 (dB)")

    fig.tight_layout()
    plotting.save(fig, "e4_robustness")

    def line(label, r):
        k, g = interior_peak(r["snr"]), int(np.argmax(r["snr"]))
        edge = "" if g == k else f"  (global max at the D = {r['D'][g]:.4f} edge, "
        edge += "" if g == k else f"{r['snr'][g]:.2f} dB, intra-well not SR)"
        print(f"    {label}: SR peak D = {r['D'][k]:.4f}, {r['snr'][k]:.2f} dB{edge}")

    print("\n  amplitude family (f0 = %s):" % f0_base)
    for A, r in zip(amps, a_runs):
        line(f"A = {A}", r)
    print("  frequency family (A = %s):" % cfg["signal"]["A"])
    for f, r in zip(freqs, f_runs):
        line(f"f0 = {f}", r)
    print("  heatmap SR peak D per f0:")
    for f, d in zip(hm_f0, peak_D):
        print(f"    f0 = {f:.4f}: peak D = {d:.4f}")


if __name__ == "__main__":
    main()
