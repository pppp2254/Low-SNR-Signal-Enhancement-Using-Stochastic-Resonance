"""E8: does stochastic resonance work for aperiodic input (PLAN.md 4)?

A random sub-threshold pulse train instead of a sine. Power at one frequency is
meaningless here, so the metric is the peak normalised cross-correlation
between the clean input and the output, which is what PLAN.md 2.4 reserves for
this case.

The pulse train drives the solver through input_override, evaluated on the
solver's own step grid rather than the decimated output grid. One train is
shared across every D and every trial, so each point on the curve sees exactly
the same input and differs only in the noise added to it.

The curve has the same two regimes as E6. At very low noise the particle stays
in one well and the pulses displace it slightly, which still correlates with
the input, so r is highest at the lowest D. Once hopping begins it destroys
that and r dips. Past the dip there is a genuine resonance, and it sits at the
same optimum the periodic drive gives, which is the point of the experiment:
the effect is not specific to a sine.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.bistable import simulate
from srlab.metrics import xcorr_peak
from srlab.sweep import cached, interior_peak, load_config
from srlab.theory import barrier_height, critical_amplitude


def main():
    cfg = load_config()
    s, e8 = cfg["solver"], cfg["experiments"]["e8"]
    f0, fs, dt = cfg["signal"]["f0"], s["fs"], s["dt"]
    a, b = cfg["system"]["a"], cfg["system"]["b"]
    n_periods, discard = s["n_periods"], s["discard_periods"]
    n_trials, seed = e8["n_trials"], cfg["seeds"]["e8"]
    d_opt = barrier_height(a, b) / 2

    D = np.logspace(np.log10(cfg["noise"]["D_min"]),
                    np.log10(cfg["noise"]["D_max"]), e8["n_D"])
    n_discard = int(round(discard / f0 / dt))
    n_record = int(round(n_periods / f0 / dt))

    params = dict(a=a, b=b, D=D, dt=dt, n_periods=n_periods, discard_periods=discard,
                  decimate=s["decimate"], n_trials=n_trials, seed=seed,
                  amplitude=e8["amplitude"], rate=e8["pulse_rate"],
                  width=e8["pulse_width"])

    def compute():
        rng = np.random.default_rng(seed)
        t_steps = np.arange(n_discard + n_record + 1) * dt
        drive = signals.pulse_train(t_steps, e8["pulse_rate"], e8["pulse_width"],
                                    e8["amplitude"], rng)
        print(f"  pulse train: {drive.mean() / e8['amplitude']:.1%} duty over "
              f"{t_steps[-1]:.0f} time units")
        t, x = simulate(a=a, b=b, A=0.0, f0=f0, D=D, dt=dt, n_periods=n_periods,
                        discard_periods=discard, decimate=s["decimate"],
                        n_trials=n_trials, seed=seed + 1, input_override=drive)
        clean = drive[n_discard::s["decimate"]][: t.size]
        r = np.array([[xcorr_peak(clean, x[i, j])[0] for j in range(n_trials)]
                      for i in range(D.size)])
        return {"D": D, "t": t, "clean": clean.astype(np.float32),
                "r_mean": r.mean(axis=1), "r_std": r.std(axis=1),
                "traces": x[:, 0, :].astype(np.float32)}

    res = cached("e8", params, compute)
    D, t, clean = res["D"], res["t"], res["clean"]
    r_mean, r_std, traces = res["r_mean"], res["r_std"], res["traces"]
    k = interior_peak(r_mean)
    g = int(np.argmax(r_mean))
    dip = int(np.argmin(r_mean[: D.size // 2]))
    show = t <= 3000

    plotting.use_style()
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.4),
                             gridspec_kw={"height_ratios": [1.1, 1]})

    ax = axes[0]
    ax.fill_between(D, r_mean - r_std, r_mean + r_std,
                    color=plotting.COLORS["sim"], alpha=0.2, label="spread over trials")
    ax.plot(D, r_mean, "o-", ms=3.5, color=plotting.COLORS["sim"],
            label="peak cross-correlation")
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0,
               label=f"periodic optimum dU/2 = {d_opt}")
    ax.plot(D[k], r_mean[k], "*", ms=13, color=plotting.COLORS["optimal"], zorder=5,
            label=f"SR peak r = {r_mean[k]:.3f} at D = {D[k]:.3f}")
    ax.plot(D[g], r_mean[g], "o", ms=7, mfc="none", color=plotting.COLORS["theory"],
            zorder=5, label=f"intra-well max r = {r_mean[g]:.3f} at D = {D[g]:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("peak cross-correlation (dimensionless)")
    ax.set_title(f"E8: aperiodic input, pulse amplitude {e8['amplitude']} "
                 f"(A_c = {critical_amplitude(a, b):.3f})", loc="left")
    ax.legend(loc="lower center", ncol=2)

    ax = axes[1]
    ax.plot(t[show], clean[show] / e8["amplitude"], color=plotting.COLORS["clean"],
            lw=1.4, label="input pulses (scaled)")
    ax.plot(t[show], traces[k][show], color=plotting.COLORS["sim"], lw=0.9,
            label=f"output at the SR peak D = {D[k]:.3f}, r = {r_mean[k]:.3f}")
    ax.plot(t[show], traces[dip][show], color=plotting.COLORS["theory"], lw=0.7,
            alpha=0.75, label=f"output at the dip D = {D[dip]:.3f}, r = {r_mean[dip]:.3f}")
    ax.set_ylim(-2.4, 2.4)
    ax.set_xlabel("time (model units)")
    ax.set_ylabel("x (dimensionless)")
    ax.set_title("Before and after, first 3000 time units", loc="left")
    ax.legend(loc="lower right", ncol=3, fontsize=7.5)

    fig.tight_layout()
    plotting.save(fig, "e8_pulses")

    print(f"\n  intra-well maximum r = {r_mean[g]:.3f} at D = {D[g]:.4f} (sweep edge)")
    print(f"  dip                r = {r_mean[dip]:.3f} at D = {D[dip]:.4f}")
    print(f"  SR peak            r = {r_mean[k]:.3f} at D = {D[k]:.4f}")
    print(f"  periodic optimum dU/2 = {d_opt}, so the aperiodic peak is "
          f"{100 * abs(D[k] - d_opt) / d_opt:.1f} % away")
    print(f"  r at the highest D = {D[-1]:.4f}: {r_mean[-1]:.3f}")


if __name__ == "__main__":
    main()
