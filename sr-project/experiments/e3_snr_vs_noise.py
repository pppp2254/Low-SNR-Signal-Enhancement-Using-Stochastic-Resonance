"""E3: does the bistable system peak where theory says (PLAN.md 4)?

The headline result. SNR at f0 against noise intensity D, with bootstrap error
bars, the two-state theory curve, and a line at D = dU/2. Then the same sweep
with dt halved and with Heun, as the numerical-accuracy check of PLAN.md 2.3.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting
from srlab.bistable import simulate
from srlab.sweep import cached, chunked_snr, d_grid, load_config
from srlab.theory import barrier_height, critical_amplitude, two_state_snr_db

VARIANTS = {
    "e3": dict(label="baseline, Euler-Maruyama", dt=None, method="em", n_trials=None),
    "e3_dt_half": dict(label="dt halved", dt=None, method="em", n_trials=20),
    "e3_heun": dict(label="Heun", dt=None, method="heun", n_trials=20),
}


def theory_db(cfg, D, n_periods):
    return two_state_snr_db(cfg["system"]["a"], cfg["system"]["b"],
                            cfg["signal"]["A"], D, cfg["signal"]["f0"], n_periods)


def run(cfg, name, dt, method, n_trials):
    s, m = cfg["solver"], cfg["metrics"]
    f0, D = cfg["signal"]["f0"], d_grid(cfg)
    fs = 1.0 / (dt * s["decimate"])
    n_samples = int(round(s["n_periods"] / f0 * fs))
    seed = cfg["seeds"][name]

    params = dict(a=cfg["system"]["a"], b=cfg["system"]["b"], A=cfg["signal"]["A"],
                  f0=f0, D=D, dt=dt, n_periods=s["n_periods"],
                  discard_periods=s["discard_periods"], decimate=s["decimate"],
                  n_trials=n_trials, seed=seed, method=method,
                  n_side=m["n_side"], guard=m["guard"], n_boot=m["n_boot"])

    def compute():
        def make(lo, hi):
            # Chunked so the raw ensemble never exceeds the memory budget; a
            # spawned seed keeps each chunk independent and reproducible.
            _, x = simulate(a=params["a"], b=params["b"], A=params["A"], f0=f0,
                            D=D[lo:hi], dt=dt, n_periods=s["n_periods"],
                            discard_periods=s["discard_periods"],
                            decimate=s["decimate"], n_trials=n_trials,
                            seed=seed + lo, method=method)
            return x

        snr, band_lo, band_hi = chunked_snr(
            make, D.size, n_trials * n_samples * 4, fs, f0,
            n_side=m["n_side"], guard=m["guard"], n_boot=m["n_boot"], seed=seed)
        return {"D": D, "snr": snr, "lo": band_lo, "hi": band_hi}

    print(f"  {name}: dt={dt}, method={method}, {n_trials} trials")
    return cached(name, params, compute)


def main():
    cfg = load_config()
    a, b, A = cfg["system"]["a"], cfg["system"]["b"], cfg["signal"]["A"]
    s = cfg["solver"]
    barrier = barrier_height(a, b)
    d_opt = barrier / 2

    base = run(cfg, "e3", s["dt"], "em", cfg["ensemble"]["n_trials_final"])
    half = run(cfg, "e3_dt_half", cfg["sweeps"]["dt_check"], "em", 20)
    heun = run(cfg, "e3_heun", s["dt"], "heun", 20)

    D, snr, lo, hi = base["D"], base["snr"], base["lo"], base["hi"]
    k = int(np.argmax(snr))
    th = theory_db(cfg, D, s["n_periods"])
    k_th = int(np.argmax(th))

    plotting.use_style()
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.fill_between(D, lo, hi, color=plotting.COLORS["sim"], alpha=0.2,
                    label="bootstrap 16-84%")
    ax.plot(D, snr, "o-", ms=3.5, color=plotting.COLORS["sim"],
            label=f"simulation, {cfg['ensemble']['n_trials_final']} trials")
    ax.plot(D, th, color=plotting.COLORS["theory"], lw=1.2,
            label="two-state theory (PLAN.md 2.1)")
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0,
               label=f"theory peak, D = dU/2 = {d_opt}")
    ax.plot(D[k], snr[k], "*", ms=13, color=plotting.COLORS["optimal"], zorder=5,
            label=f"measured peak, D = {D[k]:.3f}, {snr[k]:.1f} dB")
    ax.set_xscale("log")
    # Theory runs to -50 dB at low D, where the two-state model has no intra-well
    # motion to describe. Clip the view to the measured range so the peak is legible.
    ax.set_ylim(snr.min() - 5, max(snr.max(), th[k_th]) + 3)
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("output SNR at f0 (dB)")
    ax.set_title(f"E3: SNR against noise, a = b = 1, A = {A} "
                 f"(A_c = {critical_amplitude(a, b):.3f}), f0 = {cfg['signal']['f0']}",
                 loc="left")
    ax.legend(loc="lower center", ncol=2)
    fig.tight_layout()
    plotting.save(fig, "e3_snr_vs_noise")

    fig, ax = plt.subplots(figsize=(7, 4.2))
    for r, label, color in (
        (base, f"baseline, dt = {s['dt']}, EM, {cfg['ensemble']['n_trials_final']} trials",
         plotting.COLORS["sim"]),
        (half, f"dt = {cfg['sweeps']['dt_check']}, EM, 20 trials", plotting.COLORS["alt1"]),
        (heun, f"dt = {s['dt']}, Heun, 20 trials", plotting.COLORS["alt2"]),
    ):
        ax.fill_between(r["D"], r["lo"], r["hi"], color=color, alpha=0.18)
        ax.plot(r["D"], r["snr"], "o-", ms=3, color=color, label=label)
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("output SNR at f0 (dB)")
    ax.set_title("E3 numerical checks: halved step and Heun against the baseline",
                 loc="left")
    ax.legend(loc="lower center")
    fig.tight_layout()
    plotting.save(fig, "e3_numerical_checks")

    overlap = [
        float(max(base["lo"][i], r["lo"][i]) <= min(base["hi"][i], r["hi"][i]))
        for r in (half, heun) for i in range(D.size)
    ]
    n = D.size
    print(f"\n  measured peak: D = {D[k]:.4f}, SNR = {snr[k]:.2f} dB "
          f"[{lo[k]:.2f}, {hi[k]:.2f}]")
    print(f"  theory peak:   D = {d_opt} (dU/2), curve argmax D = {D[k_th]:.4f}")
    print(f"  peak error vs dU/2: {100 * abs(D[k] - d_opt) / d_opt:.1f}%")
    print(f"  theory level at its peak: {th[k_th]:.2f} dB, "
          f"offset from measurement: {snr[k] - th[k_th]:+.2f} dB")
    print(f"  dt-halved band overlaps baseline at {int(sum(overlap[:n]))}/{n} D values")
    print(f"  Heun band overlaps baseline at {int(sum(overlap[n:]))}/{n} D values")


if __name__ == "__main__":
    main()
