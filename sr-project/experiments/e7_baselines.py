"""E7: is stochastic resonance actually worth using (PLAN.md 4)?

One low-SNR input, five detectors. The input is a weak sine plus its own noise
of intensity D_in, which the receiver cannot choose. Every method sees an input
built the same way from the same parameters.

  1. bistable, no added noise      the system as-is, driven by the input alone
  2. threshold detector, no added  the simplest nonlinear detector
  3. linear band-pass              the f0 bin of the raw input. For a sine of
                                   known frequency in white Gaussian noise this
                                   is the matched filter, so it is the honest
                                   benchmark rather than a strawman
  4. SR at the optimal D           noise added until the total reaches dU/2
  5. parameter-tuned SR            input noise left alone, a tuned instead

Two numbers per method, because they answer different questions. SNR is the
ensemble-averaged periodogram ratio over a long record, at the baseline
amplitude. AUC is single-record detection over a short one, which is what a
detector actually faces, and it runs at the weaker experiments.e6.A_det: at
A = 0.3 every method scores AUC 1.000 and the column carries no information.
The two columns are therefore at different amplitudes, which the table says.

PLAN.md 8 asks this to be reported plainly: added noise helps a nonlinear
threshold-type detector, and is not expected to beat an optimal linear filter
for a sine in white Gaussian noise.
"""

import csv

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.bistable import simulate
from srlab.metrics import bootstrap_snr, f0_power, pd_at_pfa, roc_auc, snr_db
from srlab.sweep import ROOT, cached, load_config
from srlab.theory import barrier_height, critical_amplitude


def input_sigma(D_in, fs):
    """Sample standard deviation of noise with intensity D_in at rate fs.

    One-sided PSD is 2 D_in over 0 to fs/2, so the variance is D_in fs.
    """
    return np.sqrt(D_in * fs)


def raw_input(cfg, D_in, n_periods, n_trials, seed, A=None):
    """The noisy input itself, the thing every method receives."""
    s = cfg["solver"]
    f0, fs = cfg["signal"]["f0"], s["fs"]
    A = cfg["signal"]["A"] if A is None else A
    t = signals.time_axis(f0, fs, n_periods)
    rng = np.random.default_rng(seed)
    noise = input_sigma(D_in, fs) * rng.standard_normal((n_trials, t.size))
    return (signals.sine(A, f0, t) + noise)[None].astype(np.float32)


def bistable_run(cfg, D, a, n_periods, n_trials, seed, A=None):
    s = cfg["solver"]
    _, x = simulate(a=a, b=a, A=cfg["signal"]["A"] if A is None else A,
                    f0=cfg["signal"]["f0"], D=[D], dt=s["dt"], n_periods=n_periods,
                    discard_periods=s["discard_periods"], decimate=s["decimate"],
                    n_trials=n_trials, seed=seed)
    return x


def main():
    cfg = load_config()
    s, e6, e7 = cfg["solver"], cfg["experiments"]["e6"], cfg["experiments"]["e7"]
    A, f0, fs = cfg["signal"]["A"], cfg["signal"]["f0"], s["fs"]
    D_in, a0 = e7["D_in"], cfg["system"]["a"]
    d_opt = barrier_height(a0, a0) / 2
    n_trials, n_long = cfg["ensemble"]["n_trials"], s["n_periods"]
    n_short, n_runs = e6["n_periods"], e6["n_runs"]
    pfa, seed = cfg["metrics"]["pfa"], cfg["seeds"]["e7"]
    m = cfg["metrics"]

    # E5 found SNR declining monotonically with a inside the sub-threshold
    # region, so the best tunable a is the smallest one that keeps A below A_c.
    a_tuned = float(A / 0.385)
    A_det = e6["A_det"]

    params = dict(D_in=D_in, A=A, A_det=A_det, f0=f0, a_tuned=a_tuned, d_opt=d_opt,
                  n_long=n_long, n_short=n_short, n_trials=n_trials,
                  n_runs=n_runs, seed=seed, pfa=pfa)

    def compute():
        out = {}
        sigma = input_sigma(D_in, fs)
        from srlab.threshold import threshold_detector

        def short_pair(make):
            """H1 and H0 statistics from single short records, at A_det."""
            return f0_power(make(A_det, seed + 11), fs, f0)[0], \
                   f0_power(make(0.0, seed + 22), fs, f0)[0]

        # 1. bistable, no added noise
        out["m1_long"] = bistable_run(cfg, D_in, a0, n_long, n_trials, seed)
        out["m1_h1"], out["m1_h0"] = short_pair(
            lambda amp, sd: bistable_run(cfg, D_in, a0, n_short, n_runs, sd, A=amp))

        # 2. threshold detector, no added noise
        t_long = signals.time_axis(f0, fs, n_long)
        t_short = signals.time_axis(f0, fs, n_short)
        out["m2_long"] = threshold_detector(signals.sine(A, f0, t_long), sigma, 1.0,
                                            n_trials, seed + 1)
        out["m2_h1"] = f0_power(threshold_detector(signals.sine(A_det, f0, t_short),
                                                   sigma, 1.0, n_runs, seed + 33),
                                fs, f0)[0]
        out["m2_h0"] = f0_power(threshold_detector(signals.sine(0.0, f0, t_short), sigma,
                                                   1.0, n_runs, seed + 44), fs, f0)[0]

        # 3. linear band-pass on the raw input
        out["m3_long"] = raw_input(cfg, D_in, n_long, n_trials, seed + 2)
        out["m3_h1"] = f0_power(raw_input(cfg, D_in, n_short, n_runs, seed + 55,
                                          A=A_det), fs, f0)[0]
        rng = np.random.default_rng(seed + 66)
        noise_only = (input_sigma(D_in, fs)
                      * rng.standard_normal((n_runs, t_short.size)))[None].astype(np.float32)
        out["m3_h0"] = f0_power(noise_only, fs, f0)[0]

        # 4. SR at the optimal D
        out["m4_long"] = bistable_run(cfg, d_opt, a0, n_long, n_trials, seed + 3)
        out["m4_h1"], out["m4_h0"] = short_pair(
            lambda amp, sd: bistable_run(cfg, d_opt, a0, n_short, n_runs, sd, A=amp))

        # 5. parameter-tuned SR
        out["m5_long"] = bistable_run(cfg, D_in, a_tuned, n_long, n_trials, seed + 4)
        out["m5_h1"], out["m5_h0"] = short_pair(
            lambda amp, sd: bistable_run(cfg, D_in, a_tuned, n_short, n_runs, sd, A=amp))

        # Reduce here. Caching the raw ensembles made a 37 MB npz for numbers
        # that fit in a few hundred bytes.
        cols = {k: [] for k in ("snr", "lo", "hi", "auc", "pd")}
        for key in ("m1", "m2", "m3", "m4", "m5"):
            snr, lo, hi = bootstrap_snr(out[f"{key}_long"], fs, f0, n_side=m["n_side"],
                                        guard=m["guard"], n_boot=m["n_boot"], seed=seed)
            cols["snr"].append(float(np.ravel(snr)[0]))
            cols["lo"].append(float(np.ravel(lo)[0]))
            cols["hi"].append(float(np.ravel(hi)[0]))
            cols["auc"].append(roc_auc(out[f"{key}_h0"], out[f"{key}_h1"]))
            cols["pd"].append(pd_at_pfa(out[f"{key}_h0"], out[f"{key}_h1"], pfa))
        return {k: np.array(v) for k, v in cols.items()}

    r = cached("e7", params, compute)

    labels = [
        ("m1", f"bistable, no added noise (D = {D_in})"),
        ("m2", f"threshold detector, no added noise (theta = 1)"),
        ("m3", "linear band-pass on the raw input"),
        ("m4", f"SR at the optimal D = {d_opt}"),
        ("m5", f"parameter-tuned SR (a = b = {a_tuned:.2f})"),
    ]
    rows = [{"method": label, "snr_db": float(r["snr"][i]), "snr_lo": float(r["lo"][i]),
             "snr_hi": float(r["hi"][i]), "auc": float(r["auc"][i]),
             "pd_at_pfa": float(r["pd"][i])}
            for i, (_, label) in enumerate(labels)]

    out_csv = ROOT / "results" / "e7_baselines.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"  [csv] results/{out_csv.name}")

    short = ["bistable\nno added", "threshold\nno added", "linear\nband-pass",
             f"SR at\nD = {d_opt}", "parameter\ntuned SR"]
    colors = [plotting.COLORS[c] for c in ("alt1", "alt2", "theory", "sim", "optimal")]

    plotting.use_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axes[0]
    vals = [r_["snr_db"] for r_ in rows]
    err = np.array([[r_["snr_db"] - r_["snr_lo"] for r_ in rows],
                    [r_["snr_hi"] - r_["snr_db"] for r_ in rows]])
    ax.bar(short, vals, color=colors, yerr=err, capsize=3)
    ax.set_ylabel(f"output SNR at f0 (dB), {n_long} periods")
    ax.set_title(f"SNR per method, input noise D_in = {D_in}", loc="left")
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)

    ax = axes[1]
    aucs = [r_["auc"] for r_ in rows]
    ax.bar(short, aucs, color=colors)
    ax.axhline(0.5, color=plotting.COLORS["clean"], ls="--", lw=1.0, label="chance")
    ax.set_ylim(0.4, 1.02)
    ax.set_ylabel(f"ROC AUC (dimensionless), {n_short} periods")
    ax.set_title(f"Detection per method at A = {A_det}, {n_runs} runs per hypothesis",
                 loc="left")
    for i, v in enumerate(aucs):
        ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    ax.legend(loc="lower right")

    fig.suptitle(f"E7: five detectors on the same low-SNR input, A = {A} "
                 f"(A_c = {critical_amplitude(a0, a0):.3f})", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e7_baselines")

    best_snr = max(rows, key=lambda x: x["snr_db"])
    best_auc = max(rows, key=lambda x: x["auc"])
    print()
    for r_ in rows:
        print(f"  {r_['method']:48s} SNR {r_['snr_db']:7.2f} dB   "
              f"AUC {r_['auc']:.3f}   Pd@{pfa:.0%} {r_['pd_at_pfa']:.3f}")
    print(f"\n  SNR column at A = {A}, AUC column at A = {A_det}")
    print(f"  best SNR: {best_snr['method']} ({best_snr['snr_db']:.2f} dB)")
    print(f"  best AUC: {best_auc['method']} ({best_auc['auc']:.3f})")
    if "linear" in best_auc["method"]:
        print("  The linear band-pass wins on detection, as PLAN.md 8 expects.")


if __name__ == "__main__":
    main()
