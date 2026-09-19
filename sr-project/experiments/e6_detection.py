"""E6: does noise improve detection, not just SNR (PLAN.md 4)?

Hypothesis test on one short record. H0 is noise only, H1 is the weak signal
plus noise, 500 runs each. The statistic is the output periodogram power at f0
from a single 16-period record, which is short enough that detection is hard
and the AUC curve has room to show a maximum.

SNR and detectability are different claims. E3 shows noise maximising an
ensemble-averaged SNR; this asks whether a single record becomes easier to
classify, which is what a detector actually does.

Two systems, because the answer differs between them.

The bistable system has two regimes and they compete. Below about D = 0.02 the
particle stays in one well and responds linearly to the drive; that coherent
intra-well tone makes detection trivial, AUC 1.000. Once hopping starts it
destroys that tone and AUC collapses to a dip near D = 0.05. Past the dip the
curve recovers to a genuine stochastic-resonance maximum at dU/2, then falls
again. So noise does help detection, but only relative to the dip: the
intra-well regime at very low noise still beats the SR peak outright.

PLAN.md expects AUC near 0.5 at low D, which is what the two-state model
predicts. That is wrong here for the same reason it is wrong in E3, E4 and E5:
the two-state model has no intra-well motion.

The threshold detector has no intra-well response. Below threshold it emits
nothing at all, so at low noise both hypotheses give an empty record and AUC
sits at 0.5 exactly. That is where noise has something to contribute.

At the baseline amplitude A = 0.3 every AUC saturates at 1.000, so the
detection sweeps run at experiments.e6.A_det instead. The saturated curve is
kept and plotted, since it is the reason for the change.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting, signals
from srlab.bistable import simulate
from srlab.metrics import f0_power, pd_at_pfa, roc_auc
from srlab.threshold import threshold_detector
from srlab.sweep import ROOT, cached, load_config
from srlab.theory import barrier_height

MAX_CHUNK_BYTES = 3.0e8


def ensemble_stats(cfg, D, A, seed):
    """Per-record f0 power for every D, chunked to bound memory."""
    s, e6 = cfg["solver"], cfg["experiments"]["e6"]
    f0, fs = cfg["signal"]["f0"], s["fs"]
    n_runs, n_periods = e6["n_runs"], e6["n_periods"]
    n_samples = int(round(n_periods / f0 * fs))
    per_value = n_runs * n_samples * 4
    per_chunk = max(1, int(MAX_CHUNK_BYTES // per_value))

    out = np.empty((D.size, n_runs))
    for start in range(0, D.size, per_chunk):
        stop = min(start + per_chunk, D.size)
        _, x = simulate(a=cfg["system"]["a"], b=cfg["system"]["b"], A=A, f0=f0,
                        D=D[start:stop], dt=s["dt"], n_periods=n_periods,
                        discard_periods=s["discard_periods"],
                        decimate=s["decimate"], n_trials=n_runs, seed=seed + start)
        out[start:stop] = f0_power(x, fs, f0)
        del x
    return out


def threshold_stats(cfg, sigma, A, seed):
    """Per-record f0 power from the threshold detector, one sigma at a time."""
    s, e6 = cfg["solver"], cfg["experiments"]["e6"]
    f0, fs = cfg["signal"]["f0"], s["fs"]
    t = signals.time_axis(f0, fs, e6["n_periods"])
    sig = signals.sine(A, f0, t)
    out = np.empty((sigma.size, e6["n_runs"]))
    for i, sg in enumerate(sigma):
        y = threshold_detector(sig, sg, e6["thr_theta"], e6["n_runs"], seed + i)
        out[i] = f0_power(y, fs, f0)[0]
    return out


def curve(h0, h1, pfa):
    n = h0.shape[0]
    return (np.array([roc_auc(h0[i], h1[i]) for i in range(n)]),
            np.array([pd_at_pfa(h0[i], h1[i], pfa) for i in range(n)]))


def roc_points(h0, h1):
    thresholds = np.unique(np.concatenate([h0, h1]))[::-1]
    return ([(h0 > t).mean() for t in thresholds],
            [(h1 > t).mean() for t in thresholds])


def main():
    cfg = load_config()
    e6 = cfg["experiments"]["e6"]
    A, A_det, f0 = cfg["signal"]["A"], e6["A_det"], cfg["signal"]["f0"]
    fs = cfg["solver"]["fs"]
    d_opt = barrier_height(cfg["system"]["a"], cfg["system"]["b"]) / 2
    D = np.logspace(np.log10(cfg["noise"]["D_min"]),
                    np.log10(cfg["noise"]["D_max"]), e6["n_D"])
    sigma = np.logspace(np.log10(e6["thr_sigma_min"]),
                        np.log10(e6["thr_sigma_max"]), e6["n_sigma"])
    pfa, seed = cfg["metrics"]["pfa"], cfg["seeds"]["e6"]
    base = dict(D=D, f0=f0, n_runs=e6["n_runs"], n_periods=e6["n_periods"],
                dt=cfg["solver"]["dt"], decimate=cfg["solver"]["decimate"],
                discard_periods=cfg["solver"]["discard_periods"], pfa=pfa)

    def bistable_at(amp, tag, sd):
        def compute():
            print(f"  bistable A = {amp}: H0 then H1, {e6['n_runs']} runs "
                  f"at each of {D.size} D values")
            return {"D": D, "h0": ensemble_stats(cfg, D, 0.0, sd),
                    "h1": ensemble_stats(cfg, D, amp, sd + 5000)}
        return cached(tag, dict(base, A=amp, seed=sd), compute)

    # A = 0.3 is kept because its flat AUC is the finding, not a failed run.
    sat = bistable_at(A, "e6", seed)
    bis = bistable_at(A_det, "e6_det", seed + 100)

    def compute_thr():
        print(f"  threshold detector A = {e6['thr_A']}, theta = {e6['thr_theta']}, "
              f"{sigma.size} sigma values")
        return {"sigma": sigma,
                "h0": threshold_stats(cfg, sigma, 0.0, seed + 300),
                "h1": threshold_stats(cfg, sigma, e6["thr_A"], seed + 400)}

    thr = cached("e6_threshold", dict(base, sigma=sigma, A=e6["thr_A"],
                                      theta=e6["thr_theta"], seed=seed + 300),
                 compute_thr)

    auc_sat, _ = curve(sat["h0"], sat["h1"], pfa)
    auc_bis, pd_bis = curve(bis["h0"], bis["h1"], pfa)
    auc_thr, pd_thr = curve(thr["h0"], thr["h1"], pfa)
    kb, kt = int(np.argmax(auc_bis)), int(np.argmax(auc_thr))

    plotting.use_style()
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4))

    ax = axes[0, 0]
    for target, c in zip(e6["roc_at"], ("alt1", "sim", "alt2")):
        i = int(np.argmin(np.abs(D - target)))
        fpr, tpr = roc_points(bis["h0"][i], bis["h1"][i])
        ax.plot(fpr, tpr, color=plotting.COLORS[c],
                label=f"D = {D[i]:.3f}, AUC = {auc_bis[i]:.3f}")
    ax.plot([0, 1], [0, 1], color=plotting.COLORS["clean"], ls="--", lw=1.0,
            label="chance")
    ax.set_xlabel("false alarm rate (dimensionless)")
    ax.set_ylabel("detection rate (dimensionless)")
    ax.set_title(f"Bistable, A = {A_det}: ROC at three D", loc="left")
    ax.legend(loc="lower right")

    ax = axes[0, 1]
    ax.plot(D, auc_bis, "o-", ms=3, color=plotting.COLORS["sim"],
            label=f"AUC, A = {A_det}")
    ax.plot(D, pd_bis, "s-", ms=3, color=plotting.COLORS["alt1"],
            label=f"Pd at Pfa = {pfa:.0%}")
    ax.plot(D, auc_sat, "^--", ms=3, color=plotting.COLORS["theory"],
            label=f"AUC, A = {A} (saturated)")
    ax.axhline(0.5, color=plotting.COLORS["clean"], ls="--", lw=1.0)
    ax.axvline(d_opt, color=plotting.COLORS["optimal"], ls="--", lw=1.0,
               label=f"SNR optimum D = {d_opt}")
    ax.set_xscale("log")
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("AUC and detection rate (dimensionless)")
    ax.set_title("Bistable: SR peak at dU/2, but low noise still wins", loc="left")
    ax.legend(loc="lower left", fontsize=7.5)

    ax = axes[1, 0]
    for target, c in zip(e6["roc_sigma_at"], ("alt1", "sim", "alt2")):
        i = int(np.argmin(np.abs(sigma - target)))
        fpr, tpr = roc_points(thr["h0"][i], thr["h1"][i])
        ax.plot(fpr, tpr, color=plotting.COLORS[c],
                label=f"sigma = {sigma[i]:.3f}, AUC = {auc_thr[i]:.3f}")
    ax.plot([0, 1], [0, 1], color=plotting.COLORS["clean"], ls="--", lw=1.0,
            label="chance")
    ax.set_xlabel("false alarm rate (dimensionless)")
    ax.set_ylabel("detection rate (dimensionless)")
    ax.set_title(f"Threshold detector, A = {e6['thr_A']}: ROC at three sigma",
                 loc="left")
    ax.legend(loc="lower right")

    ax = axes[1, 1]
    ax.plot(sigma, auc_thr, "o-", ms=3, color=plotting.COLORS["sim"], label="AUC")
    ax.plot(sigma, pd_thr, "s-", ms=3, color=plotting.COLORS["alt1"],
            label=f"Pd at Pfa = {pfa:.0%}")
    ax.axhline(0.5, color=plotting.COLORS["clean"], ls="--", lw=1.0, label="chance")
    ax.plot(sigma[kt], auc_thr[kt], "*", ms=13, color=plotting.COLORS["optimal"],
            zorder=5, label=f"best AUC {auc_thr[kt]:.3f} at sigma = {sigma[kt]:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("noise standard deviation sigma (dimensionless)")
    ax.set_ylabel("AUC and detection rate (dimensionless)")
    ax.set_title("Threshold detector: noise does help detection", loc="left")
    ax.legend(loc="lower center", fontsize=7.5)

    fig.suptitle(f"E6: detection from one {e6['n_periods']}-period record, "
                 f"{e6['n_runs']} runs per hypothesis", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e6_detection")

    out = ROOT / "results" / "e6_detection.csv"
    with out.open("w") as fh:
        fh.write("system,noise,auc,pd_at_pfa\n")
        for d, a_, p_ in zip(D, auc_bis, pd_bis):
            fh.write(f"bistable_A{A_det},{d:.6g},{a_:.4f},{p_:.4f}\n")
        for d, a_ in zip(D, auc_sat):
            fh.write(f"bistable_A{A},{d:.6g},{a_:.4f},\n")
        for sg, a_, p_ in zip(sigma, auc_thr, pd_thr):
            fh.write(f"threshold,{sg:.6g},{a_:.4f},{p_:.4f}\n")
    print(f"  [csv] results/{out.name}")

    print(f"\n  bistable, A = {A} (baseline): AUC {auc_sat.min():.3f} to "
          f"{auc_sat.max():.3f} -- saturated, no information")
    # The global maximum sits in the intra-well regime. The stochastic
    # resonance is the recovery after hopping destroys that tone.
    dip = int(np.argmin(auc_bis[: D.size // 2]))
    rec = dip + int(np.argmax(auc_bis[dip:]))
    dip_p = int(np.argmin(pd_bis[: D.size // 2]))
    rec_p = dip_p + int(np.argmax(pd_bis[dip_p:]))
    print(f"  bistable, A = {A_det}: AUC {auc_bis[0]:.3f} at D = {D[0]:.4f} "
          f"(intra-well), dip {auc_bis[dip]:.3f} at D = {D[dip]:.4f}, "
          f"SR recovery {auc_bis[rec]:.3f} at D = {D[rec]:.4f} "
          f"(dU/2 = {d_opt}), then {auc_bis[-1]:.3f} at D = {D[-1]:.4f}")
    print(f"    Pd at Pfa = {pfa:.0%}: dip {pd_bis[dip_p]:.3f} at D = {D[dip_p]:.4f} "
          f"recovering to {pd_bis[rec_p]:.3f} at D = {D[rec_p]:.4f}, "
          f"a gain of {pd_bis[rec_p] - pd_bis[dip_p]:+.3f}")
    print(f"  threshold detector: AUC {auc_thr[0]:.3f} at sigma = {sigma[0]:.3f}, "
          f"peak {auc_thr[kt]:.3f} at sigma = {sigma[kt]:.3f}, "
          f"{auc_thr[-1]:.3f} at sigma = {sigma[-1]:.3f}")
    print(f"  Pd at Pfa = {pfa:.0%}: threshold detector {pd_thr[0]:.3f} -> "
          f"{pd_thr[kt]:.3f} -> {pd_thr[-1]:.3f}")


if __name__ == "__main__":
    main()
