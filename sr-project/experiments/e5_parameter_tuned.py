"""E5: tuning the system when the noise is fixed (PLAN.md 4, 2.5).

Correction to PLAN.md 2.5. It gives the optimum as dU = 2 D_in, so a = 8 D_in.
That is the condition for a D sweep, where the peak comes from balancing the
1/D^2 prefactor against exp(-dU/D). Sweeping a at fixed D, the only a
dependence is SNR proportional to a exp(-a/(4 D_in)), whose stationary point is
a = 4 D_in, that is dU = D_in. Confirmed numerically against two_state_snr:
its maximum sits at exactly 0.20, 0.40, 0.80 and 1.20 for D_in = 0.05, 0.1,
0.2 and 0.3. The figure marks a = 4 D_in.

The proposal's real setting: the input is already noisy and its noise cannot be
chosen, so the system is tuned instead. With b = a the wells stay at +-1 and the
barrier is dU = a/4, so sweeping a moves the barrier. The two-state optimum
dU = 2 D_in then predicts a = 8 D_in.

The sweep runs through the existing solver without modification. With b = a the
drift is a(x - x^3), so passing a as a column vector and D as an array of the
same D_in repeated makes axis 0 an a axis instead of a D axis.

Why not input_override: routing the input noise through the drive term would
need it pre-divided by dt, since the solver multiplies the drive by dt while a
noise increment must scale as sqrt(dt). That works for Euler-Maruyama but
breaks under Heun, and costs a step-rate array per trial. Passing D = D_in
gives identical physics with the correct scaling; the only thing lost is a
noise realisation shared across a values, which averaging over trials covers.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting
from srlab.bistable import simulate
from srlab.sweep import cached, chunked_snr, interior_peak, load_config
from srlab.theory import critical_amplitude


def input_snr_db(A, D_in, fs):
    """Broadband input SNR: signal power over noise power in the sampled band.

    Signal power is A^2 / 2. The one-sided noise PSD is 2 D_in, integrated over
    0 to fs/2, giving D_in fs. PLAN.md's "input SNR near -15 dB" is this
    broadband measure, not the per-bin ratio snr_db reports; per bin the same
    input sits far above 0 dB because the signal occupies one bin.
    """
    return 10 * np.log10((A**2 / 2) / (D_in * fs))


def run(cfg, D_in, index, A, tag):
    s, m = cfg["solver"], cfg["metrics"]
    f0, fs = cfg["signal"]["f0"], s["fs"]
    lo_a, hi_a = cfg["sweeps"]["a"]
    n_a = cfg["experiments"]["e5"]["n_a"]
    a = np.logspace(np.log10(lo_a), np.log10(hi_a), n_a)
    n_trials = cfg["ensemble"]["n_trials"]
    seed = cfg["seeds"]["e5"] + 100 * index
    n_samples = int(round(s["n_periods"] / f0 * fs))

    params = dict(a=a, D_in=D_in, A=A, f0=f0, dt=s["dt"], n_periods=s["n_periods"],
                  discard_periods=s["discard_periods"], decimate=s["decimate"],
                  n_trials=n_trials, seed=seed, n_side=m["n_side"],
                  guard=m["guard"], n_boot=m["n_boot"])

    def compute():
        def make(i, j):
            col = a[i:j].reshape(-1, 1)
            _, x = simulate(a=col, b=col, A=A, f0=f0,
                            D=np.full(j - i, float(D_in)), dt=s["dt"],
                            n_periods=s["n_periods"],
                            discard_periods=s["discard_periods"],
                            decimate=s["decimate"], n_trials=n_trials, seed=seed + i)
            return x

        snr, lo, hi = chunked_snr(make, n_a, n_trials * n_samples * 4, fs, f0,
                                  n_side=m["n_side"], guard=m["guard"],
                                  n_boot=m["n_boot"], seed=seed)
        return {"a": a, "snr": snr, "lo": lo, "hi": hi}

    print(f"  e5{tag} A={A} D_in={D_in}: input SNR {input_snr_db(A, D_in, fs):.1f} dB "
          f"(broadband), A/D_in={A / D_in:.2f}, {n_a} values of a")
    return cached(f"e5{tag}_D{D_in}", params, compute)


def panel(ax, cfg, A, D_ins, runs, title):
    """One condition: SNR against a for each D_in, with a = 8 D_in marked."""
    a_sub = A / 0.385  # A_c = 0.385 a when b = a
    fs = cfg["solver"]["fs"]
    for d, r, c in zip(D_ins, runs, ("alt1", "sim", "alt2")):
        col = plotting.COLORS[c]
        k = interior_peak(r["snr"]) or int(np.argmax(r["snr"]))
        ax.fill_between(r["a"], r["lo"], r["hi"], color=col, alpha=0.18)
        ax.plot(r["a"], r["snr"], "o-", ms=3, color=col,
                label=f"D_in = {d} ({input_snr_db(A, d, fs):.0f} dB in), "
                      f"peak a = {r['a'][k]:.2f} vs {4 * d:.1f}")
        ax.axvline(4 * d, color=col, ls=":", lw=1.2)
        ax.plot(r["a"][k], r["snr"][k], "*", ms=12, color=col, zorder=5)

    lo_a = cfg["sweeps"]["a"][0]
    if a_sub > lo_a:
        ax.axvspan(lo_a, a_sub, color=plotting.COLORS["theory"], alpha=0.07)
        ax.text(lo_a * 1.05, ax.get_ylim()[0],
                f" A above threshold (a < {a_sub:.2f})", fontsize=8,
                color=plotting.COLORS["theory"], va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel("system parameter a, with b = a (dimensionless)")
    ax.set_ylabel("output SNR at f0 (dB)")
    ax.set_title(title, loc="left")
    ax.legend(loc="lower center", fontsize=7.5)


def report(label, A, D_ins, runs):
    print(f"\n  {label} (A = {A}):")
    for d, r in zip(D_ins, runs):
        g = int(np.argmax(r["snr"]))
        k = interior_peak(r["snr"])
        edge = "" if k is not None else " NO INTERIOR PEAK"
        k = k if k is not None else g
        edge += "" if g == k else f" (global max at a = {r['a'][g]:.2f}, edge rise)"
        print(f"    D_in = {d}: peak a = {r['a'][k]:.3f}, {r['snr'][k]:.2f} dB, "
              f"predicted 4 D_in = {4 * d:.2f}, ratio {r['a'][k] / (4 * d):.2f}{edge}")


def main():
    cfg = load_config()
    e5, A_plan = cfg["experiments"]["e5"], cfg["signal"]["A"]
    A_lin = e5["A_linear"]

    D_a = e5["D_in"]
    runs_a = [run(cfg, d, i, A_plan, "") for i, d in enumerate(D_a)]
    D_b = e5["D_in_linear"]
    runs_b = [run(cfg, d, i, A_lin, "_lin") for i, d in enumerate(D_b)]

    plotting.use_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    panel(axes[0], cfg, A_plan, D_a, runs_a,
          f"A = {A_plan} as PLAN.md specifies: optimum falls in the "
          f"super-threshold region")
    panel(axes[1], cfg, A_lin, D_b, runs_b,
          f"A = {A_lin}, where A/D_in <= 1: optimum is testable")
    fig.suptitle("E5: tuning a at fixed input noise, dotted lines at a = 4 D_in "
                 "(PLAN.md 2.5 says 8 D_in; see docstring)", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e5_parameter_tuned")

    report("condition A, PLAN amplitude", A_plan, D_a, runs_a)
    report("condition B, linear response", A_lin, D_b, runs_b)


if __name__ == "__main__":
    main()
