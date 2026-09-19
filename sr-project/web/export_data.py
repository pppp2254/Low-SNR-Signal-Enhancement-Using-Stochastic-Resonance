"""Export results/*.npz to web/data.js for the dashboard.

The dashboard shows the real matplotlib figures, so this only carries the
numbers the images do not: the readout strip, the peak summary and the E3
table. Writes a JSON literal onto window.__DATA__ so the page needs no fetch
and works both from the Elysia server and as a published artifact.
"""

import json

import numpy as np

from srlab.sweep import ROOT, interior_peak, load_config
from srlab.theory import barrier_height, critical_amplitude, two_state_snr_db

OUT = ROOT / "web" / "data.js"

FIGURES = [
    ("e1_regimes", "E1", "Three noise regimes in time, with spectra"),
    ("e2_threshold", "E2", "Threshold detector, SNR against noise"),
    ("e3_snr_vs_noise", "E3", "SNR against noise intensity, with two-state theory"),
    ("e3_numerical_checks", "E3", "Halved step and Heun against the baseline"),
    ("e4_robustness", "E4", "Dependence on amplitude and frequency, and the heatmap"),
    ("e5_parameter_tuned", "E5", "Tuning a at fixed input noise"),
]


def r(x, sig=5):
    a = np.asarray(x, dtype=float)
    return [float(f"{v:.{sig}g}") for v in np.atleast_1d(a).ravel()]


def load(name):
    p = ROOT / "results" / f"{name}.npz"
    return dict(np.load(p)) if p.exists() else None


def peak_of(v, key):
    """Interior peak if there is one, else the global maximum."""
    k = interior_peak(v["snr"])
    edge = k is None
    k = int(np.argmax(v["snr"])) if edge else int(k)
    return float(v[key][k]), float(v["snr"][k]), edge


def main():
    cfg = load_config()
    a, b = cfg["system"]["a"], cfg["system"]["b"]
    A, f0 = cfg["signal"]["A"], cfg["signal"]["f0"]
    s = cfg["solver"]
    barrier = float(barrier_height(a, b))

    data = {
        "meta": {
            "a": a, "b": b, "A": A, "f0": f0,
            "Ac": round(float(critical_amplitude(a, b)), 4),
            "barrier": barrier, "dOpt": barrier / 2,
            "dt": s["dt"], "fs": s["fs"], "nPeriods": s["n_periods"],
            "discard": s["discard_periods"],
            "nTrials": cfg["ensemble"]["n_trials"],
            "nTrialsFinal": cfg["ensemble"]["n_trials_final"],
        },
        "figures": [{"file": f"figures/{n}.png", "id": e, "title": t}
                    for n, e, t in FIGURES],
        "stats": [],
        "summary": [],
    }

    e3 = load("e3")
    if e3 is not None:
        D = np.array(r(e3["D"]))
        th = two_state_snr_db(a, b, A, D, f0, s["n_periods"])
        k, kt = int(np.argmax(e3["snr"])), int(np.argmax(th))
        data["e3"] = {"D": r(D), "snr": r(e3["snr"]), "lo": r(e3["lo"]),
                      "hi": r(e3["hi"]), "theory": r(th), "peak": k}
        data["stats"] = [
            ["measured peak D", f"{D[k]:.4f}"],
            ["peak SNR", f"{e3['snr'][k]:.2f}", "dB"],
            ["theory peak D", f"{barrier / 2:.4f}"],
            ["peak error", f"{100 * abs(D[k] - barrier / 2) / (barrier / 2):.1f}", "%"],
            ["D grid step", f"{100 * (D[1] / D[0] - 1):.1f}", "%"],
            ["theory offset", f"{e3['snr'][k] - th[kt]:+.2f}", "dB"],
        ]
        for name, label in (("e3_dt_half", "dt halved"), ("e3_heun", "Heun")):
            v = load(name)
            if v is not None:
                data["summary"].append(
                    ["E3", f"{label}, peak D",
                     f"{v['D'][int(np.argmax(v['snr']))]:.4f}",
                     f"max deviation {np.abs(v['snr'] - e3['snr']).max():.2f} dB"])

    e1 = load("e1")
    if e1 is not None:
        for d, snr in zip(cfg["experiments"]["e1"]["D"], e1["snr"]):
            data["summary"].insert(0, ["E1", f"D = {d}", f"{snr:.2f} dB", ""])

    e2 = load("e2")
    if e2 is not None:
        k = int(np.argmax(e2["snr"]))
        data["stats"] += [["E2 peak sigma", f"{e2['sigma'][k]:.3f}"],
                          ["E2 peak SNR", f"{e2['snr'][k]:.2f}", "dB"]]
        data["summary"].append(["E2", "peak sigma", f"{e2['sigma'][k]:.3f}",
                                f"{e2['snr'][k]:.2f} dB"])

    for A_val in cfg["sweeps"]["A"]:
        v = load(f"e4_A{A_val}")
        if v is not None:
            d, snr, edge = peak_of(v, "D")
            data["summary"].append(["E4", f"A = {A_val}, SR peak D", f"{d:.4f}",
                                    f"{snr:.2f} dB" + (" (no interior peak)" if edge else "")])
    for f_val in cfg["sweeps"]["f0"]:
        v = load(f"e4_f{f_val}")
        if v is not None:
            d, snr, edge = peak_of(v, "D")
            data["summary"].append(["E4", f"f0 = {f_val}, SR peak D", f"{d:.4f}",
                                    f"{snr:.2f} dB" + (" (no interior peak)" if edge else "")])

    e5 = cfg["experiments"]["e5"]
    for A_val, D_ins, tag in ((A, e5["D_in"], ""),
                              (e5["A_linear"], e5["D_in_linear"], "_lin")):
        for d_in in D_ins:
            v = load(f"e5{tag}_D{d_in}")
            if v is not None:
                a_pk, snr, edge = peak_of(v, "a")
                data["summary"].append(
                    ["E5", f"A = {A_val}, D_in = {d_in}, peak a", f"{a_pk:.3f}",
                     f"predicted 4 D_in = {4 * d_in:.2f}"])

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text("window.__DATA__ = " + json.dumps(data, separators=(",", ":")) + ";\n")
    print(f"  wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size / 1024:.0f} KB), "
          f"{len(data['summary'])} summary rows, {len(data['figures'])} figures")


if __name__ == "__main__":
    main()
