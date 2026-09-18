"""Export results/*.npz to web/data.js for the dashboard.

Writes a JSON literal onto window.__DATA__ so the page needs no fetch and works
both from the Elysia server and as a published artifact.
"""

import json
import pathlib

import numpy as np

from srlab.sweep import ROOT, load_config
from srlab.theory import barrier_height, critical_amplitude, two_state_snr_db

OUT = ROOT / "web" / "data.js"


def r(x, sig=5):
    a = np.asarray(x, dtype=float)
    return [float(f"{v:.{sig}g}") for v in np.atleast_1d(a).ravel()]


def load(name):
    p = ROOT / "results" / f"{name}.npz"
    return dict(np.load(p)) if p.exists() else None


def main():
    cfg = load_config()
    a, b = cfg["system"]["a"], cfg["system"]["b"]
    A, f0 = cfg["signal"]["A"], cfg["signal"]["f0"]
    s = cfg["solver"]
    data = {"meta": {
        "a": a, "b": b, "A": A, "f0": f0,
        "Ac": round(float(critical_amplitude(a, b)), 4),
        "barrier": float(barrier_height(a, b)),
        "dOpt": float(barrier_height(a, b)) / 2,
        "dt": s["dt"], "fs": s["fs"], "nPeriods": s["n_periods"],
        "discard": s["discard_periods"], "decimate": s["decimate"],
        "nTrials": cfg["ensemble"]["n_trials"],
        "nTrialsFinal": cfg["ensemble"]["n_trials_final"],
        "theta": cfg["experiments"]["e2"]["theta"], "e2A": cfg["experiments"]["e2"]["A"],
    }}

    e1 = load("e1")
    if e1 is not None:
        data["e1"] = {
            "D": cfg["experiments"]["e1"]["D"],
            "t": r(e1["t"], 6),
            "traces": [r(row, 4) for row in e1["trace"]],
            "freq": r(e1["freq"], 6),
            "spec": [r(row, 4) for row in e1["spec"]],
            "snr": r(e1["snr"], 4),
            "k0": int(e1["k0"]),
        }

    e2 = load("e2")
    if e2 is not None:
        data["e2"] = {k: r(e2[k], 5) for k in ("sigma", "snr", "lo", "hi")}
        data["e2"]["peak"] = int(np.argmax(e2["snr"]))

    for name in ("e3", "e3_dt_half", "e3_heun"):
        v = load(name)
        if v is None:
            continue
        data[name] = {k: r(v[k], 5) for k in ("D", "snr", "lo", "hi")}
        data[name]["peak"] = int(np.argmax(v["snr"]))

    if "e3" in data:
        D = np.array(data["e3"]["D"])
        th = two_state_snr_db(a, b, A, D, f0, s["n_periods"])
        data["e3"]["theory"] = r(th, 5)
        data["e3"]["theoryPeak"] = int(np.argmax(th))

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text("window.__DATA__ = " + json.dumps(data, separators=(",", ":")) + ";\n")
    print(f"  wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size / 1024:.0f} KB), "
          f"sets: {sorted(k for k in data if k != 'meta')}")


if __name__ == "__main__":
    main()
