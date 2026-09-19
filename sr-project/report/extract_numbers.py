"""Extract every number the report cites, straight from results/.

Not named numbers.py: that shadows the stdlib numbers module, which numpy
imports internally, and breaks every import in this directory.

PLAN.md WP9 requires that the report pull each value from
results/*.npz or *.json rather than from memory. This script is the single
place that reads them, and it writes report/numbers.md. Any figure quoted in
report.md should appear here; if it does not, it is unverified and belongs
behind a TODO.
"""

import json

import numpy as np

from srlab.metrics import pd_at_pfa, roc_auc
from srlab.sweep import ROOT, interior_peak, load_config
from srlab.theory import barrier_height, critical_amplitude, two_state_snr_db

OUT = ROOT / "report" / "numbers.md"


def load(name):
    p = ROOT / "results" / f"{name}.npz"
    return dict(np.load(p, allow_pickle=False)) if p.exists() else None


def rows_theory(cfg):
    a, b, A = cfg["system"]["a"], cfg["system"]["b"], cfg["signal"]["A"]
    barrier = float(barrier_height(a, b))
    return [
        ("critical amplitude A_c", f"{float(critical_amplitude(a, b)):.4f}", "theory.py"),
        ("barrier dU", f"{barrier:.4f}", "theory.py"),
        ("predicted SNR optimum dU/2", f"{barrier / 2:.4f}", "theory.py"),
        ("baseline amplitude A", f"{A}", "configs/default.yaml"),
        ("A as a fraction of A_c", f"{A / float(critical_amplitude(a, b)):.3f}",
         "derived"),
    ]


def rows_e1(cfg):
    v = load("e1")
    if v is None:
        return []
    return [(f"E1 SNR at D = {d}", f"{s:.2f} dB", "results/e1.npz")
            for d, s in zip(cfg["experiments"]["e1"]["D"], v["snr"])]


def rows_e2():
    v = load("e2")
    if v is None:
        return []
    k = int(np.argmax(v["snr"]))
    n_floor = int((v["snr"] <= -29.99).sum())
    return [
        ("E2 peak sigma", f"{v['sigma'][k]:.3f}", "results/e2.npz"),
        ("E2 peak SNR", f"{v['snr'][k]:.2f} dB", "results/e2.npz"),
        ("E2 sigma values at the -30 dB floor", f"{n_floor}", "results/e2.npz"),
        ("E2 SNR at the largest sigma", f"{v['snr'][-1]:.2f} dB", "results/e2.npz"),
    ]


def rows_e3(cfg):
    v = load("e3")
    if v is None:
        return []
    D, snr = v["D"], v["snr"]
    k = int(np.argmax(snr))
    barrier = float(barrier_height(cfg["system"]["a"], cfg["system"]["b"]))
    th = two_state_snr_db(cfg["system"]["a"], cfg["system"]["b"], cfg["signal"]["A"],
                          D, cfg["signal"]["f0"], cfg["solver"]["n_periods"])
    kt = int(np.argmax(th))
    out = [
        ("E3 measured peak D", f"{D[k]:.4f}", "results/e3.npz"),
        ("E3 peak SNR", f"{snr[k]:.2f} dB", "results/e3.npz"),
        ("E3 peak bootstrap band", f"[{v['lo'][k]:.2f}, {v['hi'][k]:.2f}] dB",
         "results/e3.npz"),
        ("E3 peak error against dU/2",
         f"{100 * abs(D[k] - barrier / 2) / (barrier / 2):.1f} %", "derived"),
        ("E3 D grid step", f"{100 * (D[1] / D[0] - 1):.1f} %", "derived"),
        ("E3 theory level at its peak", f"{th[kt]:.2f} dB", "theory.py"),
        ("E3 theory offset from measurement", f"{snr[k] - th[kt]:+.2f} dB", "derived"),
        ("E3 SNR at the lowest D", f"{snr[0]:.2f} dB", "results/e3.npz"),
        ("E3 SNR at the highest D", f"{snr[-1]:.2f} dB", "results/e3.npz"),
    ]
    for name, label in (("e3_dt_half", "dt halved"), ("e3_heun", "Heun")):
        w = load(name)
        if w is not None:
            out += [
                (f"E3 {label} peak D", f"{w['D'][int(np.argmax(w['snr']))]:.4f}",
                 f"results/{name}.npz"),
                (f"E3 {label} maximum deviation",
                 f"{np.abs(w['snr'] - snr).max():.2f} dB", "derived"),
            ]
    return out


def rows_e4(cfg):
    out = []
    for A in cfg["sweeps"]["A"]:
        v = load(f"e4_A{A}")
        if v is None:
            continue
        k = interior_peak(v["snr"])
        g = int(np.argmax(v["snr"]))
        out.append((f"E4 SR peak D at A = {A}", f"{v['D'][k]:.4f}",
                    f"results/e4_A{A}.npz"))
        if g != k:
            out.append((f"E4 global max at A = {A} (intra-well, not SR)",
                        f"{v['snr'][g]:.2f} dB at D = {v['D'][g]:.4f}", "derived"))
    for f in cfg["sweeps"]["f0"]:
        v = load(f"e4_f{f}")
        if v is None:
            continue
        k = interior_peak(v["snr"])
        out.append((f"E4 SR peak D at f0 = {f}", f"{v['D'][k]:.4f}",
                    f"results/e4_f{f}.npz"))
    return out


def rows_e5(cfg):
    e5 = cfg["experiments"]["e5"]
    out = []
    for A, D_ins, tag in ((cfg["signal"]["A"], e5["D_in"], ""),
                          (e5["A_linear"], e5["D_in_linear"], "_lin")):
        for d in D_ins:
            v = load(f"e5{tag}_D{d}")
            if v is None:
                continue
            k = interior_peak(v["snr"]) or int(np.argmax(v["snr"]))
            out.append((f"E5 peak a at A = {A}, D_in = {d}",
                        f"{v['a'][k]:.3f} against predicted {4 * d:.2f}",
                        f"results/e5{tag}_D{d}.npz"))
    return out


def rows_e6(cfg):
    out = []
    pfa = cfg["metrics"]["pfa"]
    v = load("e6_det")
    if v is not None:
        D = v["D"]
        auc = np.array([roc_auc(v["h0"][i], v["h1"][i]) for i in range(D.size)])
        pd = np.array([pd_at_pfa(v["h0"][i], v["h1"][i], pfa) for i in range(D.size)])
        dip = int(np.argmin(auc[: D.size // 2]))
        rec = dip + int(np.argmax(auc[dip:]))
        dip_p = int(np.argmin(pd[: D.size // 2]))
        rec_p = dip_p + int(np.argmax(pd[dip_p:]))
        out += [
            ("E6 bistable AUC at the lowest D", f"{auc[0]:.3f}", "results/e6_det.npz"),
            ("E6 bistable AUC dip", f"{auc[dip]:.3f} at D = {D[dip]:.4f}", "derived"),
            ("E6 bistable AUC recovery",
             f"{auc[rec]:.3f} at D = {D[rec]:.4f}", "derived"),
            ("E6 bistable AUC at the highest D", f"{auc[-1]:.3f}",
             "results/e6_det.npz"),
            ("E6 bistable Pd dip", f"{pd[dip_p]:.3f} at D = {D[dip_p]:.4f}", "derived"),
            ("E6 bistable Pd recovery",
             f"{pd[rec_p]:.3f} at D = {D[rec_p]:.4f}", "derived"),
            ("E6 bistable Pd gain over the dip",
             f"{pd[rec_p] - pd[dip_p]:+.3f}", "derived"),
        ]
    w = load("e6")
    if w is not None:
        auc = np.array([roc_auc(w["h0"][i], w["h1"][i]) for i in range(w["D"].size)])
        out.append(("E6 bistable AUC range at the baseline A",
                    f"{auc.min():.3f} to {auc.max():.3f}", "results/e6.npz"))
    t = load("e6_threshold")
    if t is not None:
        sg = t["sigma"]
        auc = np.array([roc_auc(t["h0"][i], t["h1"][i]) for i in range(sg.size)])
        pd = np.array([pd_at_pfa(t["h0"][i], t["h1"][i], pfa) for i in range(sg.size)])
        k = int(np.argmax(auc))
        out += [
            ("E6 threshold AUC at the lowest sigma", f"{auc[0]:.3f}",
             "results/e6_threshold.npz"),
            ("E6 threshold peak AUC", f"{auc[k]:.3f} at sigma = {sg[k]:.3f}",
             "results/e6_threshold.npz"),
            ("E6 threshold AUC at the highest sigma",
             f"{auc[-1]:.3f} at sigma = {sg[-1]:.0f}", "results/e6_threshold.npz"),
            ("E6 threshold Pd at the lowest sigma", f"{pd[0]:.3f}", "derived"),
            ("E6 threshold peak Pd", f"{pd.max():.3f}", "derived"),
            ("E6 threshold Pd at the highest sigma", f"{pd[-1]:.3f}", "derived"),
        ]
    return out


def rows_e7():
    p = ROOT / "results" / "e7_baselines.csv"
    if not p.exists():
        return []
    import csv

    out = []
    with p.open() as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        out.append((f"E7 {r['method']}",
                    f"SNR {float(r['snr_db']):.2f} dB, AUC {float(r['auc']):.3f}, "
                    f"Pd {float(r['pd_at_pfa']):.3f}", "results/e7_baselines.csv"))
    snr = {r["method"]: float(r["snr_db"]) for r in rows}
    best = max(snr, key=snr.get)
    sr = next(k for k in snr if k.startswith("SR at"))
    none = next(k for k in snr if k.startswith("bistable"))
    out += [
        ("E7 best method by SNR", f"{best} at {snr[best]:.2f} dB", "derived"),
        ("E7 gain from adding noise", f"{snr[sr] - snr[none]:+.2f} dB", "derived"),
        ("E7 linear advantage over SR", f"{snr[best] - snr[sr]:+.2f} dB", "derived"),
    ]
    return out


def check_report(flat):
    """Every bolded value in report.md must appear in the extracted set.

    PLAN.md WP9 requires every number in the text to come from results/. Bold
    is the convention this report uses for measured values, so this catches a
    figure that was typed from memory or left stale after a rerun.
    """
    import re

    src = ROOT / "report" / "report.md"
    if not src.exists():
        return
    haystack = " | ".join(flat.values())
    missing = []
    for claim in re.findall(r"\*\*([^*]+)\*\*", src.read_text()):
        claim = claim.strip()
        if not re.search(r"\d", claim):
            continue
        # Compare on the numeric tokens, so wording may differ from the table.
        nums = re.findall(r"-?\d+\.?\d*", claim)
        if not all(n in haystack for n in nums):
            missing.append(claim)
    if missing:
        print(f"  [check] {len(missing)} bolded claim(s) not found in results:")
        for c in missing:
            print(f"      {c}")
    else:
        print("  [check] every bolded claim in report.md traces to results/")


def main():
    cfg = load_config()
    groups = [
        ("Model and theory", rows_theory(cfg)),
        ("E1 regimes", rows_e1(cfg)),
        ("E2 threshold detector", rows_e2()),
        ("E3 SNR against noise", rows_e3(cfg)),
        ("E4 robustness", rows_e4(cfg)),
        ("E5 parameter tuning", rows_e5(cfg)),
        ("E6 detection", rows_e6(cfg)),
        ("E7 baselines", rows_e7()),
    ]
    lines = ["# Numbers cited in the report", "",
             "Generated by `report/extract_numbers.py` from `results/`. Do not edit by hand.",
             "Every value quoted in `report.md` must appear here.", ""]
    flat = {}
    for title, rows in groups:
        if not rows:
            continue
        lines += [f"## {title}", "", "| Quantity | Value | Source |",
                  "| --- | --- | --- |"]
        for name, value, source in rows:
            lines.append(f"| {name} | {value} | `{source}` |")
            flat[name] = value
        lines.append("")
    OUT.write_text("\n".join(lines))
    (ROOT / "report" / "numbers.json").write_text(json.dumps(flat, indent=2))
    print(f"  [report] report/numbers.md, numbers.json ({len(flat)} quantities)")
    check_report(flat)


if __name__ == "__main__":
    main()
