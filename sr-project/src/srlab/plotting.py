"""Shared figure style (PLAN.md 3).

Same colours and fonts everywhere. Axis labels in English with units, or
"dimensionless" where there are none.
"""

import pathlib

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

__all__ = ["use_style", "save", "COLORS"]

COLORS = {
    "sim": "#1f4e79",
    "theory": "#c0392b",
    "clean": "#7f8c8d",
    "optimal": "#e67e22",
    "alt1": "#2e8b57",
    "alt2": "#8e44ad",
}


def use_style():
    """Apply the project style to matplotlib's global rcParams."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "legend.fontsize": 8,
            "lines.linewidth": 1.3,
        }
    )


def save(fig, name, out_dir=None):
    """Write figures/<name>.png at 300 dpi and figures/<name>.pdf."""
    from .sweep import ROOT

    out = pathlib.Path(out_dir) if out_dir else ROOT / "figures"
    out.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}")
    print(f"  [figure] figures/{name}.png, .pdf")
