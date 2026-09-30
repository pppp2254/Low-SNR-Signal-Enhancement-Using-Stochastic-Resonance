"""One hundred solver iterations, for the E4 and E5 parameter families.

At 100 steps of dt = 0.01 the window is one model time unit, so each iteration
gets about 20 pixels and every one is visible. That resolution shows what the
sweeps actually change, by splitting the Euler-Maruyama update

    x[k+1] - x[k] = (a x - b x^3) dt + A cos(2 pi f0 t) dt + sqrt(2 D dt) xi

into its three contributions. E4 varies the drive amplitude A, so only the
middle term moves. E5 varies a with b = a, so only the first term moves.

Frequency, E4's other family, is deliberately absent. Its periods run from 20
to 200 time units, so over a single time unit the drive is effectively constant
for every f0 in the sweep and the panels would be indistinguishable. Frequency
is a property of the envelope, not of an iteration.
"""

import matplotlib.pyplot as plt
import numpy as np

from srlab import plotting
from srlab.bistable import simulate
from srlab.sweep import cached, load_config
from srlab.theory import barrier_height

N_STEPS = 100
START = 4000  # steps into the record, past the transient


def run_full_rate(cfg, a, b, A, D, seed):
    """A short undecimated record: every solver step is kept."""
    s = cfg["solver"]
    t, x = simulate(a=a, b=b, A=A, f0=cfg["signal"]["f0"], D=[D], dt=s["dt"],
                    n_periods=2, discard_periods=s["discard_periods"],
                    decimate=1, n_trials=1, seed=seed)
    return t, x[0, 0]


def decompose(cfg, t, x, a, b, A, D):
    """Split each step into well, drive and noise, exactly as the solver does."""
    dt, f0 = cfg["solver"]["dt"], cfg["signal"]["f0"]
    well = (a * x[:-1] - b * x[:-1] ** 3) * dt
    drive = A * np.cos(2 * np.pi * f0 * t[:-1]) * dt
    noise = np.diff(x) - well - drive
    return well, drive, noise


def panel(cfg, axes, row, label, t, x, parts, colour, expected):
    """Draws the window and returns its statistics, so the printed summary and
    the legend cannot disagree about which steps they describe."""
    sl = slice(START, START + N_STEPS)
    steps = np.arange(N_STEPS)
    well, drive, noise = (p[sl] for p in parts)

    ax = axes[row, 0]
    ax.plot(steps, x[sl], "-o", ms=3, lw=0.9, color=colour)
    ax.set_ylabel(f"x\n{label}")
    ax.set_title(f"{label}: {N_STEPS} iterations, {N_STEPS * cfg['solver']['dt']:.0f} "
                 f"time unit", loc="left", fontsize=9)

    ax = axes[row, 1]
    ax.axhline(0, color=plotting.COLORS["clean"], lw=0.8)
    ax.plot(steps, well, lw=1.0, color=plotting.COLORS["alt1"],
            label=f"well, rms {np.abs(well).mean():.4f}")
    ax.plot(steps, drive, lw=1.0, color=plotting.COLORS["theory"],
            label=f"drive, rms {np.abs(drive).mean():.4f}")
    ax.plot(steps, noise, lw=0.9, color=plotting.COLORS["sim"], alpha=0.85,
            label=f"noise, rms {np.abs(noise).mean():.4f}")
    ax.set_ylabel("step contribution")
    ax.set_title(f"expected noise scale sqrt(2 D dt) = {expected:.4f}",
                 loc="left", fontsize=9)
    ax.legend(loc="upper right", fontsize=7)
    return [float(np.abs(p).mean()) for p in (well, drive, noise)]


def figure(cfg, cases, title, name):
    plotting.use_style()
    fig, axes = plt.subplots(len(cases), 2, figsize=(12.5, 2.6 * len(cases)),
                             squeeze=False)
    for row, (label, a, b, A, D, colour, seed) in enumerate(cases):
        t, x = run_full_rate(cfg, a, b, A, D, seed)
        parts = decompose(cfg, t, x, a, b, A, D)
        expected = float(np.sqrt(2 * D * cfg["solver"]["dt"]))
        well, drive, noise = panel(cfg, axes, row, label, t, x, parts,
                                   plotting.COLORS[colour], expected)
        print(f"  {label:26s} well {well:.5f}  drive {drive:.5f}  "
              f"noise {noise:.5f}  (expected noise {expected:.5f}, "
              f"noise/well {noise / well:.0f}x)")
    for ax in axes[-1]:
        ax.set_xlabel("solver iteration")
    fig.suptitle(title, y=1.0)
    fig.tight_layout()
    plotting.save(fig, name)


def main():
    cfg = load_config()
    a0, b0 = cfg["system"]["a"], cfg["system"]["b"]
    d_opt = barrier_height(a0, b0) / 2
    seed = cfg["seeds"]["e1"]

    print(f"E4 amplitudes at D = {d_opt}: only the drive term changes")
    figure(cfg,
           [(f"A = {A}", a0, b0, A, d_opt, c, seed + i)
            for i, (A, c) in enumerate(zip(cfg["sweeps"]["A"],
                                           ("alt1", "sim", "alt2")))],
           f"E4 at iteration scale: three drive amplitudes, a = b = {a0}, "
           f"D = {d_opt}", "e4_hundred_steps")

    e5 = cfg["experiments"]["e5"]
    A_lin, D_in = e5["A_linear"], e5["D_in_linear"][1]
    a_values = [cfg["sweeps"]["a"][0], 1.0, cfg["sweeps"]["a"][1]]
    print(f"\nE5 system parameters at A = {A_lin}, D_in = {D_in}: "
          f"only the well term changes")
    figure(cfg,
           [(f"a = b = {a}", a, a, A_lin, D_in, c, seed + 10 + i)
            for i, (a, c) in enumerate(zip(a_values, ("alt1", "sim", "alt2")))],
           f"E5 at iteration scale: three system parameters, A = {A_lin}, "
           f"D_in = {D_in}", "e5_hundred_steps")


if __name__ == "__main__":
    main()
