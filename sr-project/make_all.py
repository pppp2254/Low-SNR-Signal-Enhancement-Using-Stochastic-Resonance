"""Regenerate every result, figure and report number (PLAN.md WP10).

Runs each experiment in order, then the report number check and the dashboard
export. Scripts cache on their exact parameters, so a second run only redraws;
pass --fresh to delete results/*.json first and rerun every simulation.

Heavy on a small machine: E10 peaks near 800 MB. On a 4 GB box run it as
    systemd-run --user --scope -p MemoryMax=1200M -p MemorySwapMax=0 python make_all.py
"""

import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent

# e10_frames and e10_denoise read e10_image's results; extract_numbers and export_data read all.
STEPS = [
    "experiments/e1_regimes.py",
    "experiments/e1_every_step.py",
    "experiments/e2_threshold.py",
    "experiments/e3_snr_vs_noise.py",
    "experiments/e4_robustness.py",
    "experiments/e45_hundred_steps.py",
    "experiments/e5_parameter_tuned.py",
    "experiments/e6_detection.py",
    "experiments/e7_baselines.py",
    "experiments/e8_pulses.py",
    "experiments/e9_residence_times.py",
    "experiments/e10_grating.py",
    "experiments/e10_image.py",
    "experiments/e10_frames.py",
    "experiments/e10_denoise.py",
    "experiments/e11_audio.py",
    "report/extract_numbers.py",
    "web/export_data.py",
]


def main():
    if "--fresh" in sys.argv:
        for p in (ROOT / "results").glob("*.json"):
            p.unlink()
    t_all = time.time()
    for step in STEPS:
        t = time.time()
        print(f"\n[make_all] {step}", flush=True)
        subprocess.run([sys.executable, str(ROOT / step)], cwd=ROOT, check=True)
        print(f"[make_all] {step} done in {time.time() - t:.0f} s", flush=True)
    print(f"\nall done in {(time.time() - t_all) / 60:.1f} min")


if __name__ == "__main__":
    main()
