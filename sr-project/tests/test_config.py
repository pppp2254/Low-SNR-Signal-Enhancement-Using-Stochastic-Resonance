"""WP0 check: the package imports and configs/default.yaml matches PLAN.md 2.2.

These assertions guard the two failure modes that PLAN.md section 9 lists first:
a signal that is not sub-threshold, and f0 not landing on an FFT bin.
"""
import math
from pathlib import Path

import yaml

CONFIG = Path(__file__).resolve().parents[1] / "configs" / "default.yaml"


def load():
    with CONFIG.open() as f:
        return yaml.safe_load(f)


def test_import_srlab():
    import srlab

    assert srlab.__version__


def test_baseline_values_match_plan():
    c = load()
    assert (c["system"]["a"], c["system"]["b"]) == (1.0, 1.0)
    assert c["signal"]["A"] == 0.3
    assert c["signal"]["f0"] == 0.01
    assert c["solver"]["dt"] == 0.01
    assert c["solver"]["decimate"] == 10
    assert c["solver"]["n_periods"] == 128
    assert c["solver"]["discard_periods"] == 10
    assert c["ensemble"]["n_trials"] == 20
    assert (c["noise"]["D_min"], c["noise"]["D_max"], c["noise"]["n_D"]) == (0.01, 2.0, 40)


def test_signal_is_sub_threshold():
    c = load()
    a, b = c["system"]["a"], c["system"]["b"]
    critical = math.sqrt(4 * a**3 / (27 * b))
    assert c["signal"]["A"] < critical, "A must stay below A_c or there is no SR peak"
    assert c["experiments"]["e2"]["A"] < c["experiments"]["e2"]["theta"]


def test_f0_lands_on_an_fft_bin():
    c = load()
    fs = 1.0 / (c["solver"]["dt"] * c["solver"]["decimate"])
    assert fs == c["solver"]["fs"]
    n_samples = c["solver"]["n_periods"] * fs / c["signal"]["f0"]
    assert n_samples == int(n_samples)
    bin_index = c["signal"]["f0"] * n_samples / fs
    assert bin_index == c["solver"]["n_periods"]
