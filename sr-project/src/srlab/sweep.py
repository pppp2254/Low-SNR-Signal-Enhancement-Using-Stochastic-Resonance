"""Generic sweep runner with caching (PLAN.md section 3).

Planned functions
-----------------
run_sweep(fn, params, out_dir="results", name=..., overwrite=False)
    Runs fn over a parameter grid and saves raw arrays to results/<name>.npz
    plus the exact parameters to results/<name>.json, so a figure can be
    rebuilt from results/ without rerunning the simulation. Skips the run when
    a cached result with identical parameters exists.
load_config(path="configs/default.yaml")
"""
