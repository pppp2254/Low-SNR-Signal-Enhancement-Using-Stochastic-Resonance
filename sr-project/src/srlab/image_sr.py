"""Low-light image extension, E10 (PLAN.md section 6).

Works on the V channel of HSV only; colour channels are left untouched.

Planned functions
-----------------
darken(img, scale, noise_sigma, rng)
    Build a dark test image from a well-lit one so a ground truth exists.
approach_a(v, sigma, theta, n_frames, seed)
    Threshold SR: add Gaussian noise, threshold each pixel, average n_frames
    binary frames. Reference [5]; reuses threshold.py.
approach_b(v, n_iter, a, b, dt)
    Dynamic SR: iterate the discretised bistable equation per pixel. After
    Chouhan, Jha and Biswas (2013); read the paper for a, b, the step size and
    the stopping rule, and cite the equation numbers here.
baselines(v)
    Histogram equalisation, CLAHE, gamma (0.4 to 0.5).
image_metrics(ref, out)
    PSNR, SSIM (with ground truth); entropy, mean brightness, RMS contrast
    (without).
"""
