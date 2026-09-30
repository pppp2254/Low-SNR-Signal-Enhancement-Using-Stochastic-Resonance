"""E11: does stochastic resonance apply to real sound?

Input is data/audio/input.wav, or a path given on the command line. With
neither, a synthetic three-note melody stands in so the script always runs.
Record your own with:

    arecord -f S16_LE -r 16000 -c 1 -d 5 data/audio/input.wav

The sound is made sub-threshold, then pushed through the two systems this
project studies. Each is scored by fidelity_db against the clean sound, and
every interesting case is written as a WAV so it can be listened to.

Route 1, threshold array. N comparators, each adding its own noise, outputs
averaged. Memoryless, so it runs at audio rate directly. With no noise the
sub-threshold sound never crosses and the output is silence; this is the
effect audio engineers call dither.

Route 2, the bistable model. It resonates near one frequency, so time is
stretched until the sound's dominant frequency lands on f0_model. Broadband
sound cannot all sit at one frequency, which is the limitation being tested.

A pure tone separates the two ways "works" can be meant. The same bistable
route on a 277 Hz tone gives SNR at 277 Hz a shallow maximum, 28.15 dB at
D = 0.15 against 23.78 dB at D = 0.8, while waveform fidelity stays below
0.30 dB at every D. The bump is small with only four trials per point; the
contrast between the two metrics is the result, not the bump. The hopping output is a random two-state signal:
the tone stands far above the background in its own FFT bin, but most of the
output's power is broadband switching. SR helps detect a tone; it does not
reproduce a sound.

The linear baseline is analytic: N linear sensors with the same noise,
averaged, give fidelity 10 log10(1 + N var(s) / sigma^2). It is what the
threshold array competes against, and noise can only hurt it.
"""

import hashlib
import sys
from math import gcd

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy.signal import resample_poly

from srlab import plotting
from srlab.bistable import simulate
from srlab.metrics import fidelity_db, snr_db
from srlab.sweep import ROOT, cached, load_config
from srlab.threshold import threshold_detector

OUT = ROOT / "results" / "audio"
DEFAULT_INPUT = ROOT / "data" / "audio" / "input.wav"


def load_sound(path, fs, max_seconds):
    """Mono, resampled to fs, trimmed, DC removed, peak 1."""
    rate, x = wavfile.read(path)
    if np.issubdtype(x.dtype, np.integer):
        x = x / np.iinfo(x.dtype).max
    x = np.asarray(x, dtype=float)
    if x.ndim == 2:
        x = x.mean(axis=1)
    g = gcd(int(rate), fs)
    x = resample_poly(x, fs // g, int(rate) // g)[: int(max_seconds * fs)]
    x = x - x.mean()
    peak = np.abs(x).max()
    if peak == 0:
        raise ValueError(f"{path} is silent")
    return x / peak


def synthetic_sound(fs, seconds=3.0):
    """Three notes with harmonics and a plucked envelope."""
    t = np.arange(int(seconds * fs)) / fs
    x = np.zeros_like(t)
    notes = (220.0, 277.18, 329.63)
    seg = t.size // len(notes)
    for i, f in enumerate(notes):
        tt = t[i * seg:(i + 1) * seg] - t[i * seg]
        env = np.minimum(1.0, tt / 0.02) * np.exp(-1.5 * tt)
        x[i * seg:(i + 1) * seg] = env * sum(np.sin(2 * np.pi * k * f * tt) / k
                                             for k in (1, 2, 3))
    x -= x.mean()
    return x / np.abs(x).max()


def dominant_frequency(x, fs):
    spec = np.abs(np.fft.rfft(x * np.hanning(x.size)))
    freq = np.fft.rfftfreq(x.size, 1 / fs)
    band = freq > 50
    return float(freq[band][np.argmax(spec[band])])


def to_wav(name, y, fs):
    """Peak-normalised to 0.9 for listening. An all-zero output stays silent."""
    y = np.asarray(y, dtype=float)
    y = y - y.mean()
    peak = np.abs(y).max()
    if peak:
        y = y / peak * 0.9
    OUT.mkdir(parents=True, exist_ok=True)
    wavfile.write(OUT / name, fs, (y * 32767).astype(np.int16))


def fitted(reference, output):
    """Output rescaled by least squares onto the reference, for overlaying."""
    o = output - output.mean()
    if not (o @ o):
        return np.zeros_like(reference)
    return o * (o @ (reference - reference.mean())) / (o @ o)


def threshold_route(cfg, s, fs):
    e = cfg["experiments"]["e11"]
    sigma = np.logspace(np.log10(e["sigma_min"]), np.log10(e["sigma_max"]), e["n_sigma"])
    seed = cfg["seeds"]["e11"]
    fid = {}
    for N in e["units"]:
        fid[N] = np.array([
            fidelity_db(s, threshold_detector(s, sg, e["theta"], N,
                                              seed + 1000 * N + i)[0].mean(axis=0))
            for i, sg in enumerate(sigma)])
        print(f"  threshold, {N:3d} units: best {fid[N].max():6.2f} dB "
              f"at sigma = {sigma[int(np.argmax(fid[N]))]:.3f}")
    n_max = max(e["units"])
    best = sigma[int(np.argmax(fid[n_max]))]
    outs = {
        "no_noise": threshold_detector(s, 0.0, e["theta"], n_max, seed)[0].mean(axis=0),
        "best_1": threshold_detector(s, sigma[int(np.argmax(fid[1]))], e["theta"], 1,
                                     seed + 1)[0].mean(axis=0),
        "best_max": threshold_detector(s, best, e["theta"], n_max, seed + 2)[0].mean(axis=0),
        "too_much": threshold_detector(s, sigma[-1], e["theta"], n_max,
                                       seed + 3)[0].mean(axis=0),
    }
    for name, y in outs.items():
        to_wav(f"threshold_{name}.wav", y, fs)
    linear = 10 * np.log10(1 + n_max * s.var() / sigma**2)
    return sigma, fid, linear, outs["best_max"]


def bistable_route(cfg, x, fs):
    e, sv = cfg["experiments"]["e11"], cfg["solver"]
    a, b, A, dt = cfg["system"]["a"], cfg["system"]["b"], cfg["signal"]["A"], sv["dt"]
    x = x[: int(e["bistable_seconds"] * fs)]
    f_dom = dominant_frequency(x, fs)
    # Stretch time so f_dom sits at f0_model: model time per audio sample.
    steps = max(1, round(f_dom / (e["f0_model"] * fs) / dt))
    n_total = x.size * steps
    drive = A * np.interp(np.arange(n_total + 1) / steps, np.arange(x.size), x)
    D = np.array([0.0] + list(e["bistable_D"]))
    # simulate() sizes the record as n_periods / f0 / dt. One nominal period of
    # exactly n_total steps makes the record the clip, and decimate = steps
    # brings the output back to audio rate.
    print(f"  bistable: f_dom = {f_dom:.1f} Hz -> f0_model {e['f0_model']}, "
          f"{steps} solver steps per audio sample, {n_total:,} steps per trial")
    _, y = simulate(a=a, b=b, A=0.0, f0=1.0 / (n_total * dt), D=D, dt=dt, n_periods=1,
                    discard_periods=0, decimate=steps, n_trials=e["bistable_trials"],
                    seed=cfg["seeds"]["e11"] + 7, input_override=drive)
    s = A * x
    fid = np.array([np.mean([fidelity_db(s, trial) for trial in y[i]])
                    for i in range(D.size)])
    k = 1 + int(np.argmax(fid[1:]))
    to_wav("bistable_no_noise.wav", y[0, 0], fs)
    to_wav("bistable_best.wav", y[k, 0], fs)
    to_wav("bistable_too_much.wav", y[-1, 0], fs)
    # The overlay uses the no-noise output: it is the bistable's best result on
    # sound, and the with-noise optimum reproduces almost nothing.
    return D, fid, f_dom, steps, y[0, 0].astype(float), s


def tone_contrast(cfg):
    """Bistable route on a pure tone, scored two ways.

    277 Hz over a whole number of seconds is a whole number of cycles, so the
    tone sits exactly on an FFT bin and snr_db can measure it.
    """
    e, sv = cfg["experiments"]["e11"], cfg["solver"]
    fs, dt, A, f_tone = e["fs_audio"], sv["dt"], cfg["signal"]["A"], 277.0
    x = np.sin(2 * np.pi * f_tone * np.arange(e["bistable_seconds"] * fs) / fs)
    steps = round(f_tone / (e["f0_model"] * fs) / dt)
    n_total = x.size * steps
    drive = A * np.interp(np.arange(n_total + 1) / steps, np.arange(x.size), x)
    D = np.array(e["bistable_D"])
    _, y = simulate(a=cfg["system"]["a"], b=cfg["system"]["b"], A=0.0,
                    f0=1.0 / (n_total * dt), D=D, dt=dt, n_periods=1, discard_periods=0,
                    decimate=steps, n_trials=e["bistable_trials"],
                    seed=cfg["seeds"]["e11"] + 11, input_override=drive)
    snr, _, _ = snr_db(y, fs, f_tone)
    fid = np.array([np.mean([fidelity_db(A * x, tr) for tr in y[i]])
                    for i in range(D.size)])
    print(f"  pure tone: SNR at {f_tone:.0f} Hz peaks at {snr.max():.2f} dB, "
          f"fidelity never above {fid.max():.2f} dB")
    return D, np.asarray(snr), fid


def main():
    cfg = load_config()
    e = cfg["experiments"]["e11"]
    fs = e["fs_audio"]
    src = sys.argv[1] if len(sys.argv) > 1 else (DEFAULT_INPUT if DEFAULT_INPUT.exists()
                                                  else None)
    if src:
        x = load_sound(src, fs, e["max_seconds"])
        label, digest = str(src), hashlib.sha1(open(src, "rb").read()).hexdigest()
    else:
        x = synthetic_sound(fs)
        label, digest = "synthetic three-note melody", "synthetic"
    print(f"  input: {label}, {x.size / fs:.2f} s at {fs} Hz")
    s = e["peak_fraction"] * e["theta"] * x
    to_wav("input_clean.wav", x, fs)

    def compute():
        sigma, fid, linear, thr_best = threshold_route(cfg, s, fs)
        D, bfid, f_dom, steps, bis_quiet, bis_ref = bistable_route(cfg, x, fs)
        tD, tsnr, tfid = tone_contrast(cfg)
        out = {"sigma": sigma, "linear": linear, "D": D, "bistable": bfid,
               "f_dom": np.array(f_dom), "steps": np.array(steps),
               "tone_D": tD, "tone_snr": tsnr, "tone_fid": tfid}
        out.update({f"thr_{N}": v for N, v in fid.items()})
        n = int(0.04 * fs)
        mid = x.size // 6
        out["snip_ref"] = s[mid:mid + n]
        out["snip_thr"] = fitted(s, thr_best)[mid:mid + n]
        out["snip_bis_ref"] = bis_ref[mid:mid + n]
        out["snip_bis"] = fitted(bis_ref, bis_quiet)[mid:mid + n]
        return out

    r = cached("e11", dict(e, input=digest, seed=cfg["seeds"]["e11"]), compute)
    units = e["units"]

    plotting.use_style()
    fig, grid = plt.subplots(2, 2, figsize=(13, 8.6))
    axes = [grid[0, 0], grid[0, 1], grid[1, 1], grid[1, 0]]
    ax = axes[0]
    for N, c in zip(units, ("alt1", "sim", "alt2")):
        v = r[f"thr_{N}"]
        k = int(np.argmax(v))
        ax.plot(r["sigma"], v, "o-", ms=3, color=plotting.COLORS[c],
                label=f"{N} threshold unit{'s' if N > 1 else ''}, best {v[k]:.1f} dB")
    ax.plot(r["sigma"], r["linear"], "--", color=plotting.COLORS["clean"],
            label=f"linear sensor, {max(units)} units")
    ax.set_xscale("log")
    ax.set_ylim(-1, max(r[f"thr_{max(units)}"].max(), 5) + 8)
    ax.set_xlabel("added noise sigma (dimensionless)")
    ax.set_ylabel("fidelity to the clean sound (dB)")
    ax.set_title("Threshold array: noise lets the sound through", loc="left")
    ax.legend(loc="upper right", fontsize=7.5)

    ax = axes[1]
    D, bf = r["D"], r["bistable"]
    ax.plot(D[1:], bf[1:], "o-", ms=3.5, color=plotting.COLORS["sim"],
            label="bistable, with noise")
    ax.axhline(bf[0], color=plotting.COLORS["theory"], ls="--", lw=1.1,
               label=f"no noise, intra-well only: {bf[0]:.1f} dB")
    ax.set_ylim(-1, bf[0] + 4)
    ax.set_xscale("log")
    ticks = [0.03, 0.1, 0.3, 0.8]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("fidelity to the clean sound (dB)")
    ax.set_title(f"Bistable, {float(r['f_dom']):.0f} Hz mapped to f0 = {e['f0_model']}",
                 loc="left")
    ax.legend(loc="center right", fontsize=7.5)

    ax = axes[2]
    t_ms = np.arange(r["snip_ref"].size) / fs * 1000
    ax.plot(t_ms, r["snip_ref"] / np.abs(r["snip_ref"]).max(),
            color=plotting.COLORS["clean"], lw=2.0, label="clean sound")
    ax.plot(t_ms, r["snip_thr"] / np.abs(r["snip_ref"]).max(),
            color=plotting.COLORS["sim"], lw=0.8, label="threshold array, best")
    ax.plot(t_ms, r["snip_bis"] / np.abs(r["snip_bis_ref"]).max(),
            color=plotting.COLORS["alt2"], lw=1.1, label="bistable, no noise")
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("amplitude, rescaled")
    ax.set_title("40 ms of waveform", loc="left")
    ax.legend(loc="lower right", fontsize=7.5)

    ax = axes[3]
    ax.plot(r["tone_D"], r["tone_snr"], "o-", ms=3.5, color=plotting.COLORS["sim"],
            label="SNR at 277 Hz")
    ax.plot(r["tone_D"], r["tone_fid"], "s--", ms=3.5, color=plotting.COLORS["theory"],
            label="waveform fidelity")
    ax.set_xscale("log")
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("noise intensity D (dimensionless)")
    ax.set_ylabel("dB")
    ax.set_title("Pure 277 Hz tone: resonance in SNR, not in the waveform", loc="left")
    ax.legend(loc="center right", fontsize=7.5)

    fig.suptitle(f"E11: stochastic resonance on sound, input = {label}", y=1.0)
    fig.tight_layout()
    plotting.save(fig, "e11_audio")

    print()
    for N in units:
        v = r[f"thr_{N}"]
        k = int(np.argmax(v))
        print(f"  threshold {N:3d} units: 0.00 dB with no noise -> "
              f"{v[k]:.2f} dB at sigma = {r['sigma'][k]:.3f} -> {v[-1]:.2f} dB at "
              f"sigma = {r['sigma'][-1]:.0f}")
    kb = 1 + int(np.argmax(bf[1:]))
    print(f"  bistable: no noise {bf[0]:.2f} dB, best with noise {bf[kb]:.2f} dB "
          f"at D = {D[kb]}, most noise {bf[-1]:.2f} dB")
    kt = int(np.argmax(r["tone_snr"]))
    print(f"  pure tone: SNR {r['tone_snr'][0]:.2f} -> {r['tone_snr'][kt]:.2f} dB at "
          f"D = {r['tone_D'][kt]} -> {r['tone_snr'][-1]:.2f} dB; fidelity "
          f"{r['tone_fid'].min():.2f} to {r['tone_fid'].max():.2f} dB")
    print(f"  listen: {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
