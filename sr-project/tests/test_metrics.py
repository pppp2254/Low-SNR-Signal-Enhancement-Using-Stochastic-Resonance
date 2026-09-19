"""WP2: SNR, correlation and detection metrics (PLAN.md section 2.4)."""
import numpy as np
import pytest

from srlab import metrics, signals

FS, F0, N_PERIODS = 10.0, 0.01, 32


def noisy_sine(A, sigma, n_trials, seed=0):
    rng = np.random.default_rng(seed)
    t = signals.time_axis(F0, FS, N_PERIODS)
    return signals.sine(A, F0, t) + sigma * rng.standard_normal((n_trials, t.size))


def analytic_snr_db(A, sigma, n_samples):
    """For a sine in white noise measured per FFT bin.

    Signal bin power is A**2 N / (2 fs); the one-sided noise PSD is
    2 sigma**2 / fs; their ratio is A**2 N / (4 sigma**2), independent of fs.
    """
    return 10 * np.log10(A**2 * n_samples / (4 * sigma**2))


@pytest.mark.parametrize("A,sigma", [(0.05, 1.0), (0.02, 1.0), (0.1, 2.0)])
def test_known_snr_is_recovered_within_1_db(A, sigma):
    x = noisy_sine(A, sigma, n_trials=50)
    snr, S, N = metrics.snr_db(x, FS, F0)
    assert snr == pytest.approx(analytic_snr_db(A, sigma, x.shape[-1]), abs=1.0)
    assert S > 0 and N > 0


def test_pure_noise_is_near_or_below_zero_db():
    for seed in range(5):
        rng = np.random.default_rng(seed)
        x = rng.standard_normal((50, signals.n_samples_for(F0, FS, N_PERIODS)))
        assert metrics.snr_db(x, FS, F0)[0] < 1.0


def test_f0_must_land_on_a_bin():
    x = np.random.default_rng(0).standard_normal((4, 1000))
    metrics.snr_db(x, fs=10.0, f0=0.3, n_side=5)  # bin 30 exactly, fine
    with pytest.raises(ValueError, match="not an.*integer"):
        metrics.snr_db(x, fs=10.0, f0=0.3001, n_side=5)


def test_background_band_must_fit():
    x = np.random.default_rng(0).standard_normal((4, 1000))
    with pytest.raises(ValueError, match="DC bin"):
        metrics.snr_db(x, fs=10.0, f0=0.1, n_side=25)  # bin 10, band would go negative


def test_constant_record_raises_instead_of_returning_nan():
    with pytest.raises(ValueError, match="background power is zero"):
        metrics.snr_db(np.zeros((4, 1000)), fs=10.0, f0=0.3, n_side=5)


def test_snr_preserves_leading_sweep_axis():
    x = np.stack([noisy_sine(A, 1.0, 20, seed=i) for i, A in enumerate([0.02, 0.05, 0.1])])
    snr, _, _ = metrics.snr_db(x, FS, F0)
    assert snr.shape == (3,)
    assert np.all(np.diff(snr) > 0), "SNR must rise with signal amplitude"


def test_bootstrap_band_brackets_the_estimate_and_narrows_with_trials():
    wide = metrics.bootstrap_snr(noisy_sine(0.05, 1.0, 10, seed=1), FS, F0, n_boot=200)
    narrow = metrics.bootstrap_snr(noisy_sine(0.05, 1.0, 80, seed=1), FS, F0, n_boot=200)
    for snr, lo, hi in (wide, narrow):
        assert lo <= snr <= hi
    assert (narrow[2] - narrow[1]) < (wide[2] - wide[1])


def test_xcorr_peak_recovers_a_known_shift():
    """On an aperiodic signal the lag is exact; the coefficient cannot be.

    Correlation here is linear, not circular, so a lag of L costs L samples of
    overlap while the normalisation keeps the full energy of both records: the
    ceiling is about (n - L) / n.
    """
    clean = np.random.default_rng(0).standard_normal(4000)
    r, lag = metrics.xcorr_peak(clean, np.roll(clean, 40))
    assert lag == 40
    assert 0.98 < r < 1.0
    r_self, lag_self = metrics.xcorr_peak(clean, clean)
    assert r_self == pytest.approx(1.0, abs=1e-6) and lag_self == 0


def test_xcorr_peak_lag_is_ambiguous_for_a_periodic_input():
    """Why PLAN.md reserves this metric for pulses: on a sine the lag is not
    identifiable, and the triangular envelope pulls the reported peak inward."""
    t = signals.time_axis(F0, FS, N_PERIODS)
    clean = signals.sine(1.0, F0, t)
    _, lag = metrics.xcorr_peak(clean, np.roll(clean, 40))
    assert lag != 40


def test_xcorr_peak_finds_a_pulse_train_in_noise():
    """The E8 use case: a sub-threshold pulse train recovered from a noisy record."""
    rng = np.random.default_rng(3)
    t = signals.time_axis(F0, FS, N_PERIODS)
    clean = signals.pulse_train(t, rate=0.005, width=40.0, amplitude=0.3, rng=rng)
    noisy = clean + 0.2 * rng.standard_normal(t.size)
    r, lag = metrics.xcorr_peak(clean, noisy)
    assert lag == 0
    assert r > 0.4
    # An unrelated record must score far lower, and at a meaningless lag.
    r_unrelated, _ = metrics.xcorr_peak(clean, rng.standard_normal(t.size))
    assert r > 10 * abs(r_unrelated)


def test_xcorr_peak_is_near_zero_for_unrelated_signals():
    rng = np.random.default_rng(0)
    t = signals.time_axis(F0, FS, N_PERIODS)
    r, _ = metrics.xcorr_peak(signals.sine(1.0, F0, t), rng.standard_normal(t.size))
    assert abs(r) < 0.1


def test_auc_is_half_for_identical_distributions():
    rng = np.random.default_rng(0)
    h0 = rng.standard_normal(2000)
    h1 = rng.standard_normal(2000)
    assert metrics.roc_auc(h0, h1) == pytest.approx(0.5, abs=0.03)


def test_auc_rises_with_separation():
    rng = np.random.default_rng(0)
    h0 = rng.standard_normal(2000)
    aucs = [metrics.roc_auc(h0, rng.standard_normal(2000) + d) for d in (0.0, 1.0, 4.0)]
    assert aucs[0] == pytest.approx(0.5, abs=0.03)
    assert 0.7 < aucs[1] < 0.85
    assert aucs[2] > 0.99


def test_auc_handles_ties():
    """The threshold detector emits 0/1, so mid-ranks matter."""
    assert metrics.roc_auc([0, 0, 1, 1], [0, 0, 1, 1]) == pytest.approx(0.5)
    assert 0.5 < metrics.roc_auc([0, 0, 0, 1], [0, 1, 1, 1]) < 1.0


def test_pd_at_pfa():
    rng = np.random.default_rng(0)
    h0 = rng.standard_normal(4000)
    assert metrics.pd_at_pfa(h0, rng.standard_normal(4000)) == pytest.approx(0.05, abs=0.03)
    assert metrics.pd_at_pfa(h0, rng.standard_normal(4000) + 6.0) > 0.99
    # The realised false-alarm rate must not exceed the requested one.
    assert np.mean(h0 > np.quantile(h0, 0.95, method="higher")) <= 0.05


def test_float32_ensemble_matches_float64():
    """The solver returns float32. Casting the whole ensemble to float64 up
    front would double it and defeat the chunking, so the estimator must accept
    float32 directly and give the same answer."""
    x = noisy_sine(0.05, 1.0, n_trials=20).astype(np.float32)
    assert x.dtype == np.float32
    a = metrics.snr_db(x, FS, F0)[0]
    b = metrics.snr_db(x.astype(np.float64), FS, F0)[0]
    assert a == pytest.approx(b, abs=1e-5)
