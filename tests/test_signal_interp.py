import numpy as np
import pandas as pd

# Import through the package name used by the future pip distribution.
# `conftest.py` makes the current source tree discoverable when pytest is
# launched from the Documents workspace rather than from this directory.
from oswald.signal_interp import (
    AIC,
    Find_start,
    Find_start_burst,
    Get_filtered_signal,
    _aic_curve,
    _moving_mean,
)


def test_find_start_uses_midpoint_before_extrema():
    # The expected index is the midpoint projected backwards from the two
    # largest opposite-signed extrema, clipped at the beginning of the data.
    signal = np.array([0.0, 1.0, 0.0, -2.0, 0.0, 0.0])

    assert Find_start(signal) == 0


def test_find_start_burst_never_returns_negative_index():
    # A short record cannot contain all requested periods, but the algorithm
    # must still return a valid array index.
    signal = np.array([0.0, 2.0, 0.0, -1.0, 0.0])

    assert Find_start_burst(signal, 21) == 0


def test_aic_curve_matches_scalar_aic_values():
    # Compare the optimized implementation with the original public scalar
    # formula at every split point, including its NaN boundary values.
    signal = np.array([0.2, 1.0, 0.5, -0.4, 0.8, 0.1])
    curve = _aic_curve(signal)

    expected = np.array([AIC(k, signal.size, signal) for k in range(signal.size + 1)])

    np.testing.assert_allclose(curve, expected, equal_nan=True)


def test_moving_mean_matches_centered_windows():
    # The convolution helper returns one value for each complete window.
    signal = np.arange(7, dtype=float)

    np.testing.assert_allclose(
        _moving_mean(signal, 3), [np.mean(signal[i:i + 3]) for i in range(5)]
    )


def test_filtered_signal_returns_one_centered_signal_per_input():
    # Use a deterministic synthetic waveform so this test checks shape and
    # numerical validity without relying on any laboratory data files.
    times = [np.arange(200, dtype=float) * 1e-5]
    rough = [np.sin(2 * np.pi * 10_000 * times[0])]
    starts = np.array([10])
    isvp = pd.Series([True])
    input_freqs = pd.Series([10_000.0])

    filtered = Get_filtered_signal(
        times, starts, rough, isvp, input_freqs,
        freq_rangeP=np.array([100.0, 50_000.0]),
        freq_rangeS=np.array([100.0, 50_000.0]),
        min_frange=5_000.0,
    )

    assert len(filtered) == 1
    assert filtered[0].shape == rough[0].shape
    assert np.isfinite(filtered[0]).all()