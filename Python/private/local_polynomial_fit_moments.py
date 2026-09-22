"""Windowed moments for local polynomial regression (port of MATLAB/private/local_polynomial_fit_moments.m).

Internal implementation detail: not part of the public API.
"""

import numpy as np


def _fftfilt(h, y):
    """Causal FIR filtering via FFT-based convolution, equivalent to MATLAB's
    fftfilt(h, y): out[n] = sum_k h[k]*y[n-k], truncated to len(y). numpy-only
    (no scipy), matching helper_functions.py's convention."""
    h = np.asarray(h)
    y = np.asarray(y)
    n = len(y) + len(h) - 1
    nfft = 1 << (n - 1).bit_length()
    Y = np.fft.fft(y, nfft)
    H = np.fft.fft(h, nfft)
    out = np.fft.ifft(Y * H)
    return out[:len(y)]


def local_polynomial_fit_moments(y, half, order):
    y = np.asarray(y).reshape(-1)
    N = len(y)
    L = 2 * half + 1
    interior_idx = np.arange(L, N + 1)  # 1-indexed, matching the MATLAB convention below

    is_complex = np.iscomplexobj(y)
    j = np.arange(half, -half - 1, -1).astype(float)  # matches the fftfilt causal-shift convention
    S = np.zeros((order + 1, len(interior_idx)), dtype=y.dtype if is_complex else float)
    for p in range(order + 1):
        filtered = _fftfilt(j ** p, y)
        S[p, :] = filtered[interior_idx - 1] if is_complex else filtered[interior_idx - 1].real

    filtered_yy = _fftfilt(np.ones(L), y ** 2)
    Syy = filtered_yy[interior_idx - 1] if is_complex else filtered_yy[interior_idx - 1].real

    return S, Syy
