"""Realize a higher-order FDPE compensation filter (port of MATLAB/private/apply_block_fdpe_filter.m).

Internal implementation detail: not part of the public API.
"""

import numpy as np


def apply_block_fdpe_filter(Rx_symbols_derotated, coefficients, half, cfg):
    """coefficients has shape (order, N): row p (0-indexed) holds the
    window-offset^(p+1) coefficient, matching the MATLAB 1-indexed
    convention (row 1 there is offset^1)."""
    Rx_symbols_derotated = np.asarray(Rx_symbols_derotated).reshape(-1)
    N = len(Rx_symbols_derotated)
    coefficients = np.atleast_2d(coefficients)
    order = coefficients.shape[0]

    taps = cfg.compensation_filter_taps
    block_length = cfg.compensation_block_length
    overlap = int(round(cfg.compensation_block_overlap * block_length))
    hop = block_length - overlap
    guard = overlap // 2

    starts = np.arange(0, N - block_length + 1, hop)  # 0-indexed
    n = np.arange(-taps, taps + 1)
    scale = half / taps

    Rx_symbols_out = np.zeros(N, dtype=complex)
    for start_idx in starts:
        end_idx = start_idx + block_length  # exclusive
        center_idx = (start_idx + end_idx - 1) // 2  # 0-indexed center sample

        write_start = start_idx + guard
        write_end = write_start + hop  # exclusive
        if write_end > N:
            continue

        phase_ramp = np.zeros(len(n))
        for p in range(1, order + 1):
            phase_ramp = phase_ramp + coefficients[p - 1, center_idx] * (n * scale) ** p
        phase_ramp = np.flip(phase_ramp)

        comp_filter = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(np.exp(-1j * phase_ramp))))
        block_out = np.convolve(Rx_symbols_derotated[start_idx:end_idx], comp_filter, mode='same')
        Rx_symbols_out[write_start:write_end] = block_out[guard:guard + hop]

    return Rx_symbols_out
