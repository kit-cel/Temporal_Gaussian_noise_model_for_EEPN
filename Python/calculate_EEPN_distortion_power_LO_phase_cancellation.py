"""Model EEPN power after LO-PC
(port of MATLAB/calculate_EEPN_distortion_power_LO_phase_cancellation.m)."""

import numpy as np


def calculate_EEPN_distortion_power_LO_phase_cancellation(Rx_phi, cfg):
    x = Rx_phi[::cfg.oversampling_factor]
    x = np.asarray(x).reshape(-1)
    N = len(x)
    L = cfg.CD_memory + 1
    half = int(np.ceil((L - 1) / 2))

    cs1 = np.concatenate(([0.0], np.cumsum(x)))
    cs2 = np.concatenate(([0.0], np.cumsum(x ** 2)))

    idx = np.arange(1, N + 1)
    idx_start = np.maximum(1, idx - half)
    idx_end = np.minimum(N, idx + half)

    sum_x = cs1[idx_end] - cs1[idx_start - 1]
    sum_x2 = cs2[idx_end] - cs2[idx_start - 1]
    win_len = idx_end - idx_start + 1

    m1 = sum_x / win_len
    m2 = sum_x2 / win_len

    return m2 + x ** 2 - 2 * x * m1
