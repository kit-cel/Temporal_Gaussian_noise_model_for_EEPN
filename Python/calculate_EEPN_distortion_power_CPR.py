"""Model EEPN power after ideal CPR (port of MATLAB/calculate_EEPN_distortion_power_CPR.m).

Equivalent to the "Minimal Model Implementation" one-liner in the README
(movvar(Rx_phi, CD_memory+1)), just computed via cumulative sums.
"""

import numpy as np


def calculate_EEPN_distortion_power_CPR(Rx_phi, cfg):
    y = Rx_phi[::cfg.oversampling_factor]
    y = np.asarray(y).reshape(-1)
    N = len(y)
    L = cfg.CD_memory + 1
    half = int(np.ceil((L - 1) / 2))

    cs1 = np.concatenate(([0.0], np.cumsum(y)))
    cs2 = np.concatenate(([0.0], np.cumsum(y ** 2)))

    idx = np.arange(1, N + 1)
    idx_start = np.maximum(1, idx - half)
    idx_end = np.minimum(N, idx + half)

    sum_y = cs1[idx_end] - cs1[idx_start - 1]
    sum_y2 = cs2[idx_end] - cs2[idx_start - 1]
    win_len = idx_end - idx_start + 1

    m1 = sum_y / win_len
    m2 = sum_y2 / win_len

    return m2 - m1 ** 2
