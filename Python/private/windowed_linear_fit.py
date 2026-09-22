"""Local order-1 least-squares fit over the CD-memory window (port of MATLAB/private/windowed_linear_fit.m).

Internal implementation detail: not part of the public API (mirrors the
visibility that MATLAB's private/ folder enforces, even though Python's
import system does not enforce it the same way).
"""

import numpy as np


def windowed_linear_fit(y, cfg):
    y = np.asarray(y).reshape(-1)
    N = len(y)

    L = cfg.CD_memory + 1
    half = int(np.ceil((L - 1) / 2))
    i = np.arange(1, N + 1)

    s = np.maximum(1, i - half)
    e = np.minimum(N, i + half)
    m = e - s + 1

    cs_y = np.concatenate(([0.0], np.cumsum(y)))
    cs_y2 = np.concatenate(([0.0], np.cumsum(y ** 2)))
    cs_i = np.concatenate(([0.0], np.cumsum(i.astype(float))))
    cs_i2 = np.concatenate(([0.0], np.cumsum((i.astype(float)) ** 2)))
    cs_iy = np.concatenate(([0.0], np.cumsum(i * y)))

    Sy = cs_y[e] - cs_y[s - 1]
    Syy = cs_y2[e] - cs_y2[s - 1]
    Si = cs_i[e] - cs_i[s - 1]
    Sii = cs_i2[e] - cs_i2[s - 1]
    Siy = cs_iy[e] - cs_iy[s - 1]

    D = m * Sii - Si ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        a0 = (Sii * Sy - Si * Siy) / D
        a1 = (m * Siy - Si * Sy) / D

    SSE = Syy - (a0 * Sy + a1 * Siy)
    with np.errstate(divide='ignore', invalid='ignore'):
        EEPN_power = SSE / m
    EEPN_power = np.where(D == 0, 0.0, EEPN_power)

    phase_fit = {
        'phase_symbols': y,
        'phase_estimate': a0 + a1 * i,
        'slope': a1,
        'half_window': half,
        'window_start': s,
        'window_end': e,
    }
    return EEPN_power, phase_fit
