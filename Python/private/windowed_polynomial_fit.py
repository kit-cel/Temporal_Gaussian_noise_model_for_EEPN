"""Local order-m least-squares fit over the CD-memory window (port of MATLAB/private/windowed_polynomial_fit.m).

Internal implementation detail: not part of the public API.
"""

import numpy as np

from private.local_polynomial_fit_moments import local_polynomial_fit_moments


def windowed_polynomial_fit(y, order, cfg):
    y = np.asarray(y).reshape(-1)

    L = cfg.CD_memory + 1
    half = int(np.ceil((L - 1) / 2))

    y_padded = np.concatenate((np.full(half, y[0]), y, np.full(half, y[-1])))
    S, Syy = local_polynomial_fit_moments(y_padded, half, order)

    # float powers: half**powers as int would silently overflow int64 for
    # half in the thousands and order above ~5 (half**order can reach 1e36).
    powers = np.arange(order + 1, dtype=float).reshape(-1, 1)
    half_pow = float(half) ** powers
    S_norm = S / half_pow

    x = np.arange(-half, half + 1) / half
    V = np.vander(x, order + 1, increasing=True)
    G = V.T @ V

    beta_norm = np.linalg.solve(G, S_norm)
    SSE = Syy - np.sum(beta_norm * S_norm, axis=0)

    EEPN_power = np.real(SSE) / L

    phase_fit = {
        'phase_symbols': y,
        'coefficients': beta_norm / half_pow,
        'order': order,
        'half_window': half,
    }
    return EEPN_power, phase_fit
