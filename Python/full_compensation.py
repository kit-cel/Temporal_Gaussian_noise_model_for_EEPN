"""Reverse the local order-m FDPE, ideal genie-aided adaptive filtering
(port of MATLAB/full_compensation.m, phase_based only)."""

import numpy as np

from private.windowed_polynomial_fit import windowed_polynomial_fit
from private.apply_block_fdpe_filter import apply_block_fdpe_filter


def full_compensation(Rx_symbols, Rx_phi, cfg):
    order = cfg.compensation_order

    phase_sequence = Rx_phi[::cfg.oversampling_factor]
    EEPN_power, fit = windowed_polynomial_fit(phase_sequence, order, cfg)
    a0 = fit['coefficients'][0, :]
    coefficients = fit['coefficients'][1:, :]
    half = fit['half_window']

    Rx_symbols_order0_compensated = np.asarray(Rx_symbols).reshape(-1) * np.exp(-1j * a0)
    Rx_symbols_full = apply_block_fdpe_filter(Rx_symbols_order0_compensated, coefficients, half, cfg)

    phase_fit = {'coefficients': fit['coefficients'], 'order': order, 'half_window': half}
    return Rx_symbols_full, EEPN_power, phase_fit
