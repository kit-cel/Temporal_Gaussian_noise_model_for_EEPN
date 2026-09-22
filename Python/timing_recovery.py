"""Reverse the local order-1 (linear) FDPE, ideal genie-aided
(port of MATLAB/timing_recovery.m, phase_based only)."""

import numpy as np

from private.windowed_linear_fit import windowed_linear_fit
from private.apply_block_fdpe_filter import apply_block_fdpe_filter


def timing_recovery(Rx_symbols, Rx_phi, cfg):
    phase_sequence = Rx_phi[::cfg.oversampling_factor]
    EEPN_power, fit = windowed_linear_fit(phase_sequence, cfg)
    a0 = fit['phase_estimate']
    a1 = fit['slope']
    half = fit['half_window']

    Rx_symbols_order0_compensated = np.asarray(Rx_symbols).reshape(-1) * np.exp(-1j * a0)
    coefficients = a1.reshape(1, -1)
    Rx_symbols_TR = apply_block_fdpe_filter(Rx_symbols_order0_compensated, coefficients, half, cfg)

    phase_fit = {'phase_estimate': a0, 'slope': a1, 'half_window': half}
    return Rx_symbols_TR, EEPN_power, phase_fit
