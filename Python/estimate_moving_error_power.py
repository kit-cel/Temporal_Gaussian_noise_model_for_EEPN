"""Estimate error power using moving variance (port of MATLAB/estimate_moving_error_power.m)."""

from helper_functions import movvar


def estimate_moving_error_power(Rx_symbols, Tx_symbols, cfg):
    sigma = movvar(Rx_symbols - Tx_symbols, cfg.block_length_SNR_evaluation)
    D = cfg.discard_symbols_analysis
    return sigma[D:len(sigma) - D - 1]
