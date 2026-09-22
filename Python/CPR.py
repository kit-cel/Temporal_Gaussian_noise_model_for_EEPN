"""Carrier phase recovery, ideal genie-aided order-0 FDPE compensation
(port of MATLAB/CPR.m, phase_based only)."""

import numpy as np

from helper_functions import movmean


def CPR(Rx_symbols, Rx_phi, cfg):
    Rx_symbols = np.asarray(Rx_symbols).reshape(-1)
    Rx_phi_symbols = Rx_phi[::cfg.oversampling_factor]
    phase_est = movmean(Rx_phi_symbols, cfg.CD_memory + 1)
    Rx_symbols_CPR = Rx_symbols * np.exp(-1j * phase_est)
    return Rx_symbols_CPR, phase_est
