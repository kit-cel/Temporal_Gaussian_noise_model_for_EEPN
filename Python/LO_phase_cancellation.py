"""Reverse the instantaneous LO phase walk-off, genie-aided (port of MATLAB/LO_phase_cancellation.m)."""

import numpy as np


def LO_phase_cancellation(Rx_symbols, Rx_phi, cfg):
    Rx_phi_symbols = Rx_phi[::cfg.oversampling_factor]
    return np.asarray(Rx_symbols).reshape(-1) * np.exp(-1j * Rx_phi_symbols)
