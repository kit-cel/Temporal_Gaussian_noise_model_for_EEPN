"""Derive dependent parameters from user-set cfg fields (port of MATLAB/init_parameter.m)."""

import numpy as np


def init_parameter(cfg):
    cfg.lam = 1.55e-6                       # Center wavelength in m
    cfg.c_0 = 3e8                           # Vacuum speed of light in m/s

    cfg.t = np.arange(cfg.num_transmission_symbols) / cfg.symbol_rate

    cfg.df = cfg.symbol_rate / cfg.num_transmission_symbols
    cfg.f = (np.arange(cfg.num_transmission_symbols * cfg.oversampling_factor)
              - cfg.num_transmission_symbols * cfg.oversampling_factor / 2) * cfg.df

    # Calculate the temporal broadening a pulse due to chromatic dispersion
    cfg.D = cfg.D_cd * cfg.fiber_length
    cfg.Delta_lambda = (cfg.lam ** 2) / cfg.c_0 * cfg.symbol_rate
    cfg.Delta_T = cfg.D * cfg.Delta_lambda
    cfg.CD_memory = int(round(cfg.Delta_T * cfg.symbol_rate))

    # Variance of the Wiener process
    cfg.sigma2_LO = 2 * np.pi * cfg.linewidth / (cfg.oversampling_factor * cfg.symbol_rate)

    cfg.system_noise_power = 10 ** (-cfg.snr / 10)
    return cfg
