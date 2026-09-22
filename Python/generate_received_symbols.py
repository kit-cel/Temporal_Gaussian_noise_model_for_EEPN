"""Simulate the Tx/channel/Rx-frontend chain (port of MATLAB/generate_received_symbols.m)."""

import numpy as np

from helper_functions import upsample, rrc_pulse, qammod


def generate_received_symbols(cfg):
    rrcFilter = rrc_pulse(cfg.RRC_roll_off, cfg.RRC_span, cfg.oversampling_factor)
    H_cd = np.exp(1j * np.pi * cfg.lam ** 2 / cfg.c_0 * cfg.D_cd * cfg.fiber_length * cfg.f ** 2)

    Tx_int = np.random.randint(0, cfg.mod_order, size=int(cfg.num_transmission_symbols))
    Tx_symbols = qammod(cfg.mod_order, Tx_int)
    Tx_symbols_up = upsample(Tx_symbols, cfg.oversampling_factor)
    Tx_signal = np.convolve(Tx_symbols_up, rrcFilter, mode='same')

    X = np.fft.fftshift(np.fft.fft(Tx_signal, norm='ortho'))
    dispersed_signal = np.fft.ifft(np.fft.ifftshift(H_cd * X), norm='ortho')
    noise = (np.sqrt(cfg.system_noise_power / 2) * np.random.randn(len(Tx_signal))
             + 1j * np.sqrt(cfg.system_noise_power / 2) * np.random.randn(len(Tx_signal)))
    Rx_signal = dispersed_signal + noise

    Rx_delta_phi = np.sqrt(cfg.sigma2_LO) * np.random.randn(
        int(cfg.num_transmission_symbols * cfg.oversampling_factor))
    Rx_phi = np.cumsum(Rx_delta_phi) + 2 * np.pi * np.random.rand()

    Rx_signal_with_phase_noise = Rx_signal * np.exp(1j * Rx_phi)
    Y = np.fft.fftshift(np.fft.fft(Rx_signal_with_phase_noise, norm='ortho'))
    Rx_signal_cdc = np.fft.ifft(np.fft.ifftshift(Y * np.conj(H_cd)), norm='ortho')
    Rx_symbols_up = np.convolve(Rx_signal_cdc, rrcFilter, mode='same')
    Rx_symbols = Rx_symbols_up[::cfg.oversampling_factor]

    return Rx_symbols, Tx_symbols, Rx_phi, Rx_symbols_up
