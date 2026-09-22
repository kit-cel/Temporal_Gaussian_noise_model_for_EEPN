"""Model EEPN power after order-m adaptive filtering
(port of MATLAB/calculate_EEPN_distortion_power_full_compensation.m)."""

from private.windowed_polynomial_fit import windowed_polynomial_fit


def calculate_EEPN_distortion_power_full_compensation(Rx_phi, cfg):
    y = Rx_phi[::cfg.oversampling_factor]
    return windowed_polynomial_fit(y, cfg.compensation_order, cfg)
