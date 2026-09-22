"""Model EEPN power after CPR and timing recovery
(port of MATLAB/calculate_EEPN_distortion_power_timing_recovery.m)."""

from private.windowed_linear_fit import windowed_linear_fit


def calculate_EEPN_distortion_power_timing_recovery(Rx_phi, cfg):
    y = Rx_phi[::cfg.oversampling_factor]
    return windowed_linear_fit(y, cfg)
