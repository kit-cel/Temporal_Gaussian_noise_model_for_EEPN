function [EEPN_power, phase_fit] = calculate_EEPN_distortion_power_full_compensation(Rx_phi, cfg)
%CALCULATE_EEPN_DISTORTION_POWER_FULL_COMPENSATION Model EEPN power after
%order-m adaptive filtering.
%   Rx_phi is the real, unwrapped LO phase vector at the oversampled rate.
%   cfg.oversampling_factor determines the symbol-rate phase samples.
%   cfg.compensation_order sets the polynomial order m of the local LO
%   phase fit (Ntilde in the paper notation). m=0 reproduces
%   calculate_EEPN_distortion_power_CPR and m=1 reproduces
%   calculate_EEPN_distortion_power_timing_recovery, here evaluated
%   through the shared private/windowed_polynomial_fit.m engine instead
%   of their closed-form solutions. See that file for the fit itself
%   (FFT-based moments, ill-conditioning-avoiding normalized offsets,
%   edge-replicated boundary padding).
%
%   EEPN_power is a column vector with one value per symbol. Background
%   noise is excluded, and no samples are discarded from the output.
%
%   Optional phase_fit supplies the quantities used by the compensation:
%     phase_symbols - Symbol-rate LO phase (column vector).
%     coefficients  - (order+1)-by-N local polynomial coefficients, one
%                     column per symbol, in powers of the window offset
%                     k in [-half_window, half_window] (not the absolute
%                     symbol index).
%     order         - Polynomial order m.
%     half_window   - Half-width of the centered window, in symbols.

    y = Rx_phi(1:cfg.oversampling_factor:end);
    [EEPN_power, phase_fit] = windowed_polynomial_fit(y, cfg.compensation_order, cfg);
end
