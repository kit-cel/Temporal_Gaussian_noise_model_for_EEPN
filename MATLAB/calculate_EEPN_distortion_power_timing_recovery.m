function [EEPN_power, phase_fit] = calculate_EEPN_distortion_power_timing_recovery(Rx_phi, cfg)
%CALCULATE_EEPN_DISTORTION_POWER_TIMING_RECOVERY Model EEPN power after CPR and timing recovery.
%   Rx_phi is the real, unwrapped LO phase vector at the oversampled rate.
%   cfg.oversampling_factor determines the symbol-rate phase samples.
%   Each output is the mean squared residual of a local least-squares line
%   fit to the LO phase (see private/windowed_linear_fit.m). In the
%   existing normalized EEPN model, this is the remaining distortion
%   power after ideal CPR and timing recovery.
%   The centered CD-memory window has 2*ceil(cfg.CD_memory/2)+1 samples
%   in the interior and shrinks at the sequence boundaries.
%   EEPN_power is a column vector with one value per symbol. Background
%   noise is excluded, and no samples are discarded from the output.
%
%   Optional phase_fit supplies the quantities used by the compensation:
%     phase_symbols  - Symbol-rate LO phase (column vector).
%     phase_estimate - Local fitted phase at each symbol, in radians.
%     slope          - Local fitted slope, in radians per symbol.
%     half_window    - Half-width of the centered window, in symbols.
%     window_start   - First symbol index of each fitting window.
%     window_end     - Last symbol index of each fitting window.
%   The existing degenerate-window convention sets power to zero when the
%   fit determinant is zero; fitted phase and slope remain undefined there.

    y = Rx_phi(1:cfg.oversampling_factor:end);
    [EEPN_power, phase_fit] = windowed_linear_fit(y, cfg);
end
