function [Rx_symbols_TR, EEPN_power, phase_fit] = timing_recovery(Rx_symbols, Rx_phi, cfg)
%TIMING_RECOVERY Reverse the local order-1 (linear) FDPE (ideal, genie-aided).
%   Rx_symbols is the symbol-rate receive signal (after CD compensation,
%   matched filtering and downsampling). Rx_phi is the real, unwrapped LO
%   phase vector at the oversampled rate.
%   cfg.compensation_filter_taps, cfg.compensation_block_length and
%   cfg.compensation_block_overlap configure the finite-tap FIR
%   realization of the linear-phase (i.e. timing) correction; see
%   private/apply_block_fdpe_filter.m.
%
%   Per the receiver-structure table in [2], timing recovery is a single
%   receiver structure of order 1: it jointly compensates a0 AND a1,
%   theta^comp_t(f) = a0_t + a1_t*f. This is realized here as one joint
%   private/windowed_linear_fit.m fit of the true LO phase Rx_phi
%   (genie-aided, not realizable by an actual receiver; see [2] for
%   practical data-aided/blind estimators).
%
%   Rx_symbols_TR is a column vector of symbol-rate compensated symbols;
%   samples inside the first/last block are left at zero (see
%   apply_block_fdpe_filter) and fall inside the region that
%   cfg.discard_symbols_analysis removes from every analysis in this
%   codebase. EEPN_power is the model's instantaneous distortion power
%   (matches calculate_EEPN_distortion_power_timing_recovery). phase_fit
%   exposes the a0/a1 values used for the derotation and filter
%   construction.

    phase_sequence = Rx_phi(1:cfg.oversampling_factor:end);
    [EEPN_power, fit] = windowed_linear_fit(phase_sequence, cfg);
    a0 = fit.phase_estimate;
    a1 = fit.slope;
    half = fit.half_window;

    Rx_symbols_order0_compensated = Rx_symbols(:) .* exp(-1j*a0(:));
    coefficients = a1(:).';
    Rx_symbols_TR = apply_block_fdpe_filter(Rx_symbols_order0_compensated, coefficients, half, cfg);

    if nargout > 2
        phase_fit.phase_estimate = a0(:);
        phase_fit.slope = a1(:);
        phase_fit.half_window = half;
    end
end
