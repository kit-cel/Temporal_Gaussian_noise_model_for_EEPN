function [Rx_symbols_full, EEPN_power, phase_fit] = full_compensation(Rx_symbols, Rx_phi, cfg)
%FULL_COMPENSATION Reverse the local order-m FDPE (ideal, genie-aided
%adaptive filtering).
%   Rx_symbols is the symbol-rate receive signal (after CD compensation,
%   matched filtering and downsampling). Rx_phi is the real, unwrapped LO
%   phase vector at the oversampled rate. cfg.compensation_order sets the
%   polynomial order m (Ntilde in the paper notation, [2]).
%   cfg.compensation_filter_taps, cfg.compensation_block_length and
%   cfg.compensation_block_overlap configure the finite-tap FIR
%   realization of the order 1..m terms; see private/apply_block_fdpe_filter.m.
%
%   As in timing_recovery.m (the order-1 special case), this is one
%   receiver structure of order m that jointly compensates a0..am,
%   theta^comp_t(f) = sum_{n=0}^m a^(n)_t*f^n (see [2]). This is realized
%   here as one joint private/windowed_polynomial_fit.m fit of the true
%   LO phase Rx_phi (genie-aided, not realizable by an actual receiver;
%   see [2] for practical data-aided/blind estimators). "Adaptive
%   filtering" with m = Ntilde_AF is the "higher-order compensation" case
%   of the receiver-structure table in [2].
%
%   Rx_symbols_full is a column vector of symbol-rate compensated
%   symbols; samples inside the first/last block are left at zero (see
%   apply_block_fdpe_filter) and fall inside the region that
%   cfg.discard_symbols_analysis removes from every analysis in this
%   codebase. EEPN_power is the model's instantaneous distortion power
%   (matches calculate_EEPN_distortion_power_full_compensation). phase_fit
%   exposes the coefficients used for the derotation and filter
%   construction.

    order = cfg.compensation_order;

    phase_sequence = Rx_phi(1:cfg.oversampling_factor:end);
    [EEPN_power, fit] = windowed_polynomial_fit(phase_sequence, order, cfg);
    a0 = fit.coefficients(1, :).';
    coefficients = fit.coefficients(2:end, :);
    half = fit.half_window;

    Rx_symbols_order0_compensated = Rx_symbols(:) .* exp(-1j*a0(:));
    Rx_symbols_full = apply_block_fdpe_filter(Rx_symbols_order0_compensated, coefficients, half, cfg);

    if nargout > 2
        phase_fit.coefficients = [a0(:).'; coefficients];
        phase_fit.order = order;
        phase_fit.half_window = half;
    end
end
