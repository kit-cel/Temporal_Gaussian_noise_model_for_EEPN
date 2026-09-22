function sigma = estimate_moving_error_power(Rx_symbols, Tx_symbols, cfg)
%ESTIMATE_MOVING_ERROR_POWER Estimate error power using moving variance.
%   Power is estimated as the variance of Rx_symbols-Tx_symbols within each
%   moving window, with the local mean removed by movvar.
%   Rx_symbols and Tx_symbols must be aligned vectors of the same size.
%   cfg.block_length_SNR_evaluation sets the moving variance window length.
%   The existing analysis convention discards cfg.discard_symbols_analysis
%   values at the start and one more than that at the end.

    sigma = movvar(Rx_symbols-Tx_symbols, cfg.block_length_SNR_evaluation);
    sigma(1:cfg.discard_symbols_analysis) = [];
    sigma(end-cfg.discard_symbols_analysis:end) = [];
end
