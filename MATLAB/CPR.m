function [Rx_symbols_CPR, phase_est] = CPR(Rx_symbols, Rx_phi, cfg)
%CPR Carrier phase recovery (ideal, genie-aided order-0 FDPE compensation).
%   Rx_symbols is the symbol-rate receive signal (after CD compensation,
%   matched filtering and downsampling). Rx_phi is the real, unwrapped LO
%   phase vector at the oversampled rate; cfg.oversampling_factor
%   downsamples it to the symbol rate used here.
%
%   This is the ideal, genie-aided phase estimate: the local mean of the
%   true LO phase, averaged over the CD-memory window (cfg.CD_memory+1
%   symbols) -- the same window as calculate_EEPN_distortion_power_CPR's
%   model, so this is the direct empirical counterpart of that model
%   rather than a receiver-realistic design. Not realizable by an actual
%   receiver, which only has access to the received symbols, not Rx_phi;
%   see [2] for practical (data-aided/blind) estimators. CPR corresponds
%   to compensating the zeroth-order FDPE term a^(0)_t in the
%   receiver-structure table of [2].
%
%   Rx_symbols_CPR is a column vector of symbol-rate compensated symbols.
%   Optional phase_est returns the estimated phase offset per symbol.

    Rx_symbols = Rx_symbols(:);
    Rx_phi_symbols = Rx_phi(1:cfg.oversampling_factor:end);
    phase_est = movmean(Rx_phi_symbols(:), cfg.CD_memory+1);
    Rx_symbols_CPR = Rx_symbols .* exp(-1j*phase_est);
end
