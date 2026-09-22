function Rx_symbols_LO_PC = LO_phase_cancellation(Rx_symbols, Rx_phi, cfg)
%LO_PHASE_CANCELLATION Reverse the instantaneous LO phase walk-off (genie-aided).
%   Rx_symbols is the symbol-rate receive signal (after CD compensation,
%   matched filtering and downsampling). Rx_phi is the real, unwrapped LO
%   phase vector at the oversampled rate; cfg.oversampling_factor
%   downsamples it to the symbol rate used here, matching the other
%   compensation functions (CPR, timing_recovery, full_compensation),
%   which all take the already-downsampled Rx_symbols.
%
%   LO-PC removes the LO's instantaneous phase itself (not a windowed
%   fit), so it is frequency-flat and is realized exactly as a
%   sample-wise derotation; no filtering is required. See
%   calculate_EEPN_distortion_power_LO_phase_cancellation for the
%   corresponding theoretical distortion power. This is the "no
%   compensation" case of the receiver-structure table in [2] (order -,
%   removes only the LO phase phi_t itself, none of the resulting FDPE).
%
%   Rx_symbols_LO_PC is a column vector of symbol-rate compensated symbols.

    Rx_phi_symbols = Rx_phi(1:cfg.oversampling_factor:end);
    Rx_symbols_LO_PC = Rx_symbols(:) .* exp(-1j*Rx_phi_symbols(:));
end
