function EEPN_power = calculate_EEPN_distortion_power_LO_phase_cancellation(Rx_phi, cfg)
%CALCULATE_EEPN_DISTORTION_POWER_LO_PHASE_CANCELLATION Model EEPN power after LO-PC.
%   Rx_phi is the real, unwrapped LO phase vector at the oversampled rate.
%   cfg.oversampling_factor determines the symbol-rate phase samples.
%   Each output is the mean squared phase difference relative to the phase
%   at the current symbol, evaluated over a centered CD-memory window.
%   The existing window convention uses 2*ceil(cfg.CD_memory/2)+1 samples
%   in the interior and shrinks the window at the sequence boundaries.
%   EEPN_power is a column vector with one value per symbol. Background
%   noise is excluded, and no samples are discarded from the output.
%
%   LO-PC removes the instantaneous LO phase itself, not a windowed fit,
%   so it is not part of the CPR/timing-recovery/full-compensation
%   polynomial-order family (see calculate_EEPN_distortion_power_CPR,
%   calculate_EEPN_distortion_power_timing_recovery, and
%   calculate_EEPN_distortion_power_full_compensation for orders 0, 1, m).

    x = Rx_phi(1:cfg.oversampling_factor:end);
    x = x(:);

    N = length(x);
    L = cfg.CD_memory + 1;
    half = ceil((L-1)/2);

    % Cumulative sums allow constant-time evaluation of each window sum.
    cs1 = cumsum([0; x]);
    cs2 = cumsum([0; x.^2]);

    % Shrink the centered window at the sequence boundaries.
    idx_start = max(1, (1:N)' - half);
    idx_end = min(N, (1:N)' + half);

    sum_x = cs1(idx_end + 1) - cs1(idx_start);
    sum_x2 = cs2(idx_end + 1) - cs2(idx_start);
    win_len = idx_end - idx_start + 1;

    m1 = sum_x ./ win_len;
    m2 = sum_x2 ./ win_len;

    % Mean of (window phase - current phase).^2.
    EEPN_power = m2 + x.^2 - 2*x.*m1;
end
