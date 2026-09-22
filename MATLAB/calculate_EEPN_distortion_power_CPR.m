function EEPN_power = calculate_EEPN_distortion_power_CPR(Rx_phi, cfg)
%CALCULATE_EEPN_DISTORTION_POWER_CPR Model EEPN power after ideal CPR.
%   Rx_phi is the real, unwrapped LO phase vector at the oversampled rate.
%   cfg.oversampling_factor determines the symbol-rate phase samples.
%   CPR removes a constant (zeroth-order) phase offset equal to the mean
%   LO phase within the CD-memory window (the IDR estimate; see
%   calculate_EEPN_distortion_power_full_compensation for the general
%   order-m case). Each output is the mean squared deviation of the
%   window phase from that local mean, i.e. the local phase variance.
%   The centered window has 2*ceil(cfg.CD_memory/2)+1 samples in the
%   interior and shrinks at the sequence boundaries, matching
%   calculate_EEPN_distortion_power_LO_phase_cancellation.
%   EEPN_power is a column vector with one value per symbol. Background
%   noise is excluded, and no samples are discarded from the output.
%
%   This formula is equivalent to the "Minimal Model Implementation"
%   one-liner in the README (movvar(Rx_phi, CD_memory+1)), just computed
%   via cumulative sums instead of MATLAB's movvar.

    y = Rx_phi(1:cfg.oversampling_factor:end);
    y = y(:);

    N = length(y);
    L = cfg.CD_memory + 1;
    half = ceil((L-1)/2);

    % Cumulative sums allow constant-time evaluation of each window sum.
    cs1 = cumsum([0; y]);
    cs2 = cumsum([0; y.^2]);

    % Shrink the centered window at the sequence boundaries.
    idx_start = max(1, (1:N)' - half);
    idx_end = min(N, (1:N)' + half);

    sum_y = cs1(idx_end + 1) - cs1(idx_start);
    sum_y2 = cs2(idx_end + 1) - cs2(idx_start);
    win_len = idx_end - idx_start + 1;

    m1 = sum_y ./ win_len;
    m2 = sum_y2 ./ win_len;

    % Mean of (window phase - window mean).^2.
    EEPN_power = m2 - m1.^2;
end
