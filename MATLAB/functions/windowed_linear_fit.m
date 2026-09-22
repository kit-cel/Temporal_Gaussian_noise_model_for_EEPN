function [EEPN_power, phase_fit] = windowed_linear_fit(y, cfg)
%WINDOWED_LINEAR_FIT Local order-1 least-squares fit over the CD-memory
%window, for any symbol-rate phase-like sequence.
%   y is a symbol-rate column vector: the true LO phase, or a data-
%   estimated proxy for it. cfg.CD_memory sets the window.
%
%   This is the shared engine behind
%   calculate_EEPN_distortion_power_timing_recovery and timing_recovery.m.
%
%   See calculate_EEPN_distortion_power_timing_recovery for the meaning
%   of EEPN_power and phase_fit (phase_symbols, phase_estimate, slope,
%   half_window, window_start, window_end).

    y = y(:);
    N = length(y);

    L = cfg.CD_memory+1;
    half = ceil((L-1)/2);
    i = (1:N)';

    % Shrink the centered window at the sequence boundaries.
    s = max(1, i-half);
    e = min(N, i+half);
    m = e-s+1;

    % Cumulative sums allow constant-time evaluation of each window sum.
    cs_y = cumsum([0; y]);
    cs_y2 = cumsum([0; y.^2]);
    cs_i = cumsum([0; i]);
    cs_i2 = cumsum([0; i.^2]);
    cs_iy = cumsum([0; i.*y]);

    Sy = cs_y(e+1)-cs_y(s);
    Syy = cs_y2(e+1)-cs_y2(s);
    Si = cs_i(e+1)-cs_i(s);
    Sii = cs_i2(e+1)-cs_i2(s);
    Siy = cs_iy(e+1)-cs_iy(s);

    % Solve the normal equations for a0 + a1*i in each window.
    D = m.*Sii-Si.^2;
    a0 = (Sii.*Sy-Si.*Siy)./D;
    a1 = (m.*Siy-Si.*Sy)./D;

    % Normalize by the number of phase samples, not by residual degrees of freedom.
    SSE = Syy-(a0.*Sy+a1.*Siy);
    EEPN_power = SSE./m;
    EEPN_power(D == 0) = 0;

    if nargout > 1
        phase_fit.phase_symbols = y;
        phase_fit.phase_estimate = a0+a1.*i;
        phase_fit.slope = a1;
        phase_fit.half_window = half;
        phase_fit.window_start = s;
        phase_fit.window_end = e;
    end
end
