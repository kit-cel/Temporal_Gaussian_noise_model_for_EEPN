function [EEPN_power, phase_fit] = windowed_polynomial_fit(y, order, cfg)
%WINDOWED_POLYNOMIAL_FIT Local order-m least-squares fit over the
%CD-memory window, for any symbol-rate phase-like sequence.
%   y is a symbol-rate column vector: the true LO phase, or a data-
%   estimated proxy for it. order is the polynomial order m. cfg.CD_memory
%   sets the window.
%
%   This is the shared engine behind
%   calculate_EEPN_distortion_power_full_compensation and
%   full_compensation.m, the order-m analogue of windowed_linear_fit.m
%   (order 1).
%
%   See calculate_EEPN_distortion_power_full_compensation for the meaning
%   of EEPN_power and phase_fit (phase_symbols, coefficients, order,
%   half_window). coefficients is (order+1)-by-N, row p+1 holding the
%   coefficient of the window-offset^p term (p = 0..order), in powers of
%   the raw window offset k in [-half_window, half_window].

    y = y(:);

    L = cfg.CD_memory + 1;
    half = ceil((L-1)/2);

    y_padded = [repmat(y(1), half, 1); y; repmat(y(end), half, 1)];
    [S, Syy] = local_polynomial_fit_moments(y_padded, half, order);

    % Solve the normal equations in the offset x = k/half in [-1,1] rather
    % than the raw offset k in [-half,half]: for half in the thousands and
    % order above a few, the raw-offset Vandermonde Gram matrix is
    % catastrophically ill-conditioned (entries span k^0 .. k^(2*order)).
    % S scales linearly with the kernel, so S in the normalized offset is
    % simply S ./ half.^p; no second FFT-based pass is needed.
    powers = (0:order).';
    S_norm = S ./ (half.^powers);

    x = (-half:half).'/half;
    V = x.^(0:order);
    G = V.' * V;

    beta_norm = G \ S_norm;
    SSE = Syy - sum(beta_norm .* S_norm, 1);

    EEPN_power = (SSE ./ L).';

    if nargout > 1
        phase_fit.phase_symbols = y;
        % Convert back to coefficients of the raw window offset k, to
        % match the convention consumed by apply_block_fdpe_filter.
        phase_fit.coefficients = beta_norm ./ (half.^powers);
        phase_fit.order = order;
        phase_fit.half_window = half;
    end
end
