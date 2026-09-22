function [S, Syy] = local_polynomial_fit_moments(y, half, order)
%LOCAL_POLYNOMIAL_FIT_MOMENTS Windowed moments for local polynomial regression.
%   [S, Syy] = LOCAL_POLYNOMIAL_FIT_MOMENTS(y, half, order) computes, for
%   every sample ell with a fully interior centered window
%   (half+1 <= ell <= length(y)-half), the moments
%       S(p+1, ell-half) = sum_{k=-half}^{half} k^p * y(ell+k),  p = 0..order
%       Syy(ell-half)    = sum_{k=-half}^{half} y(ell+k)^2
%   using FFT-based FIR filtering (fftfilt), which stays tractable even
%   for window half-widths of several thousand samples (a direct
%   sliding-window sum or per-window polyfit would not).
%
%   S is (order+1)-by-(length(y)-2*half), one row per moment order p and
%   one column per valid center ell, in ascending order of ell. Syy is a
%   row vector with the same number of columns.
%
%   Callers with a signal shorter than 2*half+1 samples, or that need
%   values near the sequence boundary, must pad or otherwise extend y
%   themselves; this function only ever returns fully interior windows.

    y = y(:);
    N = length(y);
    L = 2*half + 1;
    interior_idx = (L:N).';

    j = (half:-1:-half).';   % matches the fftfilt/filter causal-shift convention
    S = zeros(order+1, numel(interior_idx));
    for p = 0:order
        filtered = fftfilt(j.^p, y);
        S(p+1, :) = filtered(interior_idx).';
    end

    filtered_yy = fftfilt(ones(L, 1), y.^2);
    Syy = filtered_yy(interior_idx).';
end
