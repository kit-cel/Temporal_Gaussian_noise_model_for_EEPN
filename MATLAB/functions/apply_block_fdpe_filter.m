function Rx_symbols_out = apply_block_fdpe_filter(Rx_symbols_derotated, coefficients, half, cfg)
%APPLY_BLOCK_FDPE_FILTER Realize a higher-order FDPE compensation filter.
%   Rx_symbols_derotated are symbol-rate symbols with the local
%   zeroth-order (constant) phase term already removed by the caller
%   (e.g. by derotating with the local polynomial fit's own value at
%   each symbol; see timing_recovery.m and full_compensation.m).
%
%   coefficients is an order-by-N matrix of local polynomial fit
%   coefficients for orders p = 1..order (row p holds the coefficient of
%   the window-offset^p term; N must match length(Rx_symbols_derotated)).
%   half is the half-width, in symbols, of the CD-memory window that the
%   coefficients were fitted over (see calculate_EEPN_distortion_power_
%   timing_recovery and calculate_EEPN_distortion_power_full_compensation).
%
%   Since only cfg.compensation_filter_taps one-sided taps are used to
%   approximate the (CD-memory-wide) ideal all-pass compensation filter,
%   the coefficients are evaluated on a coarser, rescaled offset grid
%   (scale = half/taps) so that the short filter spans the same
%   frequency support as the full-width ideal filter. The fitted
%   coefficients are assumed locally constant over one
%   cfg.compensation_block_length-symbol block, and processing proceeds
%   via overlap-add with overlap cfg.compensation_block_overlap.
%
%   Rx_symbols_out has the same length as Rx_symbols_derotated. Samples
%   before the first block, in its guard region, or past the last full
%   block are left at zero; these fall inside the region that
%   cfg.discard_symbols_analysis removes from every analysis in this
%   codebase and must not be used directly.

    Rx_symbols_derotated = Rx_symbols_derotated(:);
    N = length(Rx_symbols_derotated);
    order = size(coefficients, 1);

    taps = cfg.compensation_filter_taps;
    block_length = cfg.compensation_block_length;
    overlap = round(cfg.compensation_block_overlap * block_length);
    hop = block_length - overlap;
    guard = floor(overlap/2);

    starts = 1:hop:(N - block_length + 1);
    n = (-taps:taps).';
    scale = half/taps;

    Rx_symbols_out = zeros(N, 1);
    for iter = 1:numel(starts)
        start_idx = starts(iter);
        end_idx = start_idx + block_length - 1;
        center_idx = round((start_idx + end_idx)/2);

        write_start = start_idx + guard;
        write_end = write_start + hop - 1;
        if write_end > N
            continue
        end

        phase_ramp = zeros(size(n));
        for p = 1:order
            phase_ramp = phase_ramp + coefficients(p, center_idx)*(n*scale).^p;
        end
        phase_ramp = flip(phase_ramp);

        comp_filter = fftshift(ifft(ifftshift(exp(-1j*phase_ramp))));
        block_out = conv(Rx_symbols_derotated(start_idx:end_idx), comp_filter, 'same');
        Rx_symbols_out(write_start:write_end) = block_out((guard+1):(guard+hop));
    end
end
