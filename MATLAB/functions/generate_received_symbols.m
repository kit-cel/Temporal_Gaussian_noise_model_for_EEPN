function [Rx_symbols, Tx_symbols, Rx_phi, Rx_symbols_up] = generate_received_symbols(cfg)
%GENERATE_RECEIVED_SYMBOLS Simulate the Tx/channel/Rx-frontend chain.
%   Generates random Tx symbols, pulse-shapes and passes them through a
%   chromatic-dispersion channel with AWGN, applies LO phase noise, then
%   performs CD compensation, matched filtering and downsampling back to
%   symbol rate. Requires cfg already populated by init_parameter (uses
%   cfg.f, cfg.lambda, cfg.c_0, cfg.sigma2_LO, cfg.system_noise_power)
%   plus cfg.mod_order, cfg.num_transmission_symbols, cfg.oversampling_factor,
%   cfg.RRC_roll_off, cfg.RRC_span, cfg.D_cd, cfg.fiber_length.
%
%   Outputs:
%     Rx_symbols    - symbol-rate receive signal after CD compensation,
%                     matched filtering and downsampling.
%     Tx_symbols    - the transmitted (ground-truth) symbols.
%     Rx_phi        - oversampled, unwrapped LO phase (cumulative Wiener
%                     process + random initial phase).
%     Rx_symbols_up - the same signal as Rx_symbols before downsampling
%                     (cfg.oversampling_factor samples/symbol); optional.

    rrcFilter = rcosdesign(cfg.RRC_roll_off, cfg.RRC_span, cfg.oversampling_factor).';
    H_cd = exp(1j*pi*cfg.lambda^2/cfg.c_0*cfg.D_cd*cfg.fiber_length*cfg.f.^2).';

    Tx_symbols = qammod(randi(cfg.mod_order, cfg.num_transmission_symbols, 1)-1, cfg.mod_order, "UnitAveragePower", true);
    Tx_symbols_up = upsample(Tx_symbols, cfg.oversampling_factor);
    Tx_signal = conv(Tx_symbols_up, rrcFilter, 'same');

    dispersed_signal = ifft(fftshift(H_cd.*fftshift(fft(Tx_signal))));
    noise = sqrt(cfg.system_noise_power/2)*randn(size(Tx_signal)) + 1j*sqrt(cfg.system_noise_power/2)*randn(size(Tx_signal));
    Rx_signal = dispersed_signal + noise;

    Rx_delta_phi = sqrt(cfg.sigma2_LO)*randn(cfg.num_transmission_symbols*cfg.oversampling_factor, 1);
    Rx_phi = cumsum(Rx_delta_phi) + 2*pi*rand(1);

    Rx_signal_with_phase_noise = Rx_signal.*exp(1j*Rx_phi);
    Rx_signal_cdc = ifft(fftshift(fftshift(fft(Rx_signal_with_phase_noise)).*conj(H_cd)));
    Rx_symbols_up = conv(Rx_signal_cdc, rrcFilter, 'same');
    Rx_symbols = Rx_symbols_up(1:cfg.oversampling_factor:end);
end
