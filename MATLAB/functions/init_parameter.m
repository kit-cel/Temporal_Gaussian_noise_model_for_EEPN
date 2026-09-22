function cfg = init_parameter(cfg)
%INIT_PARAMETER calcualtes the parameters which are based on others
    cfg.lambda = 1.55e-6;                       % Center Wavelength in m
    cfg.c_0 = 3e8;                              % Vacuum speed of light in m/s

    cfg.t = [0:cfg.num_transmission_symbols-1]/cfg.symbol_rate; % Time vector

    cfg.df = cfg.symbol_rate/cfg.num_transmission_symbols; % frequncy resolution
    cfg.f = [-cfg.num_transmission_symbols*cfg.oversampling_factor/2:1:cfg.num_transmission_symbols*cfg.oversampling_factor/2-1]*cfg.df; % Frequency vector

    % Calculate the temporal broadening a pulse due to chromatic dispersion
    cfg.D = cfg.D_cd*cfg.fiber_length;                  % Temporal broadening/accumulated dispersion in ns/nm
    cfg.Delta_lambda = (cfg.lambda)^2/cfg.c_0*cfg.symbol_rate; % Sepctral width in wavelength
    cfg.Delta_T = cfg.D*cfg.Delta_lambda;               % Temporal broadening/accumulated dispersion in s
    cfg.CD_memory = round((cfg.Delta_T*cfg.symbol_rate));     % Temporal broadening/CD memory in samples

    % Variance of the Wiener process
    cfg.sigma2_LO = 2*pi*cfg.linewidth*1/(cfg.oversampling_factor*cfg.symbol_rate);

    cfg.system_noise_power = 10^(-cfg.snr/10);
end
