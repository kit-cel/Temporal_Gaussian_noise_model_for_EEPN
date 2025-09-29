# Temporal Gaussian Noise Model for Equalization-Enhanced Phase Noise
## Benedikt Geiger, Fred Buchali, Vahid Aref, and Laurent Schmalen

This repository contains an implementation of the **Temporal Gaussian Noise (TGN) Model** for **Equalization-Enhanced Phase Noise (EEPN)** proposed in [1], along with a full end-to-end system simulation as a reference.  
Implementations are provided in **MATLAB** (both as a Live Script and as a plain `.m` file) and in **Python** (as a Jupyter Notebook).  

---

## Minimal Model Implementation  

The **TGN Model** can be realized in just **three lines of code**:  

**MATLAB**
```matlab
% (1) Load or generate LO phase noise realization
Rx_phi = cumsum(sqrt(sigma2_LO) * randn(size(Tx_symbols),1)) + 2*pi*rand(1);

% (2) Calculate the time-varying distortion power
sigma_time_varying = system_noise_power + movvar(Rx_phi, CD_memory + 1);

% (3) Sample AWGN from time-varying distortion power and add to transmit signal
Rx_symbols = Tx_symbols + sqrt(time_varying/2) .* ...
             (randn(size(Tx_symbols)) + 1j*randn(size(Tx_symbols)));
```

**Python**
```python
\# (1) Load or generate LO phase noise realization
Rx_phi = np.cumsum(np.sqrt(sigma2_LO) * np.random.randn(len(Tx_symbols))) + 2*np.pi*np.random.rand()

\# (2) Calculate the time-varying distortion power
sigma_time_varying = system_noise_power + movvar(Rx_phi, CD_memory + 1)

\# (3) Sample AWGN from time-varying distortion power and add to transmit signal
Rx_symbols = Tx_symbols + (np.sqrt(sigma_time_varying/2)*np.random.randn(len(Tx_symbols)) + 1j*np.sqrt(sigma_time_varying/2)*np.random.randn(len(Tx_symbols)))
```

---

## Variable Definitions  

**`sigma2_LO`** — variance of the LO Wiener process  

$$
\sigma^2_{\text{LO}} = \frac{2\pi\cdot\text{linewidth}}{\text{oversampling factor}\cdot\text{symbol rate}}
$$

**`CD_memory`** — chromatic-dispersion–induced memory (in samples)  

$$
\mathrm{CD\_memory} = D_{CD}\cdot\text{fiber length}\cdot\frac{\lambda^2}{c_0}\cdot\text{symbol rate}^2
$$

**`system_noise_power`** — noise power accounting for ASE, fiber nonlinearity, and transceiver impairments

**`Tx_symbols`** — normalized transmit symbols (e.g., 16‑QAM) 

---

## Reference  

[1] B. Geiger, F. Buchali, V. Aref, and L. Schmalen, *“A temporal Gaussian noise model for equalization-enhanced phase noise,”*  
Proc. Eur. Conf. Opt. Commun. (ECOC), Copenhagen, Denmark, Sep. 2025.  
[arXiv:2507.08470](http://arxiv.org/abs/2507.08470)  
