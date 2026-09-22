# (Generalized) Temporal Gaussian Noise Model for Equalization-Enhanced Phase Noise
## Benedikt Geiger, Fred Buchali, Vahid Aref, and Laurent Schmalen

This repository contains an implementation of the **Temporal Gaussian Noise (TGN) Model** for **Equalization-Enhanced Phase Noise (EEPN)** proposed in [1], and the **Generalized Temporal Gaussian Noise (TGN) Model** for EEPN proposed in [2], along with a full end-to-end system simulation as a reference.

## Quickstart

Two scripts reproduce the results from the two papers:
- **ECOC_25** — the Temporal Gaussian noise model, only the CPR from [1].
- **ECOC_26** — generalizes the Temporal Gaussian noise model to all **four compensation cases** from [2] and reproduces its Fig. 3: distortion power over time and CCDF, model vs. full-system simulation, for all four cases side by side.

Each is provided in MATLAB (as a plain script and as an interactive Live Script) and in Python (as a Jupyter notebook):

| | MATLAB | Python |
|---|---|---|
| **ECOC_25** — 1 case (CPR), paper [1] | [`.m`](MATLAB/ECOC25_Temporal_Gaussian_Noise_Model_EEPN.m) · [`.mlx` Live Script](MATLAB/ECOC25_Temporal_Gaussian_Noise_Model_EEPN_Matlab_Live_Script.mlx) | [`.ipynb`](Python/ECOC25_Temporal_Gaussian_Noise_Model_EEPN.ipynb) |
| **ECOC_26** — 4 cases, paper [2] | [`.m`](MATLAB/ECOC26_Generalized_Temporal_Gaussian_Noise_Model_EEPN.m) · [`.mlx` Live Script](MATLAB/ECOC26_Generalized_Temporal_Gaussian_Noise_Model_EEPN_Live_Script.mlx) | [`.ipynb`](Python/ECOC26_Generalized_Temporal_Gaussian_Noise_Model_EEPN.ipynb) |

## Overview

EEPN arises since the chromatic dispersion compensation (CDC) filter interacts with the local oscillator (LO) phase noise which results in a burst-like SNR degradation. This repository contains code to calculate this time-varying EEPN distortion power and to simulate EEPN-impaired links. In particular, the LO phase fluctuations are transformed into frequency-dependent phase noise. Depending on the employed digital signal processing at the receiver, the frequency-dependent phase noise results in a frequency-dependent phase error. In particular, we distinguish four cases:

| Case | Order Ñ | Compensates | Compensation function | TGN model function |
|---|---|---|---|---|
| No compensation (**LO-PC**) | – | instantaneous LO phase φ_t | `LO_phase_cancellation` ([MATLAB](MATLAB/LO_phase_cancellation.m) · [Python](Python/LO_phase_cancellation.py)) | `calculate_EEPN_distortion_power_LO_phase_cancellation` ([MATLAB](MATLAB/calculate_EEPN_distortion_power_LO_phase_cancellation.m) · [Python](Python/calculate_EEPN_distortion_power_LO_phase_cancellation.py)) |
| **CPR** | 0 | constant phase offset a⁽⁰⁾ | `CPR` ([MATLAB](MATLAB/CPR.m) · [Python](Python/CPR.py)) | `calculate_EEPN_distortion_power_CPR` ([MATLAB](MATLAB/calculate_EEPN_distortion_power_CPR.m) · [Python](Python/calculate_EEPN_distortion_power_CPR.py)) |
| **Timing recovery** | 1 | a⁽⁰⁾ + a⁽¹⁾·f | `timing_recovery` ([MATLAB](MATLAB/timing_recovery.m) · [Python](Python/timing_recovery.py)) | `calculate_EEPN_distortion_power_timing_recovery` ([MATLAB](MATLAB/calculate_EEPN_distortion_power_timing_recovery.m) · [Python](Python/calculate_EEPN_distortion_power_timing_recovery.py)) |
| **Adaptive filtering** (higher-order compensation) | Ñ_AF | Σₙ a⁽ⁿ⁾·fⁿ | `full_compensation` ([MATLAB](MATLAB/full_compensation.m) · [Python](Python/full_compensation.py)) | `calculate_EEPN_distortion_power_full_compensation` ([MATLAB](MATLAB/calculate_EEPN_distortion_power_full_compensation.m) · [Python](Python/calculate_EEPN_distortion_power_full_compensation.py)) |

The compensation functions realize the **ideal, genie-aided** version of each case — i.e. they use the true LO phase directly, rather than a receiver-realistic estimate of it — to give a clean, consistent reference for validating the TGN model against a full-system simulation across all four cases.

## Minimal Model

All four cases follow the same three-step usage pattern:

**MATLAB**
```matlab
% (1) Load or generate LO phase noise realization
Rx_phi = cumsum(sqrt(sigma2_LO) * randn(size(Tx_symbols),1)) + 2*pi*rand(1);

% (2) Calculate the time-varying distortion power (see table below)
sigma_time_varying = system_noise_power + <case-specific term>;

% (3) Sample AWGN from time-varying distortion power and add to transmit signal
Rx_symbols = Tx_symbols + sqrt(sigma_time_varying/2) .* ...
             (randn(size(Tx_symbols)) + 1j*randn(size(Tx_symbols)));
```

**Python**
```python
# (1) Load or generate LO phase noise realization
Rx_phi = np.cumsum(np.sqrt(sigma2_LO) * np.random.randn(len(Tx_symbols))) + 2*np.pi*np.random.rand()

# (2) Calculate the time-varying distortion power (see table below)
sigma_time_varying = system_noise_power + <case-specific term>

# (3) Sample AWGN from time-varying distortion power and add to transmit signal
Rx_symbols = Tx_symbols + (np.sqrt(sigma_time_varying/2)*np.random.randn(len(Tx_symbols)) + 1j*np.sqrt(sigma_time_varying/2)*np.random.randn(len(Tx_symbols)))
```

**Step (2), `<case-specific term>`, per case:**

| Case | MATLAB | Python |
|---|---|---|
| **LO-PC** (closed form) | `movmean(Rx_phi.^2, CD_memory+1) + Rx_phi.^2 - 2*Rx_phi.*movmean(Rx_phi, CD_memory+1)` | `movmean(Rx_phi**2, CD_memory+1) + Rx_phi**2 - 2*Rx_phi*movmean(Rx_phi, CD_memory+1)` |
| **CPR** (closed form) | `movvar(Rx_phi, CD_memory + 1)` | `movvar(Rx_phi, CD_memory + 1)` |
| **Timing recovery** (windowed fit) | `calculate_EEPN_distortion_power_timing_recovery(Rx_phi, cfg)` | *(identical call)* |
| **Adaptive filtering** (windowed fit) | `calculate_EEPN_distortion_power_full_compensation(Rx_phi, cfg)` | *(identical call)* |

For **timing recovery** and **adaptive filtering** the time-varying distortion power requires a least-squares fit over the CD-memory window ([`windowed_linear_fit.m`](MATLAB/functions/windowed_linear_fit.m) / [`windowed_polynomial_fit.m`](MATLAB/functions/windowed_polynomial_fit.m)) which is hidden behind the two functions.

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

**`compensation_order`** (Ñ_AF) — polynomial order of the adaptive-filtering / higher-order compensation case

---

## Requirements

- **Python** — see [`Python/requirements.txt`](Python/requirements.txt) (`numpy`, `matplotlib`, `jupyter`).
- **MATLAB** — base MATLAB, plus the **Signal Processing Toolbox** (`fftfilt`, used by [`local_polynomial_fit_moments.m`](MATLAB/functions/local_polynomial_fit_moments.m) for the timing-recovery and adaptive-filtering cases).

---

## Reference

If you use this code, please cite the corresponding paper(s) below rather than this repository directly.

[1] B. Geiger, F. Buchali, V. Aref, and L. Schmalen, *“A temporal Gaussian noise model for equalization-enhanced phase noise,”*  
Proc. Eur. Conf. Opt. Commun. (ECOC), Copenhagen, Denmark, Sep. 2025.  
[arXiv:2507.08470](http://arxiv.org/abs/2507.08470)

[2] B. Geiger, F. Buchali, V. Aref, and L. Schmalen, *“Modeling and Mitigation of Equalization-Enhanced Phase Noise,”*  
Proc. Eur. Conf. Opt. Commun. (ECOC), Malaga, Spain, Sep. 2026. [arXiv:2606.21468](http://arxiv.org/abs/2606.21468)

Released under the MIT License (see [LICENSE](LICENSE)).

> Parts of this repository were generated with the assistance of [Claude Code](https://claude.com/claude-code), based on an existing, validated internal implementation, and reviewed by the authors.
