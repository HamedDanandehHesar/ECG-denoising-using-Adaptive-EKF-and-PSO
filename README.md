# ECG-denoising-using-Adaptive-EKF-and-PSO
ECG denoising using Adaptive EKF and PSO

MATLAB implementation of an **adaptive Bayesian (Adaptive EKF) ECG denoising algorithm** based on **nonlinear Kalman filtering with adaptive covariance estimation**, inspired by the paper:

**Efficient Bayesian ECG denoising using adaptive covariance estimation and nonlinear Kalman filtering**  
Hamed Danandeh Hesar, Amin Danandeh Hesar  
Computers & Electrical Engineering, 2024.

This implementation removes noise from ECG signals by combining:

- ECG phase-domain modeling
- Gaussian mixture representation of ECG morphology
- Adaptive Extended Kalman Filter (AEKF)
- Adaptive Extended Kalman Smoother (AEKS)

The algorithm estimates the clean ECG signal from noisy measurements while **adapting the process and measurement covariance matrices online**.

---

# Overview

Electrocardiogram (ECG) signals are often contaminated by noise such as:

- muscle artifacts
- baseline wander
- powerline interference
- sensor noise

This project implements a **Bayesian filtering framework** that models the ECG signal as a **sum of Gaussian kernels in the phase domain** and performs denoising using:

1. **R‑peak detection**
2. **Phase-domain ECG modeling**
3. **Gaussian mixture approximation of ECG morphology**
4. **Adaptive Extended Kalman Filtering**
5. **Kalman smoothing**

---

# Main Processing Pipeline

The algorithm follows these steps:

### 1. Load ECG signal
A `.mat` file containing ECG data is loaded.

Required variables inside the file:

- `x` → ECG signal matrix  
- `fs` → sampling frequency (Hz)

The first channel is used as the ECG signal.

---

### 2. Add Noise

White Gaussian noise is added to the ECG signal to simulate a noisy measurement.

Example parameter:

```
SNR = 6 dB
```

---

### 3. R‑Peak Detection

The **Pan–Tompkins algorithm** is used to detect R‑peaks:

```
[qrs_positions] = pantompkins_qrs(signal, fs)
```

Detected peaks are visualized on the noisy ECG.

---

### 4. ECG Phase Calculation

The ECG signal is transformed into a **phase representation** based on RR intervals.

Each cardiac cycle is mapped to the phase range:

```
[-π , π]
```

This phase representation simplifies modeling of the ECG morphology.

---

### 5. Mean ECG Extraction

The algorithm estimates:

- mean ECG waveform in the phase domain
- standard deviation of ECG amplitude per phase bin

The mean ECG is then **further smoothed using wavelet denoising**.

---

### 6. Gaussian Mixture Modeling

The ECG morphology is modeled as a **sum of Gaussian kernels**:

``` math
\[
ECG(\theta) = \sum_{i=1}^{L} a_i \exp\left(-\frac{(\theta-\theta_i)^2}{2b_i^2}\right)
\]```
Where:

- \(a_i\) → Gaussian amplitude  
- \(b_i\) → Gaussian width  
- \( \theta_i \) → Gaussian center  

Parameters are estimated using **Particle Swarm Optimization (PSO)**.

Only the **strongest Gaussian components** are kept.

---

### 7. Adaptive Extended Kalman Filter (AEKF)

The state vector contains:

- ECG amplitude
- ECG phase
- instantaneous angular frequency

The filter:

- estimates the clean ECG signal
- tracks the phase evolution
- adapts covariance matrices \(Q\) and \(R\)

Adaptive covariance estimation improves robustness to noise.

---

### 8. Adaptive Covariance Estimation

Measurement noise covariance **R** and process noise covariance **Q** are updated using:

- sliding window estimation
- forgetting factor

This allows the filter to adapt to changing signal conditions.

---

### 9. Extended Kalman Smoothing (EKS)

After filtering, a **backward smoothing step** refines the ECG estimate.

The smoother produces a more accurate reconstruction than the filter alone.

---

# Input

The program requires a MATLAB `.mat` file containing:

```
x   → ECG signal matrix
fs  → sampling frequency
```

Example structure:

```
x = [ECG_signal]
fs = 360
```

When running the script, MATLAB will ask the user to select the ECG file.

---

# Output

The algorithm produces:

### 1. Denoised ECG (Adaptive EKF)

Filtered signal estimate.

### 2. Smoothed ECG (Adaptive EKS)

Improved estimate after backward smoothing.

### 3. SNR Improvement Metrics

The script computes SNR improvement:

```
AEKF_SNR
AEKS_SNR
```

---

# Visualization

The script generates several figures:

### Figure 1
Noisy ECG with detected R‑peaks.

### Figure 2
Comparison between:

- Original ECG
- Noisy ECG
- Adaptive EKF output

### Figure 3
Comparison between:

- Original ECG
- Noisy ECG
- Adaptive EKS output

---

# Key Parameters

Important parameters that can be modified:

```
SNR                % noise level
ecg_bins           % number of phase bins
MaxNumGaussian     % maximum Gaussian kernels
window_size        % covariance adaptation window
forgetting_factor  % adaptive memory factor
```

These parameters affect denoising performance.

---

# Required Functions

The following external functions must exist in the MATLAB path:

- `pantompkins_qrs.m` → R‑peak detection
- `particleswarm` → MATLAB Global Optimization Toolbox

---

# MATLAB Toolboxes Required

- Signal Processing Toolbox
- Wavelet Toolbox
- Global Optimization Toolbox

---

# Applications

This algorithm can be used for:

- Biomedical signal processing
- ECG denoising research
- wearable health monitoring
- biomedical machine learning pipelines

---

# Reference

H. Danandeh Hesar, A. Danandeh Hesar,  
**Efficient Bayesian ECG denoising using adaptive covariance estimation and nonlinear Kalman filtering**,  
Computers & Electrical Engineering, 2024.

DOI:
https://doi.org/10.1016/j.compeleceng.2024.109869

---

