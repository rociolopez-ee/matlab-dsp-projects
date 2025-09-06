# Spectral Estimation in MATLAB

## Objective
Implement and compare different spectral estimation techniques for a signal with colored noise.  
The goal is to evaluate resolution, variance, and bias trade-offs between classical and modern PSD estimation methods.

## Methods
This project includes implementations of:
- **Periodogram** (single and ensemble-averaged)
- **Blackman–Tukey method**
- **Welch’s method** (with/without overlap)
- **Minimum Variance Distortionless Response (MVDR)**
- **Multitaper method**

All algorithms are coded in MATLAB using synthetic signals (Gaussian noise passed through an IIR filter).

## How to Run
From MATLAB:
```matlab
cd spectral-estimation
run("main.m")

