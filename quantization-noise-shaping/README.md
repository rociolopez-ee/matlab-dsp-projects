# Quantization Effects in a First-Order IIR Filter

## Objective
Simulate finite word-length effects in digital filters and evaluate the impact of quantization noise on a first-order IIR system.  
The goal is to measure the signal-to-quantization-noise ratio (SQNR) under different pole locations and compare empirical results with theoretical predictions.

## Methods
This project implements:
- A test signal consisting of two sinusoidal components  
- Fixed-point quantization with 8-bit resolution  
- A first-order IIR filter with feedback, tested at different pole locations 
- Computation of:
  - Quantized vs. ideal outputs  
  - Error signal and noise variance  
  - SQNR (empirical vs. theoretical)  
- Spectrum analysis of the quantized output

## How to Run
From MATLAB:
```matlab
cd quantization-noise
run("main.m")
