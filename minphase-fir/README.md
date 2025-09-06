# Minimum-Phase FIR Filter Design in MATLAB

## Objective
Design and analyze FIR filters using multiple approaches, and transform a linear-phase equiripple FIR into a minimum-phase FIR while preserving magnitude response and reducing group delay.

## Methods
This project demonstrates three steps:

1. **Kaiser Window FIR**  
   - Use the Kaiser window method to design a low-pass FIR that meets passband/stopband ripple specifications.

2. **Equiripple FIR (Parks–McClellan)**  
   - Design an equiripple FIR filter using the minimax (Parks–McClellan) algorithm.  
   - Examine its impulse response, frequency response, and zero distribution.

3. **Minimum-Phase FIR Transformation**  
   - Modify the impulse response by lifting the center sample.  
   - Reflect zeros inside the unit circle to obtain a minimum-phase version.  
   - Normalize to compare magnitude response with the original equiripple filter.

## How to Run
From MATLAB:
```matlab
cd minphase-fir
run("main.m")
