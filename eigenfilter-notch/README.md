# FIR Filter Design Using the Eigenfilter Method with Notches

## Objective
Design a linear-phase FIR filter using the eigenfilter method, and extend it by introducing spectral notches at specified frequencies.  
The goal is to demonstrate how the eigenfilter framework can be used to satisfy custom frequency-domain constraints such as enforced zeros.

## Methods
This project implements two stages:

1. **Unconstrained Eigenfilter Design**  
   - Formulate the FIR filter design as a quadratic minimization problem.  
   - Represent the amplitude response as a cosine expansion.  
   - Construct the extended matrix \( Q_t \) and solve via eigenvalue decomposition.  
   - Reconstruct the symmetric impulse response \( h[n] \).  

2. **Constrained Eigenfilter with Notches**  
   - Introduce linear equality constraints to enforce notches at 4500 Hz and 8000 Hz.  
   - Use null-space projection to incorporate constraints without redesigning the entire filter.  
   - Compare the constrained and unconstrained responses.

## How to Run
From MATLAB:
```matlab
cd eigenfilter-notch
run("main.m")
