
<h1 align="center"> 📊 General Transfer Function Analysis in C++</h1>
<h2 align="center">## 🚀 Project Overview</h2>
<p align="center">

Prueva2.cpp is a C++ implementation that allows you to compute poles, zeros, magnitude in dB, and phase of a general transfer function. The user can also specify a custom sigma value (σ) to perform the desired adjustment to the transfer function's analysis.

This implementation uses Eigen for matrix manipulations and polynomial root-finding. Theuser can provide the numerator and denominator coefficients of the transfer function, specify angular frequency ranges, and optionally input a value for σ.

## 🛠️ Features
✅ Calculate the poles and zeros of the transfer function from polynomial coefficients, including the option to incorporate a custom sigma value if necessary.
✅ Calculate the magnitude in dB at specified angular frequencies, with the option to include a custom sigma value if needed.
✅ Compute phase response in degrees at specified angular frequencies.
✅ Supports the analysis of general transfer functions with user-provided numerator and denominator polynomials.

## 🛠️ Dependencies Required

This script relies on several external libraries:

Standard C++ Libraries:

These are part of the C++ standard and come with most compilers.
Examples: <iostream>, <string>, <cmath>, <vector>, <complex>, and <regex>.

Eigen Library:

Purpose: A C++ template library for linear algebra operations (matrix manipulation, eigenvalues, etc.).
Installation Path: eigen-3.4.0/eigen-3.4.0/Eigen/Dense
Download from: Eigen GitLab Repository

# Example_1
Introduce the numerator (e.g., s^2 + 3s + 2):
s^2 + 3s + 2
Introduce the denominator (e.g., s^2 + 5s + 6):
s^2 + 5s + 6
Introduce the minimum frequency (freqMin):
3
Introduce the maximum frequency (freqMax):
6
Introduce the sigma value (leave empty or enter 0 for s = jw):
0
Assigning coefficient[2] = 2 to companion matrix.
Assigning coefficient[1] = 3 to companion matrix.
(-1,0)
(-2,0)
Assigning coefficient[2] = 6 to companion matrix.
Assigning coefficient[1] = 5 to companion matrix.
(-2,0)
(-3,0)
Zeros (Roots of the Numerator): (-1,0) (-2,0)
Poles (Roots of the Denominator): (-2,0) (-3,0)

For angular frequency w = 18.8496 rad/s:
Transfer Function |H(jw)| : (0.9835,0.1035)
Magnitude (dB): -0.0964
Debugging: Phase (radians) = 0.1048
Debugging: Phase (degrees) = 6.0063

For angular frequency w = 25.1327 rad/s:
Transfer Function |H(jw)| : (0.9906,0.0785)
Magnitude (dB): -0.0546
Debugging: Phase (radians) = 0.0790
Debugging: Phase (degrees) = 4.5284

For angular frequency w = 31.4159 rad/s:
Transfer Function |H(jw)| : (0.9940,0.0631)
Magnitude (dB): -0.0350
Debugging: Phase (radians) = 0.0634
Debugging: Phase (degrees) = 3.6316

For angular frequency w = 37.6991 rad/s:
Transfer Function |H(jw)| : (0.9958,0.0527)
Magnitude (dB): -0.0244
Debugging: Phase (radians) = 0.0529
Debugging: Phase (degrees) = 3.0304

# Example_2

Introduce the numerator (e.g., s^2 + 3s + 2):
s^2 + 3s + 2
Introduce the denominator (e.g., s^2 + 5s + 6):
s^2 + 5s + 6
Introduce the minimum frequency (freqMin):
3
Introduce the maximum frequency (freqMax):
6
Introduce the sigma value (leave empty or enter 0 for s = jw):
4
Assigning coefficient[2] = 2 to companion matrix.
Assigning coefficient[1] = 3 to companion matrix.
(3,0)
(2,0)
Assigning coefficient[2] = 6 to companion matrix.
Assigning coefficient[1] = 5 to companion matrix.
(2,0)
(1,0)
Zeros (Roots of the Numerator): (3,0) (2,0)
Poles (Roots of the Denominator): (2,0) (1,0)

For angular frequency w = 18.8496 rad/s:
Transfer Function |H(jw)| : (0.9654,0.0932)
Magnitude (dB): -0.2658
Debugging: Phase (radians) = 0.0963
Debugging: Phase (degrees) = 5.5170

For angular frequency w = 25.1327 rad/s:
Transfer Function |H(jw)| : (0.9794,0.0738)
Magnitude (dB): -0.1559
Debugging: Phase (radians) = 0.0753
Debugging: Phase (degrees) = 4.3119

For angular frequency w = 31.4159 rad/s:
Transfer Function |H(jw)| : (0.9865,0.0607)
Magnitude (dB): -0.1018
Debugging: Phase (radians) = 0.0614
Debugging: Phase (degrees) = 3.5182

For angular frequency w = 37.6991 rad/s:
Transfer Function |H(jw)| : (0.9905,0.0513)
Magnitude (dB): -0.0715
Debugging: Phase (radians) = 0.0517
Debugging: Phase (degrees) = 2.9639
