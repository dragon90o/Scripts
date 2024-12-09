
# 📊 General Transfer Function Analysis in C++

A C++ implementation that allows you to compute the poles, zeros, magnitude in dB, and phase of a general transfer function. The code analyzes the transfer function, calculates its poles and zeros, evaluates its magnitude response (in dB), and computes the phase across specified angular frequencies.

## 🚀 Project Overview
This project provides a modular approach to analyze transfer functions. Using mathematical techniques, it determines:

Poles and Zeros: Roots of the numerator and denominator polynomials.
Magnitude Response: Measured in decibels (dB) over a range of angular frequencies.
Phase Response: The phase angle (in degrees) across a range of angular frequencies.
This functionality is useful for analyzing control systems, signal processing applications, or frequency response characteristics.

## 🛠️ Features

✅ Calculate poles and zeros of the transfer function from polynomial coefficients.
✅ Evaluate magnitude in dB at specific angular frequencies.
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
ExprTk Library:

Purpose: A symbolic mathematical engine for evaluating expressions dynamically.
Installation Path: A single header file, exprtk.hpp, needs to be placed in the correct path.
Download from: ExprTk GitHub Repository

# Example 
## transfer function: s^2 + 3s + 2/s^2 + 5s + 6
## Min frequency 3
## Max frequency 4
## sigma value 0
Introduce the numerator (e.g., s^2 + 3s + 2):
s^2 + 3s + 2
Introduce the denominator (e.g., s^2 + 5s + 6):
s^2 + 5s + 6
Introduce the minimum frequency (freqMin):
3
Introduce the maximum frequency (freqMax):
4
Introduce the sigma value (or) if it is different from 0.0, otherwise enter 0:
0
Debugging: Found 10 matches.
Debugging: term: s^2, coefficient: 1, power: 2, maxPower: 2
Debugging: term: 3s, coefficient: 3, power: 1, maxPower: 2
Debugging: term: 2, coefficient: 2, power: 0, maxPower: 2
Debugging: Coefficients vector: 1 3 2
Debugging: Found 10 matches.
Debugging: term: s^2, coefficient: 1, power: 2, maxPower: 2
Debugging: term: 5s, coefficient: 5, power: 1, maxPower: 2
Debugging: term: 6, coefficient: 6, power: 0, maxPower: 2
Debugging: Coefficients vector: 1 5 6
DEBUGING METODO FINDROOTS
Debugging: Degree of polynomial: 2
Debugging: Coefficients: 1 3 2
Assigning coefficient[2] = 2 to companion matrix.
Assigning coefficient[1] = 3 to companion matrix.
Debugging: Companion matrix:
 0 -2
 1 -3
Debugging: Roots (eigenvalues):
(-1,0)
(-2,0)
DEBUGING METODO FINDROOTS
Debugging: Degree of polynomial: 2
Debugging: Coefficients: 1 5 6
Assigning coefficient[2] = 6 to companion matrix.
Assigning coefficient[1] = 5 to companion matrix.
Debugging: Companion matrix:
 0 -6
 1 -5
Debugging: Roots (eigenvalues):
(-2,0)
(-3,0)
Zeros (Roots of the Numerator): (-1,0) (-2,0)
Poles (Roots of the Denominator): (-2,0) (-3,0)
DEBUGING FREQUENCIES
Debugging: Frequencies method
Debugging: Frequency range - Min: 3, Max: 4
Frequency: 3, Angular Frequency: 18.8496rad/s.
Frequency: 4, Angular Frequency: 25.1327rad/s.
Debugging: Frequency vector size: 2
Debugging: Angular Frequency vector size: 2

For angular frequency w = 18.8496 rad/s:
Debugging: Evaluating function: s^2 + 3s + 2
Debugging: Angular frequency (w): 18.8496
Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-355.3058,0.0000), Accumulated Result: (-355.3058,0.0000)
Debugging: Term: 3s, Coefficient: 3.0000, Power: 1, Term Value: (0.0000,56.5487), Accumulated Result: (-355.3058,56.5487)
Debugging: Term: 2, Coefficient: 2.0000, Power: 0, Term Value: (2.0000,0.0000), Accumulated Result: (-353.3058,56.5487)
Debugging: Final translated function result: (-353.3058,56.5487)
Debugging: Evaluating function: s^2 + 5s + 6
Debugging: Angular frequency (w): 18.8496
Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-355.3058,0.0000), Accumulated Result: (-355.3058,0.0000)
Debugging: Term: 5s, Coefficient: 5.0000, Power: 1, Term Value: (0.0000,94.2478), Accumulated Result: (-355.3058,94.2478)
Debugging: Term: 6, Coefficient: 6.0000, Power: 0, Term Value: (6.0000,0.0000), Accumulated Result: (-349.3058,94.2478)
Debugging: Final translated function result: (-349.3058,94.2478)
Numerator: (-353.3058,56.5487)
Denominator: (-349.3058,94.2478)
Transfer Function: (0.9835,0.1035)
Magnitude (dB): -0.0964
DEBUGING PHASE
Debugging: Transfer Function = (0.9835,0.1035)
Debugging: Real = 0.9835
Debugging: Imaginary = 0.1035
Debugging: Phase (radians) = 0.1048
Debugging: Phase (degrees) = 6.0063

For angular frequency w = 25.1327 rad/s:
Debugging: Evaluating function: s^2 + 3s + 2
Debugging: Angular frequency (w): 25.1327
Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-631.6547,0.0000), Accumulated Result: (-631.6547,0.0000)
Debugging: Term: 3s, Coefficient: 3.0000, Power: 1, Term Value: (0.0000,75.3982), Accumulated Result: (-631.6547,75.3982)
Debugging: Term: 2, Coefficient: 2.0000, Power: 0, Term Value: (2.0000,0.0000), Accumulated Result: (-629.6547,75.3982)
Debugging: Final translated function result: (-629.6547,75.3982)
Debugging: Evaluating function: s^2 + 5s + 6
Debugging: Angular frequency (w): 25.1327
Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-631.6547,0.0000), Accumulated Result: (-631.6547,0.0000)
Debugging: Term: 5s, Coefficient: 5.0000, Power: 1, Term Value: (0.0000,125.6637), Accumulated Result: (-631.6547,125.6637)
Debugging: Term: 6, Coefficient: 6.0000, Power: 0, Term Value: (6.0000,0.0000), Accumulated Result: (-625.6547,125.6637)
Debugging: Final translated function result: (-625.6547,125.6637)
Numerator: (-629.6547,75.3982)
Denominator: (-625.6547,125.6637)
Transfer Function: (0.9906,0.0785)
Magnitude (dB): -0.0546
DEBUGING PHASE
Debugging: Transfer Function = (0.9906,0.0785)
Debugging: Real = 0.9906
Debugging: Imaginary = 0.0785
Debugging: Phase (radians) = 0.0790
Debugging: Phase (degrees) = 4.5284
