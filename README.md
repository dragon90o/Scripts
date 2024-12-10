
<h1 align="center"> 📊 General Transfer Function Analysis in C++</h1>
<h2 align="center">## 🚀 Project Overview</h2>
<p align="center">

Prueva2.cpp is a full C++ implementation that allows you to compute poles, zeros, magnitude in dB, and phase of a general transfer function. The user can also specify a custom sigma value (σ) to perform the desired adjustment to the transfer function's analysis.

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

# Example 
 transfer function: s^2 + 3s + 2/s^2 + 5s + 6 - Min frequency 3- Max frequency 4 - sigma value 0 

## Introduce the numerator (e.g., s^2 + 3s + 2):
s^2 + 3s + 2
## Introduce the denominator (e.g., s^2 + 5s + 6):
s^2 + 5s + 6
## Introduce the minimum frequency (freqMin):
3
## Introduce the maximum frequency (freqMax):
4
## Introduce the sigma value (or) if it is different from 0.0, otherwise enter 0:
6
</p>
--------- 

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

## DEBUGING METODO FINDROOTS

Debugging: Degree of polynomial: 2

Debugging: Coefficients: 1 3 2

Assigning coefficient[2] = 2 to companion matrix.

Assigning coefficient[1] = 3 to companion matrix.

Debugging: Companion matrix:

0 -2 
1 -3

Root before adjustment: (-1,0), after adjustment with sigma: (5,0)

Root before adjustment: (-2,0), after adjustment with sigma: (4,0)

Debugging: Roots (final, adjusted or original):

(5,0)
(4,0)

## DEBUGING METODO FINDROOTS

Debugging: Degree of polynomial: 2

Debugging: Coefficients: 1 5 6

Assigning coefficient[2] = 6 to companion matrix.

Assigning coefficient[1] = 5 to companion matrix.

Debugging: Companion matrix:

 0 -6
 1 -5

Root before adjustment: (-2,0), after adjustment with sigma: (4,0)

Root before adjustment: (-3,0), after adjustment with sigma: (3,0)

Debugging: Roots (final, adjusted or original):

(4,0)
(3,0)

Zeros (Roots of the Numerator): (5,0) (4,0)

Poles (Roots of the Denominator): (4,0) (3,0)

## DEBUGING FREQUENCIES

Debugging: Frequencies method

Debugging: Frequency range - Min: 3, Max: 4

Frequency: 3, Angular Frequency: 18.8496rad/s.

Frequency: 4, Angular Frequency: 25.1327rad/s.

Debugging: Frequency vector size: 2

Debugging: Angular Frequency vector size: 2

For angular frequency w = 18.8496 rad/s:

Debugging: s = sigma + jw = (6.0000,18.8496)

Debugging: Magnitude of s = |s| = 19.7814

Debugging: Evaluating function: s^2 + 3s + 2

Debugging: Angular frequency (w): 18.8496

Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-319.3058,226.1947), Accumulated Result: (-319.3058,226.1947)

Debugging: Term: 3s, Coefficient: 3.0000, Power: 1, Term Value: (18.0000,56.5487), Accumulated Result: (-301.3058,282.7433)

Debugging: Term: 2, Coefficient: 2.0000, Power: 0, Term Value: (2.0000,0.0000), Accumulated Result: (-299.3058,282.7433)

Debugging: Final translated function result: (-299.3058,282.7433)

Debugging: Evaluating function: s^2 + 5s + 6

Debugging: Angular frequency (w): 18.8496

Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-319.3058,226.1947), Accumulated Result: (-319.3058,226.1947)

Debugging: Term: 5s, Coefficient: 5.0000, Power: 1, Term Value: (30.0000,94.2478), Accumulated Result: (-289.3058,320.4425)

Debugging: Term: 6, Coefficient: 6.0000, Power: 0, Term Value: (6.0000,0.0000), Accumulated Result: (-283.3058,320.4425)

Debugging: Final translated function result: (-283.3058,320.4425)
Numerator: (-299.3058,282.7433)

Denominator: (-283.3058,320.4425)

Transfer Function: (0.9587,0.0864)

## Magnitude (dB): -0.3308

DEBUGING PHASE

Debugging: Transfer Function = (0.9587,0.0864)

Debugging: Real = 0.9587

Debugging: Imaginary = 0.0864

Debugging: Phase (radians) = 0.0899

## Debugging: Phase (degrees) = 5.1498

For angular frequency w = 25.1327 rad/s:

Debugging: s = sigma + jw = (6.0000,25.1327)

Debugging: Magnitude of s = |s| = 25.8390

Debugging: Evaluating function: s^2 + 3s + 2

Debugging: Angular frequency (w): 25.1327

Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-595.6547,301.5929), Accumulated Result: (-595.6547,301.5929)

Debugging: Term: 3s, Coefficient: 3.0000, Power: 1, Term Value: (18.0000,75.3982), Accumulated Result: (-577.6547,376.9911)

Debugging: Term: 2, Coefficient: 2.0000, Power: 0, Term Value: (2.0000,0.0000), Accumulated Result: (-575.6547,376.9911)

Debugging: Final translated function result: (-575.6547,376.9911)

Debugging: Evaluating function: s^2 + 5s + 6

Debugging: Angular frequency (w): 25.1327

Debugging: Term: s^2, Coefficient: 1.0000, Power: 2, Term Value: (-595.6547,301.5929), Accumulated Result: (-595.6547,301.5929)

Debugging: Term: 5s, Coefficient: 5.0000, Power: 1, Term Value: (30.0000,125.6637), Accumulated Result: (-565.6547,427.2566)

Debugging: Term: 6, Coefficient: 6.0000, Power: 0, Term Value: (6.0000,0.0000), Accumulated Result: (-559.6547,427.2566)

Debugging: Final translated function result: (-559.6547,427.2566)
Numerator: (-575.6547,376.9911)

Denominator: (-559.6547,427.2566)

Transfer Function: (0.9747,0.0705)

## Magnitude (dB): -0.1995

DEBUGING PHASE

Debugging: Transfer Function = (0.9747,0.0705)

Debugging: Real = 0.9747

Debugging: Imaginary = 0.0705

Debugging: Phase (radians) = 0.0722

## Debugging: Phase (degrees) = 4.1387
</p>
