# Chebyshev-convergence-matlab

# Chebyshev Convergence in MATLAB

This project focuses on the decay of Chebyshev interpolation coefficients for three test functions sampled at Chebyshev extreme points on `[-1,1]`. It also studies how the oscillation parameter `omega` affects the length of the preconvergent phase, and how the observed coefficient decay compares with theoretical predictions for analytic and finitely smooth functions.

## Functions

Three functions are involved in this project:

- `f(x) = sin(omega*x + 1)`
- `g(x) = sin(omega*x + 1) / (4(x - 0.1)^2 + 1)`
- `h(x) = sin(omega*x + 1) * |(x - 0.1)^3|`

Each was chosen to exemplify one specific scenario:

- `f(x)`: Analytic and entire function.
- `g(x)`: Analytic on `[-1,1]`, but coefficient decay affected by complex singularities.
- `h(x)`: More slowly decay due to absolute value term.

## Aims

Two main aims:

1. Study how the parameter `omega` affects the length of the preconvergent phase in the Chebyshev coefficients.
2. Compare the computed coefficient decay with theoretical bounding condition for analytic and nonanalytic functions.

## Files

- `Chebyshev_Convergence_matlab.m`: Main MATLAB file to run the convergence test and bounding condition comparisons.

- `discrete_cosine_transform.m`: Helper code for computing Chebyshev coefficients through a DCT-I / FFT method.

## Process

### Part I: Varying `Omega`
This project evaluates the three functions when

- `omega = 100`
- `omega = 60`
- `omega = 20`

and plots the absolute values of their Chebyshev coefficients. This part serves for observation of preconvergent phase change as `omega` changes.

### Part II: Bounding Condition Comparisons
Next, this project fixes `omega = 20` and compares the computed coefficients with the bounding curves motivated by two standard types of convergence behavior:

- exponential decay for analytic functions
- power decay for finitely smooth functions

## How to run

Open MATLAB in this folder and run

```matlab
Chebyshev_Convergence_matlab
