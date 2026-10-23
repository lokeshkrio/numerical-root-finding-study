
# Comparative Study of Numerical Root-Finding Algorithms

## Overview

A Python-based experimental comparison of four numerical
methods for solving nonlinear equations.

## Algorithms

- Bisection
- Regula Falsi
- Newton-Raphson
- Secant
- Müller (TO BE ADDED)

## Project Structure

- `src/algorithms.py`: algorithms
- `tests/`: automated test cases
- `experiments/run_experiments.py`: benchmark runner
- `experiments/plot_results.py`: visualization
- `results/`: CSV measurements and generated plots

## Requirements

Python 3.10 or newer is recommended.

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run Tests

```bash
python -m pytest -v
```

## Run Experiments

```bash
python -m experiments.run_experiments
```

## Generate Plots

```bash
python -m experiments.plot_results
```

## Evaluation Metrics

- Iteration count
- Absolute residual
- Absolute error relative to the selected reference root
- Median execution time

## Limitations

Results depend on the selected equations, initial guesses,
tolerance, stopping criteria, and computing environment.

## Results

See `results/benchmark_results.csv` and `results/plots/`.


## Future Extension: Müller’s Method
Müller’s method use quadratic interpolation and can produce complex-valued approximations which requires additional consideration of complex arithmetic and convergence criteria.
