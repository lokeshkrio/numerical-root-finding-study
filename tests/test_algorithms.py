"""Unit tests for root-finding algorithms."""

import math
import pytest

from src.algorithms import (
    bisection,
    muller,
    newton_raphson,
    regula_falsi,
    secant,
)

TOL = 1e-6

# Test functions and derivatives
f1, df1 = (lambda x: x**2 - 2), (lambda x: 2 * x)
f2, df2 = (lambda x: x**3 - x - 2), (lambda x: 3 * x**2 - 1)
f3, df3 = (lambda x: x**2 - 9), (lambda x: 2 * x)

# Test cases: (f, a, b, expected_root)
INTERVAL_CASES = [
    (f1, 1.0, 2.0, math.sqrt(2)),
    (f2, 1.0, 2.0, 1.5213797068),
    (f3, 2.0, 4.0, 3.0),
]

NEWTON_CASES = [
    (f1, df1, 1.0, math.sqrt(2)),
    (f2, df2, 1.5, 1.5213797068),
    (f3, df3, 2.0, 3.0),
]

MULLER_CASES = [
    (f1, 0.0, 1.0, 2.0, math.sqrt(2)),
    (f2, 1.0, 1.5, 2.0, 1.5213797068),
    (f3, 2.0, 2.5, 4.0, 3.0),
]


@pytest.mark.parametrize("f,a,b,expected", INTERVAL_CASES)
def test_bisection(f, a, b, expected):
    root, iterations = bisection(f, a, b, tol=TOL)
    assert abs(root - expected) <= TOL
    assert iterations >= 0


@pytest.mark.parametrize("f,a,b,expected", INTERVAL_CASES)
def test_regula_falsi(f, a, b, expected):
    root, iterations = regula_falsi(f, a, b, tol=TOL)
    assert abs(root - expected) < TOL
    assert abs(f(root)) < TOL
    assert iterations >= 0


@pytest.mark.parametrize("f,df,x0,expected", NEWTON_CASES)
def test_newton_raphson(f, df, x0, expected):
    root, iterations = newton_raphson(f, df, x0, tol=TOL)
    assert abs(root - expected) < TOL
    assert abs(f(root)) < TOL
    assert iterations >= 0


@pytest.mark.parametrize("f,x0,x1,expected", INTERVAL_CASES)
def test_secant(f, x0, x1, expected):
    root, iterations = secant(f, x0, x1, tol=TOL)
    assert abs(root - expected) < TOL
    assert abs(f(root)) < TOL
    assert iterations >= 0


@pytest.mark.parametrize("f,x0,x1,x2,expected", MULLER_CASES)
def test_muller(f, x0, x1, x2, expected):
    root, iterations = muller(f, x0, x1, x2, tol=1e-8)
    assert abs(root - expected) < 1e-6
    assert abs(f(root)) < 1e-6
    assert iterations >= 0


# Edge cases and validation tests


@pytest.mark.parametrize(
    "method,args",
    [
        (bisection, (f1, 1.0, 2.0)),
        (regula_falsi, (f1, 1.0, 2.0)),
        (newton_raphson, (f1, df1, 1.0)),
        (secant, (f1, 1.0, 2.0)),
        (muller, (f1, 0.0, 1.0, 2.0)),
    ],
)
def test_invalid_parameters(method, args):
    with pytest.raises(ValueError, match="Tolerance must be positive"):
        method(*args, tol=0)
    with pytest.raises(ValueError, match="Maximum iterations must be an integer >= 1"):
        method(*args, max_iter=0)


def test_bracketing_validation():
    with pytest.raises(ValueError, match="must bracket a root"):
        bisection(f1, 2.0, 3.0)
    with pytest.raises(ValueError, match="must bracket a root"):
        regula_falsi(f1, 2.0, 3.0)


def test_exact_root_at_start():
    # Exact roots should return immediately with iteration count of 0
    assert bisection(f3, 3.0, 4.0) == (3.0, 0)
    assert regula_falsi(f3, 2.0, 3.0) == (3.0, 0)
    assert newton_raphson(f3, df3, 3.0) == (3.0, 0)
    assert secant(f3, 3.0, 4.0) == (3.0, 0)
    assert muller(f3, 3.0, 2.0, 4.0) == (3.0, 0)


def test_zero_division_errors():
    with pytest.raises(ZeroDivisionError):
        newton_raphson(lambda x: x**2 - 1, lambda x: 0.0, 0.0)
    with pytest.raises(ZeroDivisionError):
        secant(lambda x: 1.0, 0.0, 1.0)
