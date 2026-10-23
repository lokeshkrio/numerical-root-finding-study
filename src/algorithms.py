"""Root-finding algorithms."""

from collections.abc import Callable


def _validate_params(tol: float, max_iter: int) -> None:
    """Validates tolerance and iteration limits."""
    if tol <= 0:
        raise ValueError("Tolerance must be positive.")
    if not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("Maximum iterations must be an integer >= 1.")


def bisection(
    f: Callable[[float], float],
    a: float,
    b: float,
    tol: float = 1e-6,
    max_iter: int = 1000,
) -> tuple[float, int]:
    """Find a root in a sign-changing interval using bisection.

    Convergence criteria:
        Half-interval width <= tol and residual <= tol,
        or an exact root is found.

    Raises:
        ValueError: If the interval does not bracket a root.
        RuntimeError: If the method does not converge within max_iter.
    """
    _validate_params(tol, max_iter)

    fa, fb = f(a), f(b)

    if fa == 0:
        return a, 0
    if fb == 0:
        return b, 0

    if fa * fb > 0:
        raise ValueError("The interval must bracket a root.")

    for iteration in range(1, max_iter + 1):
        c = a + (b - a) / 2
        fc = f(c)

        if fc == 0:
            return c, iteration

        if abs(b - a) / 2 <= tol and abs(fc) <= tol:
            return c, iteration

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    raise RuntimeError("Bisection did not converge.")


def regula_falsi(
    f: Callable[[float], float],
    a: float,
    b: float,
    tol: float = 1e-6,
    max_iter: int = 1000,
) -> tuple[float, int]:
    """Find a root in a sign-changing interval using the regula falsi method.

    Convergence criteria:
        Residual <= tol, or an exact root is found.

    Raises:
        ValueError: If the interval does not bracket a root.
        RuntimeError: If the method does not converge within max_iter.
    """
    _validate_params(tol, max_iter)
    fa, fb = f(a), f(b)

    if fa == 0:
        return a, 0
    if fb == 0:
        return b, 0
    if fa * fb > 0:
        raise ValueError("The interval must bracket a root.")

    for iteration in range(1, max_iter + 1):
        denominator = fb - fa
        if denominator == 0:
            raise ZeroDivisionError("Cannot calculate the next estimate.")

        c = (a * fb - b * fa) / denominator
        fc = f(c)

        if fc == 0 or abs(fc) <= tol:
            return c, iteration

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    raise RuntimeError("Regula Falsi did not converge.")


def newton_raphson(
    f: Callable[[float], float],
    df: Callable[[float], float],
    x0: float,
    tol: float = 1e-6,
    max_iter: int = 1000,
) -> tuple[float, int]:
    """Find a root using the Newton-Raphson iterative method.

    Convergence criteria:
        Residual <= tol, or an exact root is found.

    Raises:
        ValueError: If the derivative is zero or the inputs are invalid.
        RuntimeError: If the method does not converge within max_iter.
    """
    _validate_params(tol, max_iter)
    x, fx = x0, f(x0)

    if fx == 0:
        return x, 0

    for iteration in range(1, max_iter + 1):
        dfx = df(x)
        if abs(dfx) < 1e-14:
            raise ZeroDivisionError("Derivative is zero or too small.")

        x_new = x - fx / dfx
        f_new = f(x_new)

        if f_new == 0 or (abs(x_new - x) <= tol and abs(f_new) <= tol):
            return x_new, iteration

        x, fx = x_new, f_new

    raise RuntimeError("Newton-Raphson did not converge.")


def secant(
    f: Callable[[float], float],
    x0: float,
    x1: float,
    tol: float = 1e-6,
    max_iter: int = 1000,
) -> tuple[float, int]:
    """Find a root using the secant method without requiring a derivative.

    Convergence criteria:
        Residual <= tol, or an exact root is found.

    Raises:
        ValueError: If the initial guesses are invalid.
        RuntimeError: If the method does not converge within max_iter.
    """
    _validate_params(tol, max_iter)
    f0, f1 = f(x0), f(x1)

    if f0 == 0:
        return x0, 0
    if f1 == 0:
        return x1, 0

    for iteration in range(1, max_iter + 1):
        denominator = f1 - f0
        if denominator == 0:
            raise ZeroDivisionError("Secant slope cannot be calculated.")

        x2 = x1 - f1 * (x1 - x0) / denominator
        f2 = f(x2)

        if f2 == 0 or (abs(x2 - x1) <= tol and abs(f2) <= tol):
            return x2, iteration

        x0, f0, x1, f1 = x1, f1, x2, f2

    raise RuntimeError("Secant method did not converge.")


def muller(
    f: Callable[[complex | float], complex | float],
    x0: complex | float,
    x1: complex | float,
    x2: complex | float,
    tol: float = 1e-8,
    max_iter: int = 1000,
) -> tuple[complex | float, int]:
    """Find a root using quadratic interpolation through three points.

    Convergence criteria:
        Residual <= tol, or an exact root is found.

    Raises:
        ValueError: If the initial guesses or interpolation are invalid.
        RuntimeError: If the method does not converge within max_iter.
    """
    _validate_params(tol, max_iter)
    f0, f1, f2 = f(x0), f(x1), f(x2)

    for x, fx in ((x0, f0), (x1, f1), (x2, f2)):
        if fx == 0:
            return x, 0

    for iteration in range(1, max_iter + 1):
        h0, h1 = x1 - x0, x2 - x1
        if 0 in (h0, h1, h0 + h1):
            raise ValueError("Invalid or repeated guesses.")

        delta0, delta1 = (f1 - f0) / h0, (f2 - f1) / h1
        d = (delta1 - delta0) / (h1 + h0)
        b = delta1 + h1 * d
        c = f2
        discriminant = b * b - 4 * d * c

        # Complex roots are possible in Muller's method.
        sqrt_disc = complex(discriminant) ** 0.5
        denominator = max(b + sqrt_disc, b - sqrt_disc, key=abs)

        if abs(denominator) < 1e-15:
            raise ZeroDivisionError("Unstable quadratic step.")

        dx = -2 * c / denominator
        x3 = x2 + dx
        f3 = f(x3)

        if abs(dx) <= tol * (1 + abs(x3)) and abs(f3) <= tol:
            return x3, iteration

        x0, f0, x1, f1, x2, f2 = x1, f1, x2, f2, x3, f3

    raise RuntimeError("Muller's method did not converge.")
