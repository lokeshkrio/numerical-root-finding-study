"""Run and compare root-finding algorithms."""

from collections.abc import Callable
import csv
import math
from pathlib import Path
import statistics
import sys
import time
from typing import NamedTuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.algorithms import bisection, newton_raphson, regula_falsi, secant

# --------------------------------------------------
# Configuration
# --------------------------------------------------

TOL = 1e-6
MAX_ITER = 1000
REPEATS = 20


class Experiment(NamedTuple):
    function: str
    f: Callable[[float], float]
    df: Callable[[float], float]
    root: float
    guesses: tuple[float, float]
    newton_guess: float


EXPERIMENTS = [
    Experiment(
        "x^2 - 2", lambda x: x**2 - 2, lambda x: 2 * x, math.sqrt(2), (1.0, 2.0), 1.0
    ),
    Experiment(
        "x^3 - x - 2",
        lambda x: x**3 - x - 2,
        lambda x: 3 * x**2 - 1,
        1.5213797068,
        (1.0, 2.0),
        1.5,
    ),
    Experiment(
        "cos(x) - x",
        lambda x: math.cos(x) - x,
        lambda x: -math.sin(x) - 1,
        0.7390851332,
        (0.0, 1.0),
        0.5,
    ),
    Experiment(
        "(x-1)^2*(x+2)",
        lambda x: (x - 1) ** 2 * (x + 2),
        lambda x: 3 * (x - 1) * (x + 1),
        -2.0,
        (-3.0, 0.0),
        -2.5,
    ),
    Experiment(
        "x^5 - x - 1",
        lambda x: x**5 - x - 1,
        lambda x: 5 * x**4 - 1,
        1.1673039783,
        (1.0, 2.0),
        1.2,
    ),
    Experiment(
        "exp(-x) - x",
        lambda x: math.exp(-x) - x,
        lambda x: -math.exp(-x) - 1,
        0.5671432904,
        (0.0, 1.0),
        0.5,
    ),
    Experiment(
        "sin(x) - 0.5",
        lambda x: math.sin(x) - 0.5,
        lambda x: math.cos(x),
        math.pi / 6,
        (0.0, 1.0),
        0.5,
    ),
    Experiment(
        "ln(x) + x - 2",
        lambda x: math.log(x) + x - 2,
        lambda x: 1 / x + 1,
        1.5571455990,
        (1.0, 2.0),
        1.5,
    ),
    Experiment(
        "x^3 - 2*x + 2",
        lambda x: x**3 - 2 * x + 2,
        lambda x: 3 * x**2 - 2,
        -1.7692923542,
        (-2.0, -1.0),
        -1.5,
    ),
    Experiment(
        "(x^2 - 3)/(x + 2)",
        lambda x: (x**2 - 3) / (x + 2),
        lambda x: (x**2 + 4 * x + 3) / (x + 2) ** 2,
        math.sqrt(3),
        (1.0, 2.0),
        1.5,
    ),
    Experiment(
        "x*exp(x) - 1",
        lambda x: x * math.exp(x) - 1,
        lambda x: math.exp(x) * (x + 1),
        0.5671432904,
        (0.0, 1.0),
        0.5,
    ),
]


def run_experiment(f, expected_root, algorithm, args):
    """Run an algorithm repeatedly and evaluate its result."""

    times = []
    root = None
    iterations = None
    status = "failed"
    error_message = ""

    for _ in range(REPEATS):
        start = time.perf_counter()

        try:
            root, iterations = algorithm(f, *args, tol=TOL, max_iter=MAX_ITER)
            times.append(time.perf_counter() - start)

        except Exception as exc:
            error_message = f"{type(exc).__name__}: {exc}"
            break

    residual = float("nan")
    root_error = float("nan")
    residual_ok = False
    target_root_found = False

    if not error_message:
        try:
            if not math.isfinite(root):
                raise ValueError("Algorithm returned a non-finite root")

            residual = abs(f(root))
            root_error = abs(root - expected_root)

            if not (math.isfinite(residual) and math.isfinite(root_error)):
                raise ValueError("Residual or root error is non-finite")

            # These checks are independent.
            residual_ok = residual <= TOL
            target_root_found = root_error <= TOL
            status = "completed"

        except (TypeError, ValueError, OverflowError) as exc:
            status = "invalid_result"
            error_message = f"{type(exc).__name__}: {exc}"

    return {
        "status": status,
        "root": root,
        "iterations": iterations,
        "residual": residual,
        "root_error": root_error,
        "residual_ok": residual_ok,
        "target_root_found": target_root_found,
        "median_time_seconds": (statistics.median(times) if times else float("nan")),
        "error_message": error_message,
    }


def main():
    results = []

    for exp in EXPERIMENTS:
        methods = [
            ("Bisection", bisection, exp.guesses),
            ("Regula Falsi", regula_falsi, exp.guesses),
            ("Newton-Raphson", newton_raphson, (exp.df, exp.newton_guess)),
            ("Secant", secant, exp.guesses),
        ]

        print(f"\nFunction: {exp.function}")
        for name, algorithm, args in methods:
            result = run_experiment(exp.f, exp.root, algorithm, args)
            results.append(
                {
                    "function": exp.function,
                    "algorithm": name,
                    "tolerance": TOL,
                    "repeats": REPEATS,
                    **result,
                }
            )
            print(
                f"{name:16} | {result['status']:14} | "
                f"iterations={str(result['iterations']):>4} | "
                f"residual_ok={result['residual_ok']} | "
                f"target_root={result['target_root_found']}"
            )

    output_dir = PROJECT_ROOT / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "benchmark_results.csv"

    columns = [
        "function",
        "algorithm",
        "tolerance",
        "repeats",
        "status",
        "root",
        "iterations",
        "residual",
        "root_error",
        "residual_ok",
        "target_root_found",
        "median_time_seconds",
        "error_message",
    ]

    with output_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nResults saved to: {output_file.resolve()}")


if __name__ == "__main__":
    main()
