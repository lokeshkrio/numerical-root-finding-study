"""Root-finding algorithms package."""

from src.algorithms import (
    bisection,
    muller,
    newton_raphson,
    regula_falsi,
    secant,
)

__all__ = [
    "bisection",
    "regula_falsi",
    "newton_raphson",
    "secant",
    "muller",
]
