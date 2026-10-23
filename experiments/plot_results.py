"""Create plots for the root-finding algorithm study."""

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CSV_FILE = PROJECT_ROOT / "results" / "benchmark_results.csv"
PLOT_DIR = PROJECT_ROOT / "results" / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(CSV_FILE)
successful = df[df["status"] == "completed"]

if successful.empty:
    raise ValueError("No completed results found in the CSV file.")

ALGORITHM_ORDER = ["Bisection", "Regula Falsi", "Newton-Raphson", "Secant"]


def save_bar_plot(
    data: pd.DataFrame | pd.Series,
    filename: str,
    title: str,
    ylabel: str,
    xlabel: str = "Test function",
    figsize: tuple[int, int] = (13, 6),
    rotation: int = 35,
    logy: bool = False,
    legend_title: str | None = None,
) -> None:
    """Helper to generate and save formatted bar plots."""
    ax = data.plot(kind="bar", figsize=figsize, logy=logy)
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    if legend_title:
        ax.legend(title=legend_title)
    plt.xticks(rotation=rotation, ha="right")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / filename, dpi=200)
    plt.close()


# Function-level comparison plots (Iterations, Residuals, Root Errors)
for metric, filename, title, ylabel, logy in [
    (
        "iterations",
        "iterations_by_function.png",
        "Iteration Count by Function and Algorithm",
        "Number of iterations",
        False,
    ),
    (
        "residual",
        "residuals_by_function.png",
        "Absolute Residual by Function and Algorithm",
        "Absolute residual (log scale)",
        True,
    ),
    (
        "root_error",
        "root_errors_by_function.png",
        "Absolute Root Error by Function and Algorithm",
        "Absolute root error (log scale)",
        True,
    ),
]:
    data = successful.pivot(
        index="function", columns="algorithm", values=metric
    ).reindex(columns=ALGORITHM_ORDER)
    if logy:
        data = data.clip(lower=1e-16)
    save_bar_plot(data, filename, title, ylabel, logy=logy, legend_title="Algorithm")

# Aggregate comparison plots (Mean iterations, Median execution time)
algo_group = successful.groupby("algorithm")

save_bar_plot(
    algo_group["iterations"].mean().reindex(ALGORITHM_ORDER),
    "average_iterations.png",
    "Average Iterations by Algorithm",
    "Mean iterations",
    xlabel="Algorithm",
    figsize=(8, 5),
    rotation=20,
)

save_bar_plot(
    algo_group["median_time_seconds"].median().reindex(ALGORITHM_ORDER) * 1_000_000,
    "median_execution_time.png",
    "Median Execution Time by Algorithm",
    "Median time (microseconds)",
    xlabel="Algorithm",
    figsize=(8, 5),
    rotation=20,
)

# Print numerical summary
summary = algo_group.agg(
    mean_iterations=("iterations", "mean"),
    median_time_seconds=("median_time_seconds", "median"),
    median_residual=("residual", "median"),
    median_root_error=("root_error", "median"),
).reindex(ALGORITHM_ORDER)

print("\nExperimental summary:\n", summary.to_string())
print(f"\nPlots saved in: {PLOT_DIR}")
