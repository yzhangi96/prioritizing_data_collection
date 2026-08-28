import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from scipy.stats import pearsonr


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def equation(intercept, slope):
    sign = "+" if slope >= 0 else "-"
    return f"y = {intercept:.2f} {sign} {abs(slope):.2f}x"


def plot_simulation(data, output):
    columns = ["AvgDUC", "AvgNegKL", "AvgNegScoreX"]
    labels = [
        "Avg Estimated DUC",
        "Avg Estimated Negative KL",
        "Avg Estimated Negative Domain Classifier Score",
    ]
    y = data["NegRelImp"].to_numpy(float)
    fig, axes = plt.subplots(
        1, 3, figsize=(18, 6), sharey=True, facecolor="white"
    )

    for ax, column, label in zip(axes, columns, labels):
        x = data[column].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        line_x = np.linspace(x.min(), x.max(), 100)
        ax.scatter(x, y, s=100, color="black", zorder=2)
        ax.plot(
            line_x, intercept + slope * line_x,
            color="lightgrey", linestyle="--", linewidth=2,
        )
        ax.text(
            0.95, 0.95, equation(intercept, slope),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=12, weight="bold",
        )
        ax.text(
            0.95, 0.87, f"Correlation = {pearsonr(x, y)[0]:.2f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=12, weight="bold",
        )
        ax.set_xlabel(label, fontsize=16)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        if column == "AvgNegScoreX":
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.5f"))
        ax.tick_params(axis="both", labelsize=14)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.2)

    axes[0].set_ylabel("Average Negative Relative Excess MSE", fontsize=16)
    axes[0].set_ylim(-1, 0)
    fig.suptitle(
        "Predicting model performance without outcome data",
        fontsize=24, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300, facecolor="white")
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path, default=HERE / "simulation_plot_data.csv"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "figures" / "sim.png"
    )
    args = parser.parse_args()
    data = pd.read_csv(args.input)
    plot_simulation(data, args.output)
    plt.close("all")


if __name__ == "__main__":
    main()
