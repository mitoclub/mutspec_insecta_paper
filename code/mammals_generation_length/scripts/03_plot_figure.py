#!/usr/bin/env python3
"""Generate the mammalian Figure 3A scatter plot and PGLS fitted line."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PRIMARY_ANALYSIS = "PRIMARY_all_species_log10_generation_length"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    data = pd.read_csv(args.dataset)
    results = pd.read_csv(args.results)
    primary = results.loc[
        results["analysis"].eq(PRIMARY_ANALYSIS)
    ].iloc[0]

    x = data["generation_length_d"].astype(float).to_numpy()
    y = data["pro_over_phe_plus_leutt"].astype(float).to_numpy()

    x_line = np.geomspace(x.min(), x.max(), 400)
    y_line = (
        float(primary["intercept"])
        + float(primary["slope"]) * np.log10(x_line)
    )

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.scatter(x, y, s=11, alpha=0.18)
    ax.plot(x_line, y_line, linewidth=2.3)

    ax.set_xscale("log")
    ax.set_xlabel("Generation length, days (log scale)")
    ax.set_ylabel("CYTB Pro / (Phe + LeuTT)")
    ax.set_title("Mammals")

    annotation = (
        f"PGLS, n = {int(primary['n'])}\n"
        f"β = {primary['slope']:.4f}\n"
        f"95% CI [{primary['slope_ci_low']:.4f}, "
        f"{primary['slope_ci_high']:.4f}]\n"
        f"p = {primary['slope_p']:.4g}; "
        f"Pagel λ = {primary['lambda']:.3f}"
    )
    ax.text(
        0.03,
        0.96,
        annotation,
        transform=ax.transAxes,
        va="top",
        fontsize=11,
    )

    fig.tight_layout()

    args.pdf.parent.mkdir(parents=True, exist_ok=True)
    args.png.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(args.pdf)
    fig.savefig(args.png, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
