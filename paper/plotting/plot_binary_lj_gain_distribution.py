#!/usr/bin/env python3
"""Plot type-resolved Mi-RBL gain histograms for the binary LJ mixture."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
ALPHA_TAGS = ((0.25, "0p25"), (0.31, "0p31"), (0.37, "0p37"))
ALPHA_COLORS = {0.25: "#1f77b4", 0.31: "#ff7f0e", 0.37: "#2ca02c"}


def plot_species_multiple(output_dir, curves):
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.size": 10})
    fig, ax = plt.subplots(figsize=(4.6, 3.5))
    legend_handles = []
    for alpha, groups in curves:
        color = ALPHA_COLORS.get(alpha)
        for group, linestyle, type_label in (("type_1", "-", "Type 1"),
                                             ("type_2", "--", "Type 2")):
            edges, density = groups[group]
            centers = (edges[:-1] + edges[1:]) / 2
            ax.plot(centers, density, linewidth=1.3, color=color, linestyle=linestyle)
            legend_handles.append(Line2D([], [], color=color, linewidth=1.7,
                                         linestyle=linestyle,
                                         label=rf"$\alpha={alpha:g}$, {type_label}"))
    ax.axvline(1, color="0.35", linestyle="--", linewidth=1)
    ax.set_xlabel(r"Particle gain $g_i$", fontsize=12)
    ax.set_ylabel(r"Probability density per $g_i$", fontsize=11)
    ax.set_xlim(
        min(edges[0] for _, groups in curves for edges, _ in groups.values()),
        max(edges[-1] for _, groups in curves for edges, _ in groups.values()),
    )
    ax.set_axisbelow(True)
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color=(0.85, 0.85, 0.85), alpha=1.0)
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False,
               fontsize=8.5, handlelength=1.6, handletextpad=0.4, columnspacing=1.3)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(output_dir / "binary_lj_gain_distribution.png", dpi=400)
    plt.close(fig)


def read_curves():
    curves = []
    for alpha, tag in ALPHA_TAGS:
        path = DATA_DIR / f"binlj_gain_{tag}.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        groups = {}
        for group in ("type_1", "type_2"):
            group_rows = [row for row in rows if row["group"] == group]
            edges = np.array([float(group_rows[0]["gain_left"])] +
                             [float(row["gain_right"]) for row in group_rows])
            density = np.array([float(row["density_per_gain"]) for row in group_rows])
            groups[group] = (edges, density)
        curves.append((alpha, groups))
    return curves


def main():
    plot_species_multiple(SCRIPT_DIR.parent / "figures", read_curves())


if __name__ == "__main__":
    main()
