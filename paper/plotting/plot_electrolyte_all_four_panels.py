"""Combine the electrolyte kinetic, gain, and timing results in four panels."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, MultipleLocator
import numpy as np

import electrolyte_kinetic_gamma_powerlaw as gamma_plot
import electrolyte_kinetic_timestep_linear as timestep_plot
import electrolyte_timing_data as timing_data


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_PATH = SCRIPT_DIR.parent / "figures" / "electrolyte_all_four_panels.png"
GAIN_INPUTS = {
    p: DATA_DIR / f"electrolyte_gain_P{p}.csv" for p in (5, 10)
}
KINETIC_CASES = (
    ("mirbe_a0.0_rho0.95", "RBE", 0.0, "#1f77b4", "s"),
    ("mirbe_a0.4_rho0.95", r"Mi-RBE ($\alpha=0.4$)", 0.4, "#2ca02c", "^"),
    ("mirbe_a0.5_rho0.95", r"Mi-RBE ($\alpha=0.5$)", 0.5, "#d62728", "o"),
    ("mirbe_a0.6_rho0.95", r"Mi-RBE ($\alpha=0.6$)", 0.6, "#9467bd", "v"),
)


def read_gain_curves(path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    curves = {}
    for alpha in (0.4, 0.5, 0.6):
        groups = {}
        for group in ("type_1", "type_2"):
            selected = [row for row in rows if float(row["alpha"]) == alpha
                        and row["group"] == group]
            edges = np.array([float(selected[0]["gain_left"])] +
                             [float(row["gain_right"]) for row in selected])
            density = np.array([float(row["density_per_log10_gain"]) for row in selected])
            groups[group] = (edges, density)
        curves[alpha] = groups
    return curves


def plot_gamma(ax, points_by_alpha):
    inset = ax.inset_axes([0.095, 0.585, 0.47, 0.40])
    for _, label, alpha, color, marker in KINETIC_CASES:
        points = points_by_alpha[alpha]
        line_fit = gamma_plot.fit_line(points)
        power_fit = gamma_plot.fit_power_law(points)
        damps = [point.damp for point in points]
        ax.plot(damps, [point.kinetic_error for point in points], color=color,
                marker=marker, markersize=5, markerfacecolor="none",
                linestyle="None", label=label)
        fit_damps = gamma_plot.linspace_between(min(damps), max(damps))
        ax.plot(fit_damps, [line_fit.predict(damp) for damp in fit_damps],
                color=color, linestyle="--", linewidth=1.5)

        gammas = [point.gamma for point in points]
        fit_gammas = gamma_plot.logspace_between(min(gammas), max(gammas))
        slope_label = "RBE" if alpha == 0 else rf"$\alpha={alpha:g}$"
        inset.plot([math.log10(gamma) for gamma in gammas],
                   [math.log10(point.kinetic_error) for point in points],
                   color=color, marker=marker, markersize=3.5,
                   markerfacecolor="none", linestyle="None")
        inset.plot([math.log10(gamma) for gamma in fit_gammas],
                   [math.log10(power_fit.prefactor * gamma**power_fit.exponent)
                    for gamma in fit_gammas], color=color, linestyle="--",
                   linewidth=1.1, label=f"({slope_label}): {power_fit.exponent:.2f}")

    ax.set_title(r"(a) Friction dependence, $P=10$", fontsize=10)
    ax.set_xlabel(r"$1/\gamma$")
    ax.set_ylabel(r"Per-atom kinetic energy error [$\epsilon_{\mathrm{LJ}}$]")
    ax.set_xticks([0, 1, 2, 4, 10, 20])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_ylim(0, 0.24)
    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.legend(frameon=False, fontsize=7, loc="lower right", handlelength=1.8)
    inset.set_xlabel(r"$\log_{10}\gamma$", fontsize=7, labelpad=1)
    inset.set_xlim(-1.35, 0.04)
    inset.set_ylabel(r"$\log_{10}\Delta k$", fontsize=7, labelpad=1)
    inset.tick_params(labelsize=6, length=2, pad=1)
    inset.grid(True, linewidth=0.35, alpha=0.28)
    inset.legend(title="slope", frameon=False, fontsize=5.7,
                 title_fontsize=5.7, loc="lower left",
                 bbox_to_anchor=(0, -0.04), handlelength=1.4)


def plot_timestep(ax, points_by_alpha):
    for _, label, alpha, color, marker in KINETIC_CASES:
        points = points_by_alpha[alpha]
        fit = timestep_plot.fit_line(points)
        steps = [point.timestep for point in points]
        ax.plot(steps, [point.kinetic_error for point in points], color=color,
                marker=marker, markersize=5, markerfacecolor="none",
                linestyle="None", label=f"{label}: {fit.slope:.2f}")
        fit_steps = timestep_plot.linspace_between(min(steps), max(steps))
        ax.plot(fit_steps, [fit.predict(step) for step in fit_steps],
                color=color, linestyle="--", linewidth=1.5)

    ax.set_title(r"(b) Time-step dependence, $P=10$", fontsize=10)
    ax.set_xlabel(r"Timestep $\Delta t$ [$\tau_{\mathrm{LJ}}$]")
    ax.set_ylabel(r"Per-atom kinetic energy error [$\epsilon_{\mathrm{LJ}}$]")
    ax.set_xlim(left=0)
    ax.set_xticks([0, 0.001, 0.005, 0.01])
    ax.set_xticklabels(["0", r"$10^{-3}$", "0.005", r"$10^{-2}$"])
    ax.xaxis.set_minor_locator(MultipleLocator(0.001))
    ax.set_ylim(0, 0.205)
    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.legend(title="slope", frameon=False, fontsize=7, title_fontsize=7,
              loc="upper left", handlelength=1.8)


def plot_gain(ax):
    gains = {p: read_gain_curves(path) for p, path in GAIN_INPUTS.items()}
    colors = {
        5: {0.4: "#1f77b4", 0.5: "#ff7f0e", 0.6: "#2ca02c"},
        10: {0.4: "#d62728", 0.5: "#9467bd", 0.6: "#8c564b"},
    }
    handles = []
    for p in (5, 10):
        for alpha in (0.4, 0.5, 0.6):
            color = colors[p][alpha]
            handles.append(Line2D([], [], color=color, linewidth=1.5,
                                  label=rf"$P={p},\ \alpha={alpha:g}$"))
            for group, linestyle in (("type_1", "-"), ("type_2", "--")):
                edges, density = gains[p][alpha][group]
                centers = (edges[:-1] + edges[1:]) / 2
                ax.plot(centers, density, color=color, linestyle=linestyle,
                        linewidth=1.5)
    ax.axvline(1, color="0.4", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlim(0.24, 1.07)
    ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0.3", "0.4", "0.5", "0.6", "0.8", "1.0"])
    ax.set_title("(c) Gain distribution", fontsize=10)
    ax.set_xlabel(r"Particle gain $g_i$")
    ax.set_ylabel(r"Probability density per $\log_{10}g_i$")
    ax.grid(True, linewidth=0.45, alpha=0.30)
    handles.extend((
        Line2D([], [], color="black", linestyle="-", linewidth=1.5, label="Cations"),
        Line2D([], [], color="black", linestyle="--", linewidth=1.5, label="Anions"),
    ))
    ax.legend(handles=handles, frameon=False, fontsize=6.3,
              ncol=1, loc="upper left", handlelength=1.8)


def plot_timing(ax, points_by_case):
    for case, label, color, marker in timing_data.TIMING_CASES:
        points = points_by_case[case]
        ax.plot([point.n_atoms for point in points],
                [point.seconds_per_step for point in points],
                color=color, marker=marker, markersize=5.5,
                markerfacecolor="none", linewidth=1.5, linestyle="--", label=label)
    ax.set_title(r"(d) Computational cost, $P=20$", fontsize=10)
    ax.set_xlabel(r"Number of particles $N$")
    ax.set_ylabel("CPU time per step [s]")
    ax.set_xlim(left=0)
    ax.set_xticks([1000, 10000, 50000])
    ax.set_xticklabels([r"$10^3$", r"$10^4$", r"$5\times10^4$"])
    ax.xaxis.set_minor_locator(FixedLocator([20000, 30000, 40000]))
    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.legend(frameon=False, fontsize=7, loc="upper left", handlelength=1.8)


def main():
    gamma_points = gamma_plot.read_kinetic_points(gamma_plot.DEFAULT_INPUTS,
                                                   KINETIC_CASES)
    timestep_points = timestep_plot.read_kinetic_points(timestep_plot.DEFAULT_INPUTS,
                                                         KINETIC_CASES)
    timing_points = timing_data.read_timing_points(timing_data.DEFAULT_INPUTS)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), dpi=300)
    plot_gamma(axes[0, 0], gamma_points)
    plot_timestep(axes[0, 1], timestep_points)
    plot_gain(axes[1, 0])
    plot_timing(axes[1, 1], timing_points)
    for ax in axes.flat:
        ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.07, top=0.96,
                        wspace=0.20, hspace=0.30)
    fig.savefig(OUTPUT_PATH, dpi=300)
    plt.close(fig)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
