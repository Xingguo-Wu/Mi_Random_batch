
from __future__ import annotations

import csv
from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"


# ============================================================
# 1. Data paths and parameters
# ============================================================

METHOD_INPUTS = {
    p: DATA_DIR / f"msd_method_curvesP{p}.tsv"
    for p in (5, 10, 15)
}

ABLATION_INPUTS = {
    p: DATA_DIR / f"msd_ablation_curvesP{p}.tsv"
    for p in (5, 10)
}

GAIN_INPUTS = {
    p: DATA_DIR / f"coexist_lj_gain_overall_data_P{p}.csv"
    for p in (5, 10)
}

OUTPUT_PATH = SCRIPT_DIR.parent / "figures" / "coelj_all_six_panels.png"

RBL_CASE = "mi-rbl_a0.0"
VRRBL_CASE = "vrrbl"

MI_ALPHA_CASES = {
    5: (
        ("mi-rbl_a0.58_rho0.95", "0.58"),
        ("mi-rbl_a0.59_rho0.95", "0.59"),
        ("mi-rbl_a0.60_rho0.95", "0.60"),
    ),
    10: (
        ("mi-rbl_a0.36_rho0.95", "0.36"),
        ("mi-rbl_a0.37_rho0.95", "0.37"),
        ("mi-rbl_a0.38_rho0.95", "0.38"),
    ),
    15: (
        ("mi-rbl_a0.26_rho0.95", "0.26"),
        ("mi-rbl_a0.27_rho0.95", "0.27"),
        ("mi-rbl_a0.28_rho0.95", "0.28"),
    ),
}

COLORS = {
    "standard": "black",
    "RBL": "#1f77b4",
    "VR-RBL": "#ff7f0e",
    "Mi-RBL": ("#2ca02c", "#d62728", "#9467bd"),
}
GAIN_COLORS = {
    5: ("#1f77b4", "#ff7f0e", "#2ca02c"),
    10: ("#d62728", "#9467bd", "#8c564b"),
}

MARKERS = {
    "RBL": "s",
    "VR-RBL": "D",
    "Mi-RBL": ("^", "o", "v"),
}

MARK_EVERY = 3
INSET_X_LIMIT = (-0.1, 5.7)
INSET_BOUNDS = [0.65, 0.07, 0.34, 0.32]

GAIN_INSETS = {
    5: {
        "xlim": (0.48, 0.66),
        "ylim": (2.5, 4.6),
        "xticks": (0.5, 0.6),
        "yticks": (2.5, 3.5, 4.5),
        "bounds": (0.30, 0.25, 0.26, 0.22),
    },
    10: {
        "xlim": (0.70, 0.88),
        "ylim": (5.2, 8.5),
        "xticks": (0.7, 0.8),
        "yticks": (6.0, 7.0, 8.0),
        "bounds": (0.30, 0.70, 0.26, 0.22),
    },
}


# ============================================================
# 2. MSD plotting
# ============================================================

def read_msd(path):
    curves = {}
    with path.open(encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            time, msd = curves.setdefault(row["method"], ([], []))
            time.append(float(row["time"]))
            msd.append(float(row["msd_mean"]))
    return curves


def plot_msd(ax, curves, p, title):

    series = [
        ("standard", "CL", COLORS["standard"], None, "-", 1.57),
        (RBL_CASE, "RBL", COLORS["RBL"], MARKERS["RBL"], "None", 1.0),
        (VRRBL_CASE, "VR-RBL", COLORS["VR-RBL"], MARKERS["VR-RBL"], "None", 1.0),
    ]

    for i, (case, alpha) in enumerate(MI_ALPHA_CASES[p]):
        series.append((
            case,
            rf"Mi-RBL($\alpha={alpha}$)",
            COLORS["Mi-RBL"][i],
            MARKERS["Mi-RBL"][i],
            "None",
            1.0,
        ))

    for method, label, color, marker, linestyle, linewidth in series:
        time, msd = curves[method]

        ax.plot(
            time, msd,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            marker=marker,
            markersize=5.5,
            markerfacecolor="none" if marker else color,
            markeredgewidth=0.9,
            markevery=MARK_EVERY,
            label=label,
        )

    ax.set_title(title, fontsize=10, loc="center")
    ax.set_xlabel(r"$t$ [$\tau_{\mathrm{LJ}}$]")
    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.tick_params(labelsize=8)
    ax.legend(frameon=False, fontsize=7.0, loc="upper left", handlelength=1.8)

    # Original MSD inset
    inset = ax.inset_axes(INSET_BOUNDS)

    for method, label, color, marker, linestyle, linewidth in series:
        time, msd = curves[method]

        inset.plot(
            time, msd,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            marker=marker,
            markersize=5.5,
            markerfacecolor="none" if marker else color,
            markeredgewidth=0.9,
            markevery=1,
        )

    inset.set_xlim(*INSET_X_LIMIT)
    inset.set_ylim([-0.9, 9.3])
    inset.set_yticks([0.0, 2.5, 5.0, 7.5])
    inset.tick_params(labelsize=6.2, length=2.2, pad=1.2)
    inset.grid(True, linewidth=0.35, alpha=0.28)

    for spine in inset.spines.values():
        spine.set_linewidth(0.65)


# ============================================================
# 3. Gain distribution
# ============================================================

def read_gain(path):

    bins = {}

    with path.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            alpha = float(row["alpha"])
            bins.setdefault(alpha, []).append((
                float(row["gain_left"]),
                float(row["gain_right"]),
                float(row["density_per_log10_gain"]),
            ))

    return {
        alpha: (
            [sqrt(left * right) for left, right, _ in values],
            [density for _, _, density in values],
        )
        for alpha, values in bins.items()
    }


def plot_gain(ax, gain):

    for p in (5, 10):
        for index, (_, alpha_text) in enumerate(MI_ALPHA_CASES[p]):
            alpha = float(alpha_text)
            centers, density = gain[p][alpha]

            ax.plot(
                centers,
                density,
                color=GAIN_COLORS[p][index],
                linewidth=1.6,
                label=rf"$P={p},\ \alpha={alpha_text}$",
            )

    ax.axvline(1, color="0.35", linestyle="--", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_title("(f) Gain distribution", fontsize=10, loc="center")
    ax.set_xlabel(r"Particle gain $g_i$")
    ax.set_ylabel(r"Probability density per $\log_{10}g_i$")

    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()

    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.tick_params(labelsize=11)
    ax.legend(frameon=False, fontsize=7.0, loc="upper right", handlelength=1.8)

    # Gain insets
    for p in (5, 10):
        cfg = GAIN_INSETS[p]
        inset = ax.inset_axes(cfg["bounds"])

        for index, (_, alpha_text) in enumerate(MI_ALPHA_CASES[p]):
            alpha = float(alpha_text)
            centers, density = gain[p][alpha]
            inset.plot(centers, density, color=GAIN_COLORS[p][index], linewidth=1.2)

        inset.set_xscale("log")
        inset.set_xlim(*cfg["xlim"])
        inset.set_ylim(*cfg["ylim"])
        inset.set_xticks(cfg["xticks"])
        inset.set_yticks(cfg["yticks"])
        inset.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
        inset.xaxis.set_minor_formatter(mticker.NullFormatter())
        inset.grid(True, linewidth=0.35, alpha=0.28)
        inset.tick_params(labelsize=6.2, length=2.2, pad=1.2)

        for spine in inset.spines.values():
            spine.set_linewidth(0.65)

        inset.set_title(rf"$P={p}$", fontsize=6.4, loc="left", pad=2)


# ============================================================
# 4. Main figure: 2 x 3 panels
# ============================================================

def main():

    fig, axes = plt.subplots(
        2, 3,
        figsize=(12.2, 7.0),
        dpi=300,
    )

    # Top row: formal MSD (a), (b), (c)
    for ax, p, label in zip(axes[0], (5, 10, 15), ("a", "b", "c")):
        curves = read_msd(METHOD_INPUTS[p])
        plot_msd(ax, curves, p, rf"({label}) $P={p}$")

    axes[0, 0].set_ylabel(r"$\mathrm{MSD}$ [$\sigma_{\mathrm{LJ}}^2$]")
    axes[0, 1].tick_params(labelleft=False)
    axes[0, 2].tick_params(labelleft=False)

    # Bottom row: ablation MSD (d), (e)
    for ax, p, label in zip(axes[1, :2], (5, 10), ("d", "e")):
        curves = read_msd(ABLATION_INPUTS[p])
        plot_msd(ax, curves, p, rf"({label}) Ablation, $P={p}$")

    axes[1, 0].set_ylabel(r"$\mathrm{MSD}$ [$\sigma_{\mathrm{LJ}}^2$]")
    axes[1, 1].tick_params(labelleft=False)

    # Bottom right: gain distribution (f)
    gain = {p: read_gain(GAIN_INPUTS[p]) for p in (5, 10)}
    plot_gain(axes[1, 2], gain)

    # Align the three columns and two rows.
    fig.subplots_adjust(
        left=0.075,
        right=0.992,
        bottom=0.08,
        top=0.96,
        wspace=0.10,
        hspace=0.30,
    )

    fig.savefig(OUTPUT_PATH, bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
