from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_PATH = SCRIPT_DIR.parent / "figures" / "electrolyte_charge_density.png"

DEFAULT_INPUTS = {
    5: DATA_DIR / "electrolyte_densityP5.dat",
    10: DATA_DIR / "electrolyte_densityP10.dat",
    20: DATA_DIR / "electrolyte_densityP20.dat",
}

P_VALUES = (5, 10, 20)
REFERENCE_CASE = "pppm"
REFERENCE_LABEL = "PPPM"
RBE_CASE = "mirbe_a0.0_rho0.95"
MI_ALPHA_CASES = {
    5: (
        ("mirbe_a0.4_rho0.95", "0.4"),
        ("mirbe_a0.5_rho0.95", "0.5"),
        ("mirbe_a0.6_rho0.95", "0.6"),),
    10: (
        ("mirbe_a0.4_rho0.95", "0.4"),
        ("mirbe_a0.5_rho0.95", "0.5"),
        ("mirbe_a0.6_rho0.95", "0.6"),),
    20: (
        ("mirbe_a0.1_rho0.95", "0.1"),
        ("mirbe_a0.2_rho0.95", "0.2"),
        ("mirbe_a0.3_rho0.95", "0.3"),),
}

COLORS = {
    "Reference": "black",
    "RBE": "#1f77b4",
    "Mi-RBE": ("#2ca02c", "#d62728", "#9467bd"),
}

MARKERS = {
    "RBE": "s",
    "Mi-RBE": ("^", "o", "v"),
}

MARK_EVERY_EARLY = 2
MARK_EVERY_MIDDLE = 1
MARK_EVERY_LATE = 2
MARK_EVERY_SWITCH_X1 = 0.9
MARK_EVERY_SWITCH_X2 = 1.30

X_LIMIT = (0.0, 4.0)
INSET_LIMITS = {
    5: ((1.0, 1.2), (-0.58, -0.46)),
    10: ((1.0, 1.2), (-0.58, -0.46)),
    20: ((1.0, 1.2), (-0.58, -0.46)),
}

def adaptive_markevery(r: list[float]) -> list[int]:
    indices = []
    early_count = 0
    middle_count = 0
    late_count = 0

    for i, x in enumerate(r):
        if x <= MARK_EVERY_SWITCH_X1:
            if early_count % MARK_EVERY_EARLY == 0:
                indices.append(i)
            early_count += 1
        elif x <= MARK_EVERY_SWITCH_X2:
            if middle_count % MARK_EVERY_MIDDLE == 0:
                indices.append(i)
            middle_count += 1
        else:
            if late_count % MARK_EVERY_LATE == 0:
                indices.append(i)
            late_count += 1

    return indices

def read_density(path: Path) -> dict[str, tuple[list[float], list[float]]]:
    curves: dict[str, tuple[list[float], list[float]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            case = fields[0]
            r = float(fields[1])
            rho_q = float(fields[2])
            curves.setdefault(case, ([], []))
            curves[case][0].append(r)
            curves[case][1].append(rho_q)
    return curves


def configured_series(p_value: int):
    yield REFERENCE_CASE, REFERENCE_LABEL, COLORS["Reference"], None, "-", 1.55
    yield RBE_CASE, "RBE", COLORS["RBE"], MARKERS["RBE"], "None", 1.0
    for index, (case, alpha) in enumerate(MI_ALPHA_CASES[p_value]):
        yield (
            case,
            rf"Mi-RBE($\alpha={alpha}$)",
            COLORS["Mi-RBE"][index],
            MARKERS["Mi-RBE"][index],
            "None",
            1.0,
        )


def plot_curves(ax, curves, p_value: int, add_labels: bool) -> None:
    for case, label, color, marker, linestyle, linewidth in configured_series(p_value):
        r, rho_q = curves[case]
        ax.plot(
            r,
            rho_q,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            marker=marker,
            markersize=5.0,
            markerfacecolor="none" if marker else color,
            markeredgewidth=0.9,
            markevery=adaptive_markevery(r),
            label=label if add_labels else None,
        )


def add_inset(ax, curves, p_value: int) -> None:
    xlim, ylim = INSET_LIMITS[p_value]
    inset = ax.inset_axes([0.60, 0.39, 0.36, 0.44])
    plot_curves(inset, curves, p_value, add_labels=False)
    inset.set_xlim(*xlim)
    inset.set_ylim(*ylim)
    inset.tick_params(labelsize=7, length=2.5, pad=1.5)
    inset.grid(True, linewidth=0.4, alpha=0.28)
    for spine in inset.spines.values():
        spine.set_linewidth(0.7)


def build_figure(input_paths: dict[int, Path], output_path: Path) -> None:
    fig, axes = plt.subplots(1, len(P_VALUES), figsize=(10.8, 3.35), dpi=300, sharey=True)

    for panel_label, ax, p_value in zip(("a", "b", "c"), axes, P_VALUES):
        curves = read_density(input_paths[p_value])
        plot_curves(ax, curves, p_value, add_labels=True)
        add_inset(ax, curves, p_value)

        ax.set_title(rf"({panel_label}) $P={p_value}$", fontsize=10, loc="center")
        ax.set_xlim(*X_LIMIT)
        ax.set_xlabel(r"$r$ [$\sigma_{\mathrm{LJ}}$]")
        ax.grid(True, linewidth=0.45, alpha=0.30)
        ax.tick_params(labelsize=8)
        ax.legend(frameon=False, fontsize=6.7, loc="lower right", handlelength=1.8)

    axes[0].set_ylabel(r"$\rho_q^+(r)$")
    fig.subplots_adjust(left=0.090, right=0.992, bottom=0.13, top=0.95, wspace=0.16)
    fig.savefig(output_path, bbox_inches="tight")
    print(f"Wrote {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot cation-centered electrolyte charge density panels from processed data."
    )
    parser.add_argument("--p5", type=Path, default=DEFAULT_INPUTS[5])
    parser.add_argument("--p10", type=Path, default=DEFAULT_INPUTS[10])
    parser.add_argument("--p20", type=Path, default=DEFAULT_INPUTS[20])
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_figure({5: args.p5, 10: args.p10, 20: args.p20}, args.output)


if __name__ == "__main__":
    main()
