#!/usr/bin/env python
"""Plot the three binary LJ RDF pairs from the published statistics."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_PATH = SCRIPT_DIR.parent / "figures" / "binary_lj_rdf.png"
DEFAULT_INPUTS = {
    "11": DATA_DIR / "rdf_binary_lj_11.txt",
    "12": DATA_DIR / "rdf_binary_lj_12.txt",
    "22": DATA_DIR / "rdf_binary_lj_22.txt",
}

PANEL_ORDER = ("11", "22", "12")
PANEL_LETTERS = {"11": "a", "22": "b", "12": "c"}
TARGET_MIRBL_ALPHA = 0.31
TARGET_MIRBL_RHOS = (0.91, 0.95, 0.99)
CL_COLOR = "black"
RBL_COLOR = "#1f77b4"
VRRBL_COLOR = "#ff7f0e"
MIRBL_COLORS = {0.91: "#2ca02c", 0.95: "#d62728", 0.99: "#9467bd"}
MIRBL_MARKERS = {0.91: "^", 0.95: "o", 0.99: "v"}

INSET_X_LIMITS = {
    "11": (1.070, 1.080),
    "12": (0.870, 0.880),
}
INSET_BOUNDS = [0.42, 0.46, 0.23, 0.47]
INSET_X_TICKS_BY_PAIR = {
    "11": [1.070, 1.075, 1.080],
    "12": [0.870, 0.875, 0.880],
}

MARK_EVERY_BY_PAIR = {
    "11": {"switches": (0.9, 1.20), "strides": (2, 1, 2)},
    "22": {"switches": (0.80, 0.95, 1.2, 1.38, 1.72, 1.85), "strides": (2, 1, 2, 1, 2, 1, 2)},
    "12": {"switches": (0.70, 1.13), "strides": (2, 1, 2)},
}
LEGEND_LOCS = {"11": "upper right", "22": "lower right", "12": "upper right"}


def pair_label(pair: str) -> str:
    if len(pair) == 2 and pair.isdigit():
        return f"{pair[0]}-{pair[1]}"
    return pair


def rdf_label(pair: str) -> str:
    return rf"$g_{{{pair}}}(r)$"


def alpha_rho_values(method_name: str) -> tuple[str | None, str | None]:
    number = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
    alpha = re.search(rf"(?:^|[^A-Za-z0-9])(?:a|alpha)\s*=?\s*{number}", method_name)
    rho = re.search(rf"(?:^|[^A-Za-z0-9])rho\s*=?\s*{number}", method_name)
    return (
        alpha.group(1) if alpha else None,
        rho.group(1) if rho else None,
    )


def should_plot_method(method_name: str) -> bool:
    key = method_name.lower().replace("_", "-")
    if "mi-rbl" not in key and "mirbl" not in key:
        return True

    alpha, rho = alpha_rho_values(method_name)
    if alpha is None:
        return False
    if math.isclose(float(alpha), 0.0):
        return True
    if rho is None:
        return False
    return math.isclose(float(alpha), TARGET_MIRBL_ALPHA) and any(
        math.isclose(float(rho), target_rho) for target_rho in TARGET_MIRBL_RHOS
    )


def plot_label(method_name: str) -> str:
    key = method_name.lower().replace("_", "-")
    alpha, rho = alpha_rho_values(method_name)
    if key == "standard":
        return "CL"
    if key == "vrrbl":
        return "VR-RBL"
    if key == "rbl" or (alpha is not None and math.isclose(float(alpha), 0.0)):
        return "RBL"
    if "mi-rbl" in key or "mirbl" in key:
        if alpha is not None and rho is not None:
            return rf"Mi-RBL ($\beta={rho}$)"
        return "Mi-RBL"
    return method_name


def plot_color(method_name: str) -> str | None:
    key = method_name.lower().replace("_", "-")
    alpha, rho = alpha_rho_values(method_name)
    if key == "standard":
        return CL_COLOR
    if key == "vrrbl":
        return VRRBL_COLOR
    if key == "rbl" or (alpha is not None and math.isclose(float(alpha), 0.0)):
        return RBL_COLOR
    if "mi-rbl" in key or "mirbl" in key:
        if rho is not None:
            return MIRBL_COLORS.get(float(rho))
        return None
    return None


def plot_marker(method_name: str) -> str:
    key = method_name.lower().replace("_", "-")
    alpha, rho = alpha_rho_values(method_name)
    if key == "standard":
        return "None"
    if key == "vrrbl":
        return "D"
    if key == "rbl" or (alpha is not None and math.isclose(float(alpha), 0.0)):
        return "s"
    if "mi-rbl" in key or "mirbl" in key:
        if rho is not None:
            return MIRBL_MARKERS.get(float(rho), "o")
        return "o"
    return "o"


def read_stats(path: Path) -> list[tuple[str, list[float], list[float], list[float]]]:
    series: list[tuple[str, list[float], list[float], list[float]]] = []
    current: tuple[str, list[float], list[float], list[float]] | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split()
        method = parts[0]
        if current is None or current[0] != method:
            current = (method, [], [], [])
            series.append(current)
        current[1].append(float(parts[1]))
        current[2].append(float(parts[2]))
        current[3].append(float(parts[3]))
    return series


def scaled_mark_stride(base_stride: int, marker_stride: int) -> int:
    return max(1, base_stride * max(1, marker_stride))


def adaptive_markevery(pair: str, r: list[float], marker_stride: int) -> list[int]:
    rule = MARK_EVERY_BY_PAIR[pair]
    switches = rule["switches"]
    strides = rule["strides"]
    if len(strides) != len(switches) + 1:
        raise ValueError(f"Invalid markevery rule for pair {pair}: len(strides) must be len(switches) + 1")

    scaled_strides = [scaled_mark_stride(stride, marker_stride) for stride in strides]
    counts = [0 for _ in scaled_strides]
    indices = []

    for index, radius in enumerate(r):
        segment = 0
        while segment < len(switches) and radius > switches[segment]:
            segment += 1

        if counts[segment] % scaled_strides[segment] == 0:
            indices.append(index)
        counts[segment] += 1

    return indices


def legend_key(item: tuple[object, str]) -> tuple[int, float]:
    label = item[1]
    legend_order = {"CL": 0, "RBL": 1, "VR-RBL": 2}
    if label in legend_order:
        return legend_order[label], 0.0
    beta = re.search(r"\\beta=([-+]?\d+(?:\.\d+)?)", label)
    if label.startswith("Mi-RBL") and beta:
        return 3, float(beta.group(1))
    return 4, 0.0


def add_ordered_legend(ax, pair: str) -> None:
    handles, labels = ax.get_legend_handles_labels()
    ordered = sorted(zip(handles, labels), key=legend_key)
    if not ordered:
        return
    handles, labels = zip(*ordered)
    ax.legend(
        handles,
        labels,
        frameon=False,
        fontsize=8.7,
        loc=LEGEND_LOCS[pair],
        handlelength=1.1,
        handletextpad=0.45,
        labelspacing=0.53,
        borderaxespad=0.22,
    )


def inset_ylim(
    series: list[tuple[str, list[float], list[float], list[float]]],
    x_limits: tuple[float, float],
) -> tuple[float, float]:
    y_values = []
    x_min, x_max = x_limits
    for method, r, mean_g, sem_g in series:
        if not should_plot_method(method):
            continue
        for radius, value, error in zip(r, mean_g, sem_g):
            if x_min <= radius <= x_max:
                y_values.extend((value - error, value + error))
    if not y_values:
        return 0.0, 1.0
    y_min = min(y_values)
    y_max = max(y_values)
    pad = 0.12 * (y_max - y_min) if y_max > y_min else 0.05 * max(abs(y_max), 1.0)
    return y_min - pad, y_max + pad


def draw_series(
    ax,
    pair: str,
    series: list[tuple[str, list[float], list[float], list[float]]],
    marker_stride: int,
    connect_non_cl: bool,
    add_labels: bool,
    markersize: float = 7.0,
) -> None:
    for method, r, mean_g, sem_g in series:
        if not should_plot_method(method):
            continue
        markevery = adaptive_markevery(pair, r, marker_stride)
        is_reference = method.lower().replace("_", "-") == "standard"
        line_style = "-" if is_reference or connect_non_cl else "None"
        line_width = 2.0 if is_reference else (1.20 if connect_non_cl else 0.0)
        line = ax.plot(
            r,
            mean_g,
            marker=plot_marker(method),
            markersize=0.0 if is_reference else markersize,
            markerfacecolor=None if is_reference else "none",
            markeredgewidth=0.6 if is_reference else 1.3,
            markevery=markevery,
            linestyle=line_style,
            linewidth=line_width,
            label=plot_label(method) if add_labels else None,
            color=plot_color(method),
        )[0]
        if any(sem_g) and (is_reference or connect_non_cl):
            ax.fill_between(
                r,
                [g - e for g, e in zip(mean_g, sem_g)],
                [g + e for g, e in zip(mean_g, sem_g)],
                color=line.get_color(),
                alpha=0.16,
            )


def add_inset(
    ax,
    pair: str,
    series: list[tuple[str, list[float], list[float], list[float]]],
    marker_stride: int,
    connect_non_cl: bool,
) -> None:
    x_limits = INSET_X_LIMITS.get(pair)
    if x_limits is None:
        return
    inset = ax.inset_axes(INSET_BOUNDS)
    draw_series(
        inset,
        pair,
        series,
        marker_stride,
        connect_non_cl,
        add_labels=False,
        markersize=6.5,
    )
    inset.set_xlim(*x_limits)
    inset.set_xticks(INSET_X_TICKS_BY_PAIR[pair])
    inset.set_ylim(*inset_ylim(series, x_limits))
    inset.tick_params(labelsize=6.0, length=2.0, pad=1.0)
    inset.grid(True, linewidth=0.35, alpha=0.28)
    for spine in inset.spines.values():
        spine.set_linewidth(0.65)


def plot_panel(
    ax,
    pair: str,
    series: list[tuple[str, list[float], list[float], list[float]]],
    marker_stride: int,
    connect_non_cl: bool,
    with_inset: bool,
) -> None:
    draw_series(ax, pair, series, marker_stride, connect_non_cl, add_labels=True)
    ax.set_title(rf"({PANEL_LETTERS[pair]}) {rdf_label(pair)}", fontsize=12.5, loc="center")
    ax.set_xlabel(r"$r$ [$\sigma_{11}$]", fontsize = 12.5)
    ax.set_ylabel(r"$g_{ab}(r)$" if pair == "11" else "", fontsize = 12.5)
    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.tick_params(labelsize=10)
    ax.margins(x=0.01)
    add_ordered_legend(ax, pair)
    if with_inset:
        add_inset(ax, pair, series, marker_stride, connect_non_cl)


def build_figure(
    series_by_pair: dict[str, list[tuple[str, list[float], list[float], list[float]]]],
    output_path: Path,
    dpi: int,
    marker_stride: int,
    connect_non_cl: bool,
    with_inset: bool,
) -> None:
    fig, axes = plt.subplots(1, len(PANEL_ORDER), figsize=(14.4, 4.10), dpi=dpi)
    for ax, pair in zip(axes, PANEL_ORDER):
        plot_panel(ax, pair, series_by_pair[pair], marker_stride, connect_non_cl, with_inset)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.11, top=0.95, wspace=0.14)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rdf-11", type=Path, default=DEFAULT_INPUTS["11"])
    parser.add_argument("--rdf-12", type=Path, default=DEFAULT_INPUTS["12"])
    parser.add_argument("--rdf-22", type=Path, default=DEFAULT_INPUTS["22"])
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--marker-stride", type=int, default=1)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--with-inset",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Draw zoomed inset panels for RDF pairs configured in INSET_X_LIMITS.",
    )
    parser.add_argument(
        "--connect-non-cl",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Connect non-CL RDF markers with lines; use --no-connect-non-cl for markers only.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    input_paths = {"11": args.rdf_11, "12": args.rdf_12, "22": args.rdf_22}
    series_by_pair = {pair: read_stats(input_paths[pair]) for pair in PANEL_ORDER}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    build_figure(
        series_by_pair,
        args.output,
        args.dpi,
        args.marker_stride,
        args.connect_non_cl,
        args.with_inset,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
