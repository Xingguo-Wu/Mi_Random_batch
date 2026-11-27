from __future__ import annotations

import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_PATH = SCRIPT_DIR.parent / "figures" / "binary_lj_kinetic_error.png"
SEEDS = (12345, 22345, 32345, 42345, 52345, 62345, 72345, 82345)
ALPHA_TAGS = (25, 28, 31, 34, 37)
BETA_TAGS = (91, 95, 99)

MIRBL_COLORS = {0.91: "#2ca02c", 0.95: "#d62728", 0.99: "#9467bd"}
RBL_COLOR = "#1f77b4"
VRRBL_COLOR = "#ff7f0e"
MARKERS = {0.91: "^", 0.95: "o", 0.99: "v"}


def read_stats(path: Path) -> dict[str, tuple[float, float]]:
    rows = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) == 7 and fields[1] == "kinetic":
            rows[fields[0]] = (float(fields[5]), float(fields[6]))
    return rows


def summarize(rows: dict[str, tuple[float, float]], cases: list[str]) -> tuple[float, float]:
    values = [rows[case] for case in cases]
    return statistics.fmean(row[0] for row in values), statistics.fmean(row[1] for row in values)


def load_kinetic_errors():
    cl = read_stats(DATA_DIR / "binary_standard_stats.txt")
    rbl = read_stats(DATA_DIR / "binary_mirbl_stats1.txt")
    vrrbl = read_stats(DATA_DIR / "binary_vrrbl_stats.txt")
    mirbl = read_stats(DATA_DIR / "binary_mirbl_stats2.txt")

    cl_mean, cl_std = summarize(cl, [f"cl_s{seed}" for seed in SEEDS])
    rbl_mean, rbl_std = summarize(rbl, [f"rbl_a0p0_b0p95_s{seed}" for seed in SEEDS])
    vrrbl_mean, vrrbl_std = summarize(vrrbl, [f"vrrbl_s{seed}" for seed in SEEDS])

    mean_deltas = {}
    std_deltas = {}
    for beta_tag in BETA_TAGS:
        beta = beta_tag / 100
        mean_deltas[beta] = []
        std_deltas[beta] = []
        for alpha_tag in ALPHA_TAGS:
            alpha = alpha_tag / 100
            mean, std = summarize(
                mirbl,
                [f"mirbl_a0p{alpha_tag}_b0p{beta_tag}_s{seed}" for seed in SEEDS],
            )
            mean_deltas[beta].append((alpha, mean - cl_mean))
            std_deltas[beta].append((alpha, std - cl_std))
    return (
        mean_deltas,
        rbl_mean - cl_mean,
        vrrbl_mean - cl_mean,
        std_deltas,
        rbl_std - cl_std,
        vrrbl_std - cl_std,
    )


def plot_panel(ax, deltas_by_rho, rbl_delta, vrrbl_delta, ylabel, title, inset_ylim, scientific_y=False):
    sorted_series = []
    for rho in sorted(deltas_by_rho):
        points = sorted(deltas_by_rho[rho])
        alphas = [point[0] for point in points]
        deltas = [point[1] for point in points]
        sorted_series.append((rho, alphas, deltas))
        ax.plot(
            alphas,
            deltas,
            marker=MARKERS[rho],
            linewidth=1.7,
            markersize=5.5,
            markerfacecolor="none",
            markeredgewidth=1.0,
            color=MIRBL_COLORS[rho],
            label=fr"Mi-RBL $\beta={rho:.2f}$",
        )

    xmin = ALPHA_TAGS[0] / 100
    xmax = ALPHA_TAGS[-1] / 100
    pad = 0.006
    ax.axhline(0.0, color="black", linewidth=1.25, alpha=0.9, label="CL")
    ax.hlines(vrrbl_delta, xmin - pad, xmax + pad, colors=VRRBL_COLOR, linestyles="--", linewidth=2.0, label="VR-RBL")
    ax.hlines(rbl_delta, xmin - pad, xmax + pad, colors=RBL_COLOR, linestyles=":", linewidth=2.0, label="RBL")

    inset_bounds = [0.52, 0.49, 0.44, 0.43]
    inset = ax.inset_axes(inset_bounds)
    for rho, alphas, deltas in sorted_series:
        inset.plot(
            alphas,
            deltas,
            marker=MARKERS[rho],
            linewidth=1.7,
            markersize=5.5,
            markerfacecolor="none",
            markeredgewidth=1.0,
            color=MIRBL_COLORS[rho],
        )
    inset.axhline(0.0, color="black", linewidth=1.0, alpha=0.85)
    inset.set_xlim(xmin - pad, xmax + pad)
    inset.set_ylim(*inset_ylim)
    inset.tick_params(labelsize=7, length=2.5)
    inset.grid(True, linewidth=0.4, alpha=0.30)

    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, linewidth=0.5, alpha=0.35)
    handles, labels = ax.get_legend_handles_labels()
    legend_order = {"CL": 0, "RBL": 1, "VR-RBL": 2}
    ordered = sorted(
        zip(handles, labels),
        key=lambda item: legend_order.get(item[1], 3),
    )
    ax.legend(
        *zip(*ordered),
        frameon=False,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.93),
        ncol=1,
        borderaxespad=0.0,
        handlelength=2.2,
        columnspacing=0.9,
    )

    if scientific_y:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        inset.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value * 1.0e3:g}"))
        inset.yaxis.get_offset_text().set_visible(False)
        left, bottom, width, height = inset_bounds
        ax.text(
            left - 0.01,
            bottom + height,
            r"$\times 10^{-3}$",
            transform=ax.transAxes,
            fontsize=7,
            ha="right",
            va="top",
            clip_on=False,
        )


def main() -> None:
    mean_deltas, rbl_mean_delta, vrrbl_mean_delta, std_deltas, rbl_std_delta, vrrbl_std_delta = load_kinetic_errors()

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3), dpi=180)
    plot_panel(
        axes[0],
        mean_deltas,
        rbl_mean_delta,
        vrrbl_mean_delta,
        r"Per-atom kinetic energy error [$\epsilon_{11}$]",
        "(a) Kinetic energy",
        inset_ylim=(-0.059, 0.089),
    )
    plot_panel(
        axes[1],
        std_deltas,
        rbl_std_delta,
        vrrbl_std_delta,
        r"Per-atom kinetic standard deviation error [$\epsilon_{11}$]",
        "(b) Kinetic standard deviation",
        inset_ylim=(-0.00035, 0.00121),
        scientific_y=True,
    )

    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.15, top=0.90, wspace=0.19)
    fig.savefig(OUTPUT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()




