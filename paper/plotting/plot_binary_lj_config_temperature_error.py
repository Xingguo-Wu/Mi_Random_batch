#!/usr/bin/env python3
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_PATH = SCRIPT_DIR.parent / "figures" / "binary_lj_config_temperature_error.png"
ALPHA_TAGS = (25, 28, 31, 34, 37)


def read_mean_errors(path: Path) -> dict[tuple[str, str], float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (row["method"], row["species"]): float(row["delta_t_conf_vs_cl"])
            for row in csv.DictReader(handle)
            if row["seed"] == "mean"
        }


def main() -> None:
    alphas = [tag / 100 for tag in ALPHA_TAGS]
    data = [read_mean_errors(DATA_DIR / f"binlj_tconf_p10_a0p{tag}_b0p95.csv")
            for tag in ALPHA_TAGS]
    errors = [rows["mirbl", "all"] for rows in data]
    species_errors = {
        species: [rows["mirbl", str(species)] for rows in data]
        for species in (1, 2)
    }
    baseline = data[0]
    rbl_error = baseline["rbl", "all"]
    vrrbl_error = baseline["vrrbl", "all"]
    xmin, xmax = min(alphas) - 0.006, max(alphas) + 0.006

    fig, (ax, species_ax) = plt.subplots(1, 2, figsize=(11.8, 4.3), dpi=180)
    ax.axhline(0.0, color="black", linewidth=1.25, alpha=0.9, label="CL")
    ax.hlines(rbl_error, xmin, xmax, colors="#1f77b4", linestyles=":",
              linewidth=2.0, label="RBL")
    ax.hlines(vrrbl_error, xmin, xmax, colors="#ff7f0e", linestyles="--",
              linewidth=2.0, label="VR-RBL")
    ax.plot(alphas, errors, marker="o", linewidth=1.7, markersize=5.5,
            markerfacecolor="none", markeredgewidth=1.0, color="#d62728",
            label=r"Mi-RBL $\beta=0.95$")
    inset = ax.inset_axes([0.52, 0.46, 0.44, 0.43])
    inset.plot(alphas, errors, marker="o", linewidth=1.7, markersize=5.5,
               markerfacecolor="none", markeredgewidth=1.0, color="#d62728")
    inset.axhline(0.0, color="black", linewidth=1.0, alpha=0.85)
    inset.set_xlim(xmin, xmax)
    inset.tick_params(labelsize=7, length=2.5)
    inset.grid(True, linewidth=0.4, alpha=0.30)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"Configurational temperature error [$\epsilon_{11}/k_B$]")
    ax.set_title("(c) Configurational temperature")
    ax.grid(True, linewidth=0.5, alpha=0.35)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.02, 0.93),
              borderaxespad=0.0, handlelength=2.2, columnspacing=0.9)

    ax = species_ax
    ax.axhline(0.0, color="black", linewidth=1.25, alpha=0.9, label="CL")
    for species, label, color in ((1, "Type 1", "#175987"),
                                  (2, "Type 2", "#62a0cb")):
        rbl_species_error = baseline["rbl", str(species)]
        ax.hlines(rbl_species_error, xmin, xmax, colors=color, linestyles=":",
                  linewidth=2.0, label=f"RBL {label}")
    for species, label, color in ((1, "Type 1", "#bf5f0b"),
                                  (2, "Type 2", "#ffa556")):
        vrrbl_species_error = baseline["vrrbl", str(species)]
        ax.hlines(vrrbl_species_error, xmin, xmax, colors=color, linestyles="--",
                  linewidth=2.0, label=f"VR-RBL {label}")
    inset = ax.inset_axes([0.52, 0.46, 0.44, 0.43])
    for species, label, color in ((1, "Type 1", "#a11d1e"),
                                  (2, "Type 2", "#e26869")):
        ax.plot(alphas, species_errors[species], marker="o", linewidth=1.7, markersize=5.5,
                markerfacecolor="none", markeredgewidth=1.0, color=color,
                label=fr"Mi-RBL {label} $\beta=0.95$")
        inset.plot(alphas, species_errors[species], marker="o", linewidth=1.7, markersize=5.5,
                   markerfacecolor="none", markeredgewidth=1.0, color=color)
    inset.axhline(0.0, color="black", linewidth=1.0, alpha=0.85)
    inset.set_xlim(xmin, xmax)
    inset.tick_params(labelsize=7, length=2.5)
    inset.grid(True, linewidth=0.4, alpha=0.30)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"Species configurational temperature error [$\epsilon_{11}/k_B$]")
    ax.set_title("(d) Species configurational temperature")
    ax.grid(True, linewidth=0.5, alpha=0.35)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.02, 0.915),
              borderaxespad=0.0, handlelength=2.2, columnspacing=0.9)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.15, top=0.90, wspace=0.19)

    fig.savefig(OUTPUT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
