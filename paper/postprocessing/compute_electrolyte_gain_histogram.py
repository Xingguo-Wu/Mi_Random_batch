#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import numpy as np


DUMP_NAME = "mirbe_state.dump"
GROUPS = ("type_1", "type_2")
STYLE = "mi/rbe/lag"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--p", type=int, choices=(5, 10), required=True)
    parser.add_argument("--alpha", type=float, nargs="+", required=True)
    parser.add_argument("--rho", type=float, required=True, help="Mi-RBE EMA parameter")
    parser.add_argument("--nsteps", type=int, required=True, help="Production steps")
    parser.add_argument("--start-step", type=int, required=True)
    parser.add_argument("--end-step", type=int, required=True)
    parser.add_argument("--seeds", help="Comma-separated TRAJ_SEED values; default: all matches")
    return parser.parse_args()


def find_runs(args, alpha):
    requested = set(map(int, args.seeds.split(","))) if args.seeds else None
    runs = []
    for path in args.input_dir.rglob(DUMP_NAME):
        meta_path = path.parent / "run.meta"
        if not meta_path.is_file():
            continue
        meta = dict(line.split("=", 1) for line in meta_path.read_text().splitlines()
                    if "=" in line)
        if meta.get("SYSTEM") != "electrolyte_bulk" or meta.get("KSPACE_STYLE") != STYLE:
            continue
        if meta.get("MIRBE_DUMP_GAIN") != "1":
            continue
        if (int(meta["RBE_P"]) != args.p or int(meta["NSTEP_PROD"]) != args.nsteps
                or not math.isclose(float(meta["MIRBE_ALPHA"]), alpha, abs_tol=1e-10)
                or not math.isclose(float(meta["MIRBE_RHO"]), args.rho, abs_tol=1e-10)):
            continue
        seed = int(meta["TRAJ_SEED"])
        if requested is None or seed in requested:
            runs.append((seed, path))
    if not runs:
        raise ValueError("No matching Mi-RBE gain runs")
    return sorted(runs, key=lambda item: item[0])


def selected_frames(path, start, end):
    with path.open(encoding="utf-8") as handle:
        while handle.readline():  # ITEM: TIMESTEP
            step = int(handle.readline())
            if step > end:
                return
            handle.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(handle.readline())
            handle.readline()  # ITEM: BOX BOUNDS
            for _ in range(3):
                handle.readline()
            columns = handle.readline().split()[2:]
            if step < start:
                for _ in range(natoms):
                    handle.readline()
                continue
            values = np.fromstring("".join(handle.readline() for _ in range(natoms)), sep=" ")
            values = values.reshape(natoms, len(columns))
            types = values[:, columns.index("type")].astype(int)
            gains = values[:, columns.index("d_mirbe_gain")]
            yield step, types, gains


def analyze_alpha(args, alpha):
    runs = find_runs(args, alpha)
    pooled = {group: [] for group in GROUPS}
    for seed, path in runs:
        frame_count = 0
        for step, types, gains in selected_frames(path, args.start_step, args.end_step):
            frame_count += 1
            pooled["type_1"].append(gains[types == 1])
            pooled["type_2"].append(gains[types == 2])
        if not frame_count:
            raise ValueError(f"No frames in requested window: {path}")
        print(f"seed={seed}: {frame_count} frames")

    pooled = {group: np.concatenate(values) for group, values in pooled.items()}
    min_gain = min(values.min() for values in pooled.values())
    max_gain = max(values.max() for values in pooled.values())
    edges = np.logspace(np.log10(min(min_gain, 1.0)) - 0.03,
                        np.log10(max(max_gain, 1.0)) + 0.03, 101)
    histogram_rows = []
    for group in GROUPS:
        gains = pooled[group]
        counts, _ = np.histogram(gains, bins=edges)
        density = counts / (len(gains) * np.diff(np.log10(edges)))
        for left, right, count, value in zip(edges[:-1], edges[1:], counts, density):
            histogram_rows.append({"group": group, "alpha": alpha, "gain_left": left,
                                   "gain_right": right, "count": count,
                                   "density_per_log10_gain": value})
    return histogram_rows


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = [row for alpha in args.alpha for row in analyze_alpha(args, alpha)]
    with (args.output_dir / f"electrolyte_gain_P{args.p}.csv").open("w", newline="") as handle:
        fields = ["group", "alpha", "gain_left", "gain_right", "count",
                  "density_per_log10_gain"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
