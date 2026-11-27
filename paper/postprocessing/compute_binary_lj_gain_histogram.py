#!/usr/bin/env python3
"""Compute type-resolved Mi-RBL gain histograms for the binary LJ mixture."""

import argparse
import csv
import math
from pathlib import Path

import numpy as np


TRAJECTORY = "pro_mirbl_state.dump"
SEEDS = (12345, 22345, 32345, 42345, 52345, 62345, 72345, 82345)


def find_runs(root: Path, args) -> list[tuple[int, Path]]:
    runs = {}
    for path in root.rglob(TRAJECTORY):
        meta_path = path.parent / "run.meta"
        if not meta_path.is_file():
            continue
        meta = dict(line.split("=", 1) for line in meta_path.read_text().splitlines()
                    if "=" in line)
        if meta.get("PAIR_MODE") != "mirbl":
            continue
        if (int(meta["P"]) != args.p or int(meta["PROD_STEPS"]) != args.nsteps
                or not math.isclose(float(meta["MI_ALPHA"]), args.alpha, abs_tol=1e-10)
                or not math.isclose(float(meta["MI_RHO"]), args.rho, abs_tol=1e-10)):
            continue
        seed = int(meta["PAIR_SEED_BASE"])
        if seed in SEEDS:
            runs[seed] = path
    return [(seed, runs[seed]) for seed in SEEDS]


def selected_frames(path: Path, start: int, end: int):
    with path.open(encoding="utf-8") as handle:
        while handle.readline():
            step = int(handle.readline())
            if step > end:
                return
            handle.readline()
            natoms = int(handle.readline())
            handle.readline()
            for _ in range(3):
                handle.readline()
            columns = handle.readline().split()[2:]
            if step < start:
                for _ in range(natoms):
                    handle.readline()
                continue
            values = np.fromstring("".join(handle.readline() for _ in range(natoms)), sep=" ")
            values = values.reshape(natoms, len(columns))
            yield values[:, columns.index("type")].astype(int), values[:, columns.index("d_mirbl_gain")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--rho", type=float, default=0.95, help="EMA beta in run.meta")
    parser.add_argument("--start-step", type=int, required=True)
    parser.add_argument("--end-step", type=int, required=True)
    parser.add_argument("--p", type=int, default=10)
    parser.add_argument("--nsteps", type=int, default=400000)
    args = parser.parse_args()

    pooled = {"type_1": [], "type_2": []}
    for seed, path in find_runs(args.input_dir, args):
        frame_count = 0
        for types, gains in selected_frames(path, args.start_step, args.end_step):
            pooled["type_1"].append(gains[types == 1])
            pooled["type_2"].append(gains[types == 2])
            frame_count += 1
        print(f"seed={seed}: {frame_count} frames")

    pooled = {group: np.concatenate(chunks) for group, chunks in pooled.items()}
    edges = np.linspace(min(min(values.min() for values in pooled.values()), 1.0) * 0.95,
                        max(max(values.max() for values in pooled.values()), 1.0) * 1.05, 101)
    rows = []
    for group, values in pooled.items():
        counts, _ = np.histogram(values, bins=edges)
        density = counts / (len(values) * np.diff(edges))
        for left, right, value in zip(edges[:-1], edges[1:], density):
            rows.append((group, left, right, value))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    alpha_tag = format(args.alpha, ".12g").replace(".", "p")
    output = args.output_dir / f"binlj_gain_{alpha_tag}.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("group", "gain_left", "gain_right", "density_per_gain"))
        writer.writerows(rows)
    print(output)


if __name__ == "__main__":
    main()
