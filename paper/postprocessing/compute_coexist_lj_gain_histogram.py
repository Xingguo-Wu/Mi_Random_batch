#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import numpy as np


TRAJECTORY = "pro_mirbl_state.dump"
SEEDS = "12345,22345,32345,42345,52345,62345,72345,82345"
FIELDNAMES = ("alpha", "gain_left", "gain_right", "count", "density_per_log10_gain")


def matching_runs(input_dir, alpha, rho, nsteps, seeds):
    runs = []
    for path in input_dir.rglob(TRAJECTORY):
        meta_path = path.parent / "run.meta"
        if not meta_path.is_file():
            continue
        meta = dict(line.split("=", 1) for line in meta_path.read_text().splitlines() if "=" in line)
        if meta.get("PAIR_MODE") != "mirbl":
            continue
        if (int(meta["NSTEPS"]) == nsteps
                and math.isclose(float(meta["ALPHA"]), alpha, abs_tol=1e-10)
                and math.isclose(float(meta["RHO_MIRBL"]), rho, abs_tol=1e-10)):
            runs.append((int(meta["PAIR_SEED_BASE"]), path))
    return [path for _, path in sorted(runs)]


def gain_frames(path, start, end):
    with path.open(encoding="utf-8") as source:
        while source.readline():
            step = int(source.readline())
            if step > end:
                return
            source.readline()
            count = int(source.readline())
            for _ in range(4):
                source.readline()
            columns = source.readline().split()[2:]
            gain_index = columns.index("d_mirbl_gain")
            if step < start:
                for _ in range(count):
                    source.readline()
                continue
            values = np.fromstring("".join(source.readline() for _ in range(count)), sep=" ")
            yield values.reshape(count, len(columns))[:, gain_index]


def histogram(gains):
    edges = np.logspace(
        np.log10(min(gains.min(), 1.0)) - 0.03,
        np.log10(max(gains.max(), 1.0)) + 0.03,
        101,
    )
    counts, _ = np.histogram(gains, bins=edges)
    density = counts / (len(gains) * np.diff(np.log10(edges)))
    return edges, counts, density


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--alpha", required=True, type=float, nargs="+")
    parser.add_argument("--rho", required=True, type=float)
    parser.add_argument("--nsteps", required=True, type=int)
    parser.add_argument("--start-step", required=True, type=int)
    parser.add_argument("--end-step", required=True, type=int)
    parser.add_argument("--seeds", default=SEEDS)
    args = parser.parse_args()
    seeds = [int(seed) for seed in args.seeds.split(",")]

    rows = []
    for alpha in args.alpha:
        runs = matching_runs(args.input_dir, alpha, args.rho, args.nsteps, seeds)
        samples = [frame for path in runs
                   for frame in gain_frames(path, args.start_step, args.end_step)]
        gains = np.concatenate(samples)
        edges, counts, density = histogram(gains)
        rows.extend((alpha, left, right, count, value)
                    for left, right, count, value in
                    zip(edges[:-1], edges[1:], counts, density))
        print(f"alpha={alpha:g}: {len(runs)} runs, {len(gains)} gains")
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(FIELDNAMES)
        writer.writerows(rows)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
