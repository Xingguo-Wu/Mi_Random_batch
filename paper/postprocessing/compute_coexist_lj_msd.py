#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


TRAJECTORY = "traj_unwrapped.lammpstrj"
DEFAULT_ORIGINS = tuple(range(400000, 490001, 10000))


@dataclass
class Frame:
    step: int
    ids: list[int]
    xyz: list[tuple[float, float, float]]
    com: tuple[float, float, float]


def frames(path: Path):
    with path.open(encoding="utf-8") as source:
        while line := source.readline():
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(source.readline())
            source.readline()
            count = int(source.readline())
            for _ in range(4):
                source.readline()
            columns = source.readline().split()[2:]
            indices = [columns.index(name) for name in ("id", "xu", "yu", "zu")]
            atoms = []
            for _ in range(count):
                row = source.readline().split()
                atoms.append((int(row[indices[0]]), *(float(row[i]) for i in indices[1:])))
            atoms.sort(key=lambda atom: atom[0])
            xyz = [atom[1:] for atom in atoms]
            yield Frame(step, [atom[0] for atom in atoms], xyz, center(xyz))


def center(xyz):
    n = len(xyz)
    return tuple(sum(atom[axis] for atom in xyz) / n for axis in range(3))


def msd(reference: Frame, current: Frame) -> float:
    if current.ids != reference.ids:
        raise ValueError(f"atom IDs changed between steps {reference.step} and {current.step}")
    c0, c1 = reference.com, current.com
    drift = tuple(c1[i] - c0[i] for i in range(3))
    return sum(
        sum((atom[i] - initial[i] - drift[i]) ** 2 for i in range(3))
        for initial, atom in zip(reference.xyz, current.xyz)
    ) / len(reference.ids)


def trajectory_curve(path: Path, origins: tuple[int, ...], window: int) -> dict[int, float]:
    references: dict[int, Frame] = {}
    values: dict[int, list[float]] = defaultdict(list)
    last_step = max(origins) + window
    for frame in frames(path):
        if frame.step in origins:
            references[frame.step] = frame
        for origin, reference in references.items():
            if origin <= frame.step <= origin + window:
                values[frame.step - origin].append(msd(reference, frame))
        if frame.step > last_step:
            break
    expected = set(range(0, window + 1, 100))
    if set(values) != expected or any(len(series) != len(origins) for series in values.values()):
        raise ValueError(f"{path}: incomplete 100-step MSD windows")
    return {lag: float(f"{sum(series) / len(series):.12g}") for lag, series in values.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", action="append", required=True, metavar="NAME=DIR1,DIR2,...")
    parser.add_argument("--output", required=True, type=Path, help="Final msd_*_curvesP*.tsv file")
    parser.add_argument("--t0", default=",".join(map(str, DEFAULT_ORIGINS)))
    parser.add_argument("--window-length", type=int, default=8000)
    parser.add_argument("--time-step", type=float, default=0.01)
    args = parser.parse_args()
    origins = tuple(int(value) for value in args.t0.split(","))
    rows = ["method\ttrajectory_count\tlag_step\ttime\tmsd_mean"]
    for spec in args.method:
        name, raw_paths = spec.split("=", 1)
        inputs = [Path(item) for item in raw_paths.split(",")]
        curves = []
        for input_path in inputs:
            path = input_path if input_path.is_file() else input_path / TRAJECTORY
            curves.append(trajectory_curve(path, origins, args.window_length))
        for lag in range(0, args.window_length + 1, 100):
            mean = sum(curve[lag] for curve in curves) / len(curves)
            time = lag * args.time_step
            rows.append(f"{name}\t{len(curves)}\t{lag}\t{time:g}\t{mean:.12e}")
        print(f"{name}: {len(curves)} trajectories, {len(curves[0])} lag points")
    args.output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
