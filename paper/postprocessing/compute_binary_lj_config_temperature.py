#!/usr/bin/env python3
"""Calculate KA configurational temperature from production LAMMPS dumps.

The exact lj/cut force is recomputed on each saved configuration.  The
Laplacian includes a narrow-shell estimate at each pair-specific cutoff.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


METHODS = ("cl", "rbl", "vrrbl", "mirbl")
SHELL_WIDTH = 0.02
TRAJECTORY = "traj_unwrapped.lammpstrj"
PAIR_NAMES = ("11", "12", "22")
# Match pair_coeff in binary_lj/in.binary_lj (LJ reduced units).
EPSILON = np.array([1.0, 1.5, 0.5])
SIGMA = np.array([1.00, 0.80, 0.88])
CUTOFF = np.array([2.50, 2.00, 2.20])
q6 = (SIGMA / CUTOFF) ** 6
U_PRIME_CUTOFF = 24 * EPSILON * (q6 - 2 * q6**2) / CUTOFF
SPECIES_WEIGHTS = np.array([[1.0, 0.5, 0.0], [0.0, 0.5, 1.0]])


@dataclass(frozen=True)
class Case:
    method: str
    seed: int
    directory: Path
    trajectory: Path
    meta: dict[str, str]


@dataclass(frozen=True)
class Frame:
    step: int
    bounds: np.ndarray
    types: np.ndarray
    positions: np.ndarray


@dataclass
class Result:
    case: Case
    natoms: int
    type_counts: tuple[int, int]
    box: tuple[float, float, float]
    frames: int
    first_step: int
    last_step: int
    force2: float | np.ndarray
    bulk: float | np.ndarray
    shell_counts: np.ndarray


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Root containing binary_lj run directories with run.meta and trajectories.")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--seeds", required=True, help="Comma/space-separated replica seeds.")
    p.add_argument("--alpha", required=True, type=float, help="Mi-RBL alpha; RBL has alpha=0.")
    p.add_argument("--p", type=int, default=10)
    p.add_argument("--beta", type=float, default=0.95, help="MI_RHO in run.meta.")
    p.add_argument("--nsteps", type=int, default=400000)
    p.add_argument("--start-step", type=int, default=150000)
    p.add_argument("--end-step", type=int, default=350000)
    p.add_argument("--stride", type=int, default=1000)
    return p


def read_meta(path: Path) -> dict[str, str]:
    return dict(line.split("=", 1) for line in path.read_text(encoding="utf-8").splitlines()
                if "=" in line)


def find_cases(args: argparse.Namespace, methods: list[str], seeds: list[int]) -> list[Case]:
    found: dict[tuple[str, int], Case] = {}
    for trajectory in args.input_dir.rglob(TRAJECTORY):
        meta_path = trajectory.parent / "run.meta"
        if not meta_path.is_file():
            continue
        meta = read_meta(meta_path)
        mode = meta.get("PAIR_MODE")
        if mode not in ("standard", "mirbl", "vrrbl"):
            continue
        alpha = float(meta["MI_ALPHA"]) if mode == "mirbl" else None
        method = ("cl" if mode == "standard" else "vrrbl" if mode == "vrrbl"
                  else "rbl" if alpha == 0 else "mirbl")
        seed = int(meta["SEED" if method == "cl" else "PAIR_SEED_BASE"])
        nsteps = int(meta["PROD_STEPS"])
        if method not in methods or seed not in seeds or nsteps != args.nsteps:
            continue
        if method == "mirbl" and not math.isclose(alpha, args.alpha, abs_tol=1e-10):
            continue
        if method in ("rbl", "mirbl") and not math.isclose(
            float(meta["MI_RHO"]), args.beta, abs_tol=1e-10
        ):
            continue
        if method != "cl" and int(meta["P"]) != args.p:
            continue
        key = method, seed
        if key in found:
            raise ValueError(f"Duplicate {method} seed={seed}; use a narrower --input-dir")
        found[key] = Case(method, seed, trajectory.parent, trajectory, meta)
    missing = [(method, seed) for method in methods for seed in seeds
               if (method, seed) not in found]
    return [found[(method, seed)] for method in methods for seed in seeds]


def selected_frames(path: Path, start: int, end: int, stride: int):
    with path.open(encoding="utf-8") as handle:
        while handle.readline():
            step = int(handle.readline())
            if step > end:
                return
            handle.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(handle.readline())
            handle.readline()  # ITEM: BOX BOUNDS pp pp pp
            bounds = np.array([[float(x) for x in handle.readline().split()]
                               for _ in range(3)])
            columns = handle.readline().split()[2:]  # ITEM: ATOMS ...
            selected = start <= step <= end and (step - start) % stride == 0
            if not selected:
                for _ in range(natoms):
                    handle.readline()
                continue
            values = np.fromstring("".join(handle.readline() for _ in range(natoms)), sep=" ")
            values = values.reshape(natoms, len(columns))
            yield Frame(step, bounds, values[:, columns.index("type")].astype(np.int8),
                        values[:, [columns.index(c) for c in ("xu", "yu", "zu")]])


def frame_statistics(frame: Frame, widths: list[float], by_species: bool = False
                     ) -> tuple[float | np.ndarray, float | np.ndarray, np.ndarray]:
    """Return sum |F_i|^2, interior Laplacian, and cutoff-shell pair counts."""
    box = frame.bounds[:, 1] - frame.bounds[:, 0]
    if max(CUTOFF) + max(widths) / 2 >= min(box) / 2:
        raise ValueError(f"Cutoff exceeds half the box at step {frame.step}")
    positions = np.mod(frame.positions - frame.bounds[:, 0], box)
    pairs = cKDTree(positions, boxsize=box).query_pairs(
        max(CUTOFF) + max(widths) / 2, output_type="ndarray")
    forces = np.zeros((len(frame.types), 3))
    bulk = 0.0
    bulk_by_pair = np.zeros(3) if by_species else None
    counts = np.zeros((len(widths), 3), dtype=np.int64)
    for chunk in np.array_split(pairs, max(1, math.ceil(len(pairs) / 100000))):
        if len(chunk) == 0:
            continue
        i, j = chunk[:, 0], chunk[:, 1]
        displacement = positions[i] - positions[j]
        displacement -= box * np.rint(displacement / box)
        r2 = np.einsum("ij,ij->i", displacement, displacement)
        radius = np.sqrt(r2)
        kind = frame.types[i] + frame.types[j] - 2  # 0=11, 1=12, 2=22
        cutoff = CUTOFF[kind]
        inside = radius < cutoff
        if np.any(inside):
            a, b = i[inside], j[inside]
            r2_inside = r2[inside]
            epsilon, sigma = EPSILON[kind[inside]], SIGMA[kind[inside]]
            q6 = (sigma**2 / r2_inside) ** 3
            q12 = q6**2
            pair_force = (displacement[inside]
                          * (24 * epsilon * (2 * q12 - q6) / r2_inside)[:, None])
            for axis in range(3):
                forces[:, axis] += (
                    np.bincount(a, weights=pair_force[:, axis], minlength=len(forces))
                    - np.bincount(b, weights=pair_force[:, axis], minlength=len(forces)))
            pair_bulk = 48 * epsilon * (22 * q12 - 5 * q6) / r2_inside
            bulk += float(np.sum(pair_bulk))
            if by_species:
                bulk_by_pair += np.bincount(kind[inside], weights=pair_bulk, minlength=3)
        for index, width in enumerate(widths):
            shell = (radius >= cutoff - width / 2) & (radius < cutoff + width / 2)
            counts[index] += np.bincount(kind[shell], minlength=3)
    if by_species:
        force2_by_atom = np.einsum("ij,ij->i", forces, forces)
        return (np.array([np.sum(force2_by_atom[frame.types == species])
                          for species in (1, 2)]), SPECIES_WEIGHTS @ bulk_by_pair, counts)
    return float(np.einsum("ij,ij->", forces, forces)), bulk, counts


def temperature(force2: float, bulk: float, counts: np.ndarray, width: float) -> tuple[float, float]:
    boundary = float(np.dot(counts, -2 * U_PRIME_CUTOFF / width))
    denominator = bulk + boundary
    return force2 / denominator, boundary  # k_B=1


def process_case(case: Case, args: argparse.Namespace) -> Result:
    widths = [SHELL_WIDTH]
    result = None
    for frame in selected_frames(case.trajectory, args.start_step, args.end_step, args.stride):
        if result is None:
            result = Result(case, len(frame.types),
                            (int(np.count_nonzero(frame.types == 1)),
                             int(np.count_nonzero(frame.types == 2))),
                            tuple(frame.bounds[:, 1] - frame.bounds[:, 0]),
                            0, frame.step, frame.step,
                            np.zeros(2),
                            np.zeros(2),
                            np.zeros((len(widths), 3), dtype=np.int64))
        force2, bulk, counts = frame_statistics(frame, widths, True)
        result.force2 += force2
        result.bulk += bulk
        result.shell_counts += counts
        result.frames += 1
        result.last_step = frame.step
    expected = (args.end_step - args.start_step) // args.stride + 1
    if (result is None or result.first_step != args.start_step or
            result.last_step != args.end_step or result.frames != expected):
        raise ValueError(f"Missing selected frames in {case.directory}; expected {expected}")
    return result


def fmt(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def write_outputs(results: list[Result], args: argparse.Namespace,
                  seeds: list[int]) -> list[Path]:
    tag = lambda x: format(x, ".12g").replace(".", "p")
    prefix = f"binlj_tconf_p{args.p}_a{tag(args.alpha)}_b{tag(args.beta)}"
    csv_path = args.output_dir / f"{prefix}.csv"
    json_path = args.output_dir / f"{prefix}_meta.json"
    if csv_path.exists() or json_path.exists():
        raise ValueError("Output exists; choose another --output-dir")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    by_key = {(r.case.method, r.case.seed): r for r in results}

    def components(r: Result, species: int | None):
        if species is None:
            return float(r.force2.sum()), float(r.bulk.sum()), r.shell_counts[0]
        slot = species - 1
        return (r.force2[slot], r.bulk[slot],
                r.shell_counts[0] * SPECIES_WEIGHTS[slot])

    def group_temperature(group: list[Result], species: int | None) -> float:
        parts = [components(r, species) for r in group]
        return temperature(sum(part[0] for part in parts),
                           sum(part[1] for part in parts),
                           np.sum([part[2] for part in parts], axis=0), SHELL_WIDTH)[0]

    rows = []
    for method in METHODS:
        for species in (None, 1, 2):
            label = "all" if species is None else str(species)
            for seed in seeds:
                tconf = group_temperature([by_key[method, seed]], species)
                cl_tconf = group_temperature([by_key["cl", seed]], species)
                rows.append([method, label, seed, fmt(tconf), fmt(tconf - cl_tconf)])
            tconf = group_temperature([by_key[method, seed] for seed in seeds], species)
            cl_tconf = group_temperature([by_key["cl", seed] for seed in seeds], species)
            rows.append([method, label, "mean", fmt(tconf), fmt(tconf - cl_tconf)])
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "species", "seed", "t_conf", "delta_t_conf_vs_cl"])
        writer.writerows(rows)
    record = {
        "units": "LJ reduced units, k_B=1",
        "formula": "T_conf=<sum_i |F_i^ref|^2>/<sum_i nabla_i^2 U_ref>; restrict sums for each species",
        "cutoff_term": "-2*u'(rc)*pair_count(rc-h/2 <= r < rc+h/2)/h",
        "averaging": "sum numerator and denominator separately, then take their ratio",
        "parameters": {"methods": METHODS, "seeds": seeds, "p": args.p, "alpha": args.alpha,
                       "beta": args.beta, "nsteps": args.nsteps,
                       "start_step": args.start_step, "end_step": args.end_step,
                       "stride": args.stride, "shell_width": SHELL_WIDTH,
                       "pairs": {name: {"epsilon": e, "sigma": s, "cutoff": c}
                                 for name, e, s, c in zip(PAIR_NAMES, EPSILON.tolist(),
                                                          SIGMA.tolist(), CUTOFF.tolist())}},
        "cases": [{"method": r.case.method, "seed": r.case.seed,
                   "run": r.case.directory.name,
                   "run_parameters": {key: r.case.meta[key] for key in
                                      ("PAIR_MODE", "SEED", "PAIR_SEED_BASE", "P", "MI_ALPHA",
                                       "MI_RHO", "PROD_STEPS", "RHO", "TEMP", "DT", "TDAMP")
                                      if key in r.case.meta},
                   "natoms": r.natoms, "type_counts_1_2": r.type_counts,
                   "box": r.box, "frames": r.frames,
                   "first_step": r.first_step, "last_step": r.last_step}
                  for r in results],
        "species": {"all": "mixture", "1": "A", "2": "B"},
        "pair_laplacian_and_cutoff_weights_11_12_22": SPECIES_WEIGHTS.tolist(),
    }
    json_path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n",
                         encoding="utf-8")
    return [csv_path, json_path]


def main(argv: list[str] | None = None) -> int:
    p = parser()
    args = p.parse_args(argv)
    try:
        seeds = [int(s) for s in re.split(r"[,\s]+", args.seeds.strip()) if s]
        methods = list(METHODS)
        cases = find_cases(args, methods, seeds)
        results = []
        for case in cases:
            result = process_case(case, args)
            results.append(result)
            print(f"{case.method} seed={case.seed}: {result.frames} frames", flush=True)
        reference = results[0]
        for result in results[1:]:
            if (result.natoms, result.box) != (reference.natoms, reference.box):
                raise ValueError(f"Different KA system in {result.case.directory}")
            for key in ("RHO", "TEMP", "DT", "TDAMP"):
                if float(result.case.meta[key]) != float(reference.case.meta[key]):
                    raise ValueError(f"Different {key} in {result.case.directory}")
        paths = write_outputs(results, args, seeds)
    except (ValueError, OSError) as exc:
        p.exit(2, f"error: {exc}\n")
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
