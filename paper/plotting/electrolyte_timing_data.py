
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"

DEFAULT_INPUTS = {
    n: DATA_DIR / f"electrolyte_timing_N{n}.tsv"
    for n in (3200, 6400, 12800, 25600, 51200)
}

TIMING_CASES = (
    ("elec_pppm", "PPPM", "#000000", "o"),
    ("elec_mirbe_a0p0_rho0p95", "RBE", "#1f77b4", "s"),
    ("elec_mirbe_a0p2_rho0p95", r"Mi-RBE ($\alpha=0.2$)", "#9467bd", "v"),
)


@dataclass(frozen=True)
class TimingPoint:
    n_atoms: int
    seconds_per_step: float


def read_timing_points(input_paths: dict[int, Path]) -> dict[str, list[TimingPoint]]:
    cases = {case for case, _, _, _ in TIMING_CASES}
    points_by_case: dict[str, list[TimingPoint]] = {case: [] for case, _, _, _ in TIMING_CASES}

    for n_atoms, path in input_paths.items():
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                case = row["method"]
                if case not in cases:
                    continue

                points_by_case[case].append(
                    TimingPoint(
                        n_atoms=n_atoms,
                        seconds_per_step=float(row["seconds_per_step"]),
                    )
                )

    return {case: sorted(points, key=lambda point: point.n_atoms) for case, points in points_by_case.items()}
