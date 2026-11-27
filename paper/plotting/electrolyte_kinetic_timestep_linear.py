"""Read electrolyte timestep data and fit the kinetic-error curves in the paper figure."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


TARGET_KINETIC = 1.50

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"

DEFAULT_INPUTS = {
    timestep: DATA_DIR / f"electrolyte_kinetic_dt{str(timestep).replace('.', 'p')}_summary.tsv"
    for timestep in (0.001, 0.0025, 0.005, 0.0075, 0.01)
}

@dataclass(frozen=True)
class KineticPoint:
    timestep: float
    kinetic_mean: float

    @property
    def kinetic_error(self) -> float:
        return self.kinetic_mean - TARGET_KINETIC


@dataclass(frozen=True)
class LinearFit:
    slope: float
    intercept: float

    def predict(self, x_value: float) -> float:
        return self.slope * x_value + self.intercept


def read_kinetic_points(input_paths: dict[float, Path], cases) -> dict[float, list[KineticPoint]]:
    case_config = {case: alpha for case, _, alpha, _, _ in cases}
    points_by_alpha: dict[float, list[KineticPoint]] = {alpha: [] for _, _, alpha, _, _ in cases}

    for timestep, path in input_paths.items():
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                case = row["case"]
                if case not in case_config:
                    continue
                points_by_alpha[case_config[case]].append(
                    KineticPoint(
                        timestep=float(timestep),
                        kinetic_mean=float(row["kin_mean"]),
                    )
                )

    return {alpha: sorted(points, key=lambda point: point.timestep) for alpha, points in points_by_alpha.items()}


def fit_line(points: list[KineticPoint]) -> LinearFit:
    x_values = [point.timestep for point in points]
    y_values = [point.kinetic_error for point in points]
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)

    ss_xx = sum((x - x_mean) ** 2 for x in x_values)
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values)) / ss_xx
    intercept = y_mean - slope * x_mean
    return LinearFit(slope=slope, intercept=intercept)


def linspace_between(x_min: float, x_max: float, count: int = 100) -> list[float]:
    return [x_min + (x_max - x_min) * index / (count - 1) for index in range(count)]
