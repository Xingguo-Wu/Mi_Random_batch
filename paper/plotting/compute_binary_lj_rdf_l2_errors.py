#!/usr/bin/env python
from collections import defaultdict
from math import sqrt
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PAIRS = ("11", "22", "12")
METHODS = (
    ("mi-rbl_a0.0_rho0.95", "RBL"),
    ("vrrbl", "VR-RBL"),
    ("mi-rbl_a0.31_rho0.91", "Mi-RBL (beta=0.91)"),
    ("mi-rbl_a0.31_rho0.95", "Mi-RBL (beta=0.95)"),
    ("mi-rbl_a0.31_rho0.99", "Mi-RBL (beta=0.99)"),
)


def read_means(path: Path) -> dict[str, list[float]]:
    curves = defaultdict(list)
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        method, _, mean_g, _ = line.split()
        curves[method].append(float(mean_g))
    return curves


def rmse(values: list[float], reference: list[float]) -> float:
    return sqrt(sum((value - ref) ** 2 for value, ref in zip(values, reference)) / len(reference))


def main() -> None:
    curves = {pair: read_means(DATA_DIR / f"rdf_binary_lj_{pair}.txt") for pair in PAIRS}
    print("Method\tg_11 RMSE\tg_22 RMSE\tg_12 RMSE")
    for method, label in METHODS:
        errors = [rmse(curves[pair][method], curves[pair]["standard"]) for pair in PAIRS]
        print("\t".join((label, *(f"{error:.6g}" for error in errors))))


if __name__ == "__main__":
    main()
