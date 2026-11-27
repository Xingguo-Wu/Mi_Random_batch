#!/usr/bin/env python
from collections import defaultdict
from math import sqrt
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent.parent / "data"
ALPHAS = {5: (0.4, 0.5, 0.6), 10: (0.4, 0.5, 0.6), 20: (0.1, 0.2, 0.3)}


def read_means(path: Path) -> dict[str, list[float]]:
    curves = defaultdict(list)
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        case, _, mean_density = line.split()
        curves[case].append(float(mean_density))
    return curves


def rmse(values: list[float], reference: list[float]) -> float:
    return sqrt(sum((value - ref) ** 2 for value, ref in zip(values, reference)) / len(reference))


def main() -> None:
    print("P\tMethod\tRMSE")
    for p, alphas in ALPHAS.items():
        curves = read_means(DATA_DIR / f"electrolyte_densityP{p}.dat")
        reference = curves["pppm"]
        for alpha in (0.0, *alphas):
            case = f"mirbe_a{alpha:.1f}_rho0.95"
            label = "RBE" if alpha == 0.0 else f"Mi-RBE (alpha={alpha:.1f})"
            print(f"{p}\t{label}\t{rmse(curves[case], reference):.6g}")


if __name__ == "__main__":
    main()
