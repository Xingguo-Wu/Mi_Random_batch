#!/usr/bin/env python3
import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path


PAIRS = (("11", 2), ("12", 4), ("22", 6))


@dataclass
class SeedResult:
    radius: list[float]
    rdf: dict[str, list[float]]


def has_token(name: str, token: str) -> bool:
    return f"_{token}_" in f"{name}_"


def matches(name: str, method: str, prefix: str, seed: str, nsteps: int) -> bool:
    if not name.startswith(prefix) or not has_token(name, f"n{nsteps}"):
        return False
    if method == "standard":
        return has_token(name, f"s{seed}")
    if not has_token(name, f"ps{seed}"):
        return False
    if method == "vrrbl":
        return True
    alpha, rho = re.search(r"a([\d.]+)_rho([\d.]+)", method).groups()
    alpha_tag = alpha.replace(".", "p")
    rho_tag = rho.replace(".", "p")
    return has_token(name, f"a{alpha_tag}") and (
        float(alpha) == 0.0 or has_token(name, f"rho{rho_tag}"))


def find_file(root: Path, method: str, prefix: str, seed: str, nsteps: int,
              rdf_name: str) -> Path:
    paths = [folder / rdf_name for folder in sorted(root.iterdir())
             if folder.is_dir() and matches(folder.name, method, prefix, seed, nsteps)
             and (folder / rdf_name).is_file()]
    return paths[0]


def read_blocks(path: Path):
    step = None
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 2:
            if step is not None:
                yield step, rows
            step, rows = int(parts[0]), []
        elif step is not None:
            rows.append([float(value) for value in parts])
    if step is not None:
        yield step, rows


def average_seed(path: Path, start: int, end: int, rdf_length: float) -> SeedResult:
    selected = [rows for step, rows in read_blocks(path) if start <= step <= end]
    template = [row for row in selected[0] if row[1] <= rdf_length]
    count = len(template)
    sums = {pair: [0.0] * count for pair, _ in PAIRS}
    for rows in selected:
        for index, row in enumerate(rows[:count]):
            for pair, rdf_col in PAIRS:
                sums[pair][index] += row[rdf_col]
    frames = len(selected)
    return SeedResult([row[1] for row in template],
                      {pair: [value / frames for value in values]
                       for pair, values in sums.items()})


def mean_sem(series: list[list[float]]):
    n = len(series)
    mean = [sum(values) / n for values in zip(*series)]
    std = [math.sqrt(sum((value - avg) ** 2 for value in values) / (n - 1))
           for avg, values in zip(mean, zip(*series))]
    sem = [value / math.sqrt(n) for value in std]
    return mean, sem


def write_stats(output_dir: Path, pair: str, results: dict[str, list[SeedResult]]) -> None:
    lines = ["# method r mean_g sem_g"]
    for method, seeds in results.items():
        mean, sem = mean_sem([seed.rdf[pair] for seed in seeds])
        for radius, rdf, se in zip(seeds[0].radius, mean, sem):
            values = (radius, rdf, se)
            lines.append("{} {}".format(
                method, " ".join(format(value, ".12g") for value in values)))
    (output_dir / f"rdf_binary_lj_{pair}.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-step", type=int, required=True)
    parser.add_argument("--end-step", type=int, required=True)
    parser.add_argument("--rdf-length", type=float, required=True)
    parser.add_argument("--method", action="append", required=True, metavar="NAME=PREFIX")
    parser.add_argument("--seeds", nargs="+", required=True)
    parser.add_argument("--nsteps", type=int, default=400000)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rdf-name", default="rdf.binary_lj_ka.dat")
    args = parser.parse_args()

    results = {}
    for raw in args.method:
        method, prefix = raw.split("=", 1)
        results[method] = [average_seed(
            find_file(args.input_dir, method, prefix, seed, args.nsteps, args.rdf_name),
            args.start_step, args.end_step, args.rdf_length) for seed in args.seeds]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for pair, _ in PAIRS:
        write_stats(args.output_dir, pair, results)


if __name__ == "__main__":
    main()
