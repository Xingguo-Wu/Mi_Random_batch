#!/usr/bin/env python3
"""Summarize production seconds per step for the electrolyte timing panel."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from pathlib import Path


LOOP_RE = re.compile(
    r"Loop time of\s+([0-9.eE+-]+)\s+on\s+\d+\s+procs\s+for\s+(\d+)\s+steps\s+with\s+\d+\s+atoms"
)
SEED_RE = re.compile(r"(?:^|[_-])(?:s|seed)(\d+)(?:$|[_-])", re.IGNORECASE)
REMOVE_SEED_RE = re.compile(r"(?:^|[_-])(?:s|seed)\d+(?=$|[_-])", re.IGNORECASE)
REMOVE_NSTEPS_RE = re.compile(r"(?:^|[_-])n\d+(?=$|[_-])", re.IGNORECASE)
PAPER_METHODS = (
    "elec_mirbe_a0p0_rho0p95",
    "elec_mirbe_a0p2_rho0p95",
    "elec_pppm",
)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--log-glob", default="**/log.lammps")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--summary-name", default="electrolyte_timing_method_summary.tsv")
    parser.add_argument("--per-seed-name", default="electrolyte_timing_seed_stats.tsv")
    return parser.parse_args()


def case_and_method(path: Path) -> tuple[str, str]:
    case = path.parent.name if path.name == "log.lammps" else path.stem
    method = REMOVE_NSTEPS_RE.sub("_", REMOVE_SEED_RE.sub("_", case))
    method = re.sub(r"[_-]+", "_", method).strip("_")
    return case, method


def production_seconds_per_step(path: Path) -> float:
    value = None
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = LOOP_RE.search(line)
            if match:
                value = float(match.group(1)) / int(match.group(2))
    if value is None:
        raise ValueError(f"No LAMMPS Loop time in {path}")
    return value


def main():
    args = arguments()
    rows = []
    for path in sorted(args.input_dir.glob(args.log_glob)):
        case, method = case_and_method(path)
        if method in PAPER_METHODS:
            seed_match = SEED_RE.search(case)
            seed = seed_match.group(1) if seed_match else "unknown"
            rows.append((method, seed, production_seconds_per_step(path)))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    seed_path = args.output_dir / args.per_seed_name
    with seed_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("method", "seed", "seconds_per_step"))
        writer.writerows((method, seed, f"{value:.10e}") for method, seed, value in rows)

    grouped = {method: [] for method in PAPER_METHODS}
    for method, _, value in rows:
        grouped[method].append(value)
    summary_path = args.output_dir / args.summary_name
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("method", "nseeds", "seconds_per_step"))
        writer.writerows(
            (method, len(values), f"{statistics.fmean(values):.10e}")
            for method, values in grouped.items() if values
        )
    print(f"Wrote {seed_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
