#!/usr/bin/env python3
import argparse
import csv
import statistics
from pathlib import Path


def production_kinetic_mean(path: Path) -> float:
    values = []
    kinetic_column = None
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            fields = line.split()
            if fields and fields[0] == "Step":
                kinetic_column = next(
                    (fields.index(name) for name in ("KinEng", "KE") if name in fields),
                    None,
                )
                values = []
            elif kinetic_column is not None:
                try:
                    int(fields[0])
                    values.append(float(fields[kinetic_column]))
                except (IndexError, ValueError):
                    kinetic_column = None
    return statistics.fmean(values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--log-glob", default="**/log.lammps")
    parser.add_argument("--output", type=Path, default=Path("kinetic_summary.tsv"))
    args = parser.parse_args()

    seed_means = [production_kinetic_mean(path)
                  for path in sorted(args.input_dir.glob(args.log_glob))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("nseeds", "kin_mean"))
        writer.writerow((len(seed_means), f"{statistics.fmean(seed_means):.10e}"))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
