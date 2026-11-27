#!/usr/bin/env python3
import argparse
import re
import statistics
from pathlib import Path


def read_kinetic(path: Path) -> list[tuple[int, float]]:
    records = []
    kinetic_column = None
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if parts and parts[0] == "Step":
            kinetic_column = parts.index("KinEng") if "KinEng" in parts else None
            continue
        if not parts:
            kinetic_column = None
            continue
        if kinetic_column is None:
            continue
        try:
            step = int(parts[0])
        except ValueError:
            kinetic_column = None
            continue
        records.append((step, float(parts[kinetic_column])))
    return records


def summarize(records: list[tuple[int, float]], start: int, end: int):
    selected = [(step, kinetic) for step, kinetic in records if start <= step <= end]
    values = [kinetic for _, kinetic in selected]
    return (len(values), selected[0][0], selected[-1][0],
            statistics.fmean(values), statistics.pstdev(values))


def short_case_name(name: str) -> str:
    match = re.fullmatch(r"binlj_standard_s(\d+)_n\d+", name)
    if match:
        return f"cl_s{match[1]}"
    match = re.fullmatch(r"binlj_vrrbl_ps(\d+)_n\d+", name)
    if match:
        return f"vrrbl_s{match[1]}"
    match = re.fullmatch(r"binlj_mirbl_cmp_ps(\d+)_a(\d+p\d+)_rho(\d+p\d+)_n\d+", name)
    if match:
        seed, alpha, beta = match.groups()
        method = "rbl" if float(alpha.replace("p", ".")) == 0 else "mirbl"
        return f"{method}_a{alpha}_b{beta}_s{seed}"
    raise ValueError(f"unrecognized binary LJ case: {name}")


def windows(args) -> list[tuple[int, int]]:
    if args.window_list:
        return [tuple(map(int, item.split(":")))
                for item in args.window_list.split(",")]
    return [(args.start_step, args.end_step)]


def per_log_report(case: str, total: int, summaries) -> str:
    lines = [f"log_file: {case}/log.lammps", f"records_total: {total}"]
    for start, end, used, first, last, mean, std in summaries:
        lines.extend(["", f"window: step range [{start}, {end}]",
                      f"records_used: {used}", f"selected_step_start: {first}",
                      f"selected_step_end: {last}",
                      f"{'metric':<14}{'mean':>16}{'std':>16}",
                      f"{'kinetic':<14}{mean:16.8e}{std:16.8e}"])
    return "\n".join(lines) + "\n"


def combined_report(rows) -> str:
    header = (f"{'case':<30}{'metric':<14}{'used':>8}{'step_start':>12}"
              f"{'step_end':>12}{'mean':>16}{'std':>16}")
    lines = ["combined_summary_for_multiple_logs", "", header, "-" * len(header)]
    for case, used, first, last, mean, std in rows:
        lines.append(f"{case:<30}{'kinetic':<14}{used:>8}{first:>12}{last:>12}"
                     f"{mean:16.8e}{std:16.8e}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-item", action="append", required=True,
                        help="LOG_PATH or LOG_PATH=>SUMMARY_PATH; repeat for each run")
    parser.add_argument("--start-step", type=int)
    parser.add_argument("--end-step", type=int)
    parser.add_argument("--window-list", help="Comma-separated START:END windows")
    parser.add_argument("--summary-output", type=Path)
    args = parser.parse_args()

    window_list = windows(args)
    combined_rows = []
    for item in args.log_item:
        log_text, separator, output_text = item.partition("=>")
        log_path = Path(log_text)
        output_path = Path(output_text) if separator else log_path.parent / "log_stats_summary.txt"
        records = read_kinetic(log_path)
        case = short_case_name(log_path.parent.name)
        summaries = []
        for start, end in window_list:
            used, first, last, mean, std = summarize(records, start, end)
            summaries.append((start, end, used, first, last, mean, std))
            combined_rows.append((case, used, first, last, mean, std))
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(per_log_report(case, len(records), summaries), encoding="utf-8")
        print(output_path)

    if args.summary_output:
        args.summary_output.parent.mkdir(parents=True, exist_ok=True)
        args.summary_output.write_text(combined_report(combined_rows), encoding="utf-8")
        print(args.summary_output)


if __name__ == "__main__":
    main()
