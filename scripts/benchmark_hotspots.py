#!/usr/bin/env python3

import argparse
import datetime as dt
import re
import statistics
import subprocess
import sys
from dataclasses import dataclass
from typing import Sequence


TIME_PATTERN = re.compile(r"in ([0-9]+(?:\.[0-9]+)?) seconds")


@dataclass
class BenchmarkCase:
    name: str
    command: list[str]


def parse_elapsed_seconds(cprofile_output: str) -> float:
    match = TIME_PATTERN.search(cprofile_output)
    if match is None:
        raise RuntimeError("Could not parse elapsed time from cProfile output.")
    return float(match.group(1))


def run_case(case: BenchmarkCase, iterations: int, warmups: int) -> list[float]:
    timings: list[float] = []
    for idx in range(warmups + iterations):
        proc = subprocess.run(
            case.command,
            check=True,
            capture_output=True,
            text=True,
        )
        elapsed = parse_elapsed_seconds(proc.stdout)
        if idx >= warmups:
            timings.append(elapsed)
    return timings


def summarize(timings: Sequence[float]) -> tuple[float, float, float, float]:
    return (
        statistics.mean(timings),
        statistics.median(timings),
        min(timings),
        max(timings),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run repeat cProfile benchmarks for hotspot BioKIT commands.",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=5,
        help="Number of profiled runs to aggregate per command (default: 5).",
    )
    parser.add_argument(
        "--warmups",
        type=int,
        default=1,
        help="Warmup runs per command before recording timings (default: 1).",
    )
    parser.add_argument(
        "--python",
        type=str,
        default=sys.executable,
        help="Python interpreter used for cProfile command runs.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="",
        help="Optional TSV output path for benchmark summary rows.",
    )
    args = parser.parse_args()

    if args.iterations < 1:
        raise ValueError("iterations must be >= 1")
    if args.warmups < 0:
        raise ValueError("warmups must be >= 0")

    python_bin = args.python
    cases = [
        BenchmarkCase(
            name="genome_assembly_metrics",
            command=[
                python_bin,
                "-m",
                "cProfile",
                "-s",
                "cumtime",
                "biokit-runner.py",
                "genome_assembly_metrics",
                "tests/sample_files/GCF_000146045.2_R64_genomic.fna",
            ],
        ),
        BenchmarkCase(
            name="character_frequency",
            command=[
                python_bin,
                "-m",
                "cProfile",
                "-s",
                "cumtime",
                "biokit-runner.py",
                "character_frequency",
                "tests/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
            ],
        ),
        BenchmarkCase(
            name="relative_synonymous_codon_usage",
            command=[
                python_bin,
                "-m",
                "cProfile",
                "-s",
                "cumtime",
                "biokit-runner.py",
                "relative_synonymous_codon_usage",
                "tests/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
            ],
        ),
        BenchmarkCase(
            name="translate_sequence",
            command=[
                python_bin,
                "-m",
                "cProfile",
                "-s",
                "cumtime",
                "biokit-runner.py",
                "translate_sequence",
                "tests/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
            ],
        ),
    ]

    print(
        f"Benchmarking with python={python_bin}, iterations={args.iterations}, warmups={args.warmups}"
    )
    print("")
    print("Command\tMean(s)\tMedian(s)\tMin(s)\tMax(s)")
    rows: list[str] = []
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat()

    for case in cases:
        timings = run_case(case, args.iterations, args.warmups)
        mean_v, median_v, min_v, max_v = summarize(timings)
        print(
            f"{case.name}\t{mean_v:.3f}\t{median_v:.3f}\t{min_v:.3f}\t{max_v:.3f}"
        )
        rows.append(
            "\t".join(
                [
                    timestamp,
                    case.name,
                    f"{mean_v:.6f}",
                    f"{median_v:.6f}",
                    f"{min_v:.6f}",
                    f"{max_v:.6f}",
                    str(args.iterations),
                    str(args.warmups),
                    python_bin,
                ]
            )
        )

    if args.output:
        with open(args.output, "w") as out:
            out.write(
                "timestamp_utc\tcommand\tmean_seconds\tmedian_seconds\tmin_seconds\tmax_seconds\titerations\twarmups\tpython\n"
            )
            out.write("\n".join(rows))
            out.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
