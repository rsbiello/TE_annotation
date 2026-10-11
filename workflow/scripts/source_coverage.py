#!/usr/bin/env python3
"""Partition RepeatMasker coverage into mutually exclusive library sources."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SOURCE_ORDER = {
    "curated": 0,
    "de_novo": 1,
    "ambiguous_library": 2,
    "repeatmasker_internal": 3,
}

INTERSPERSED_PREFIXES = {
    "DNA",
    "LINE",
    "LTR",
    "Other",
    "PLE",
    "RC",
    "Retroposon",
    "SINE",
    "Unknown",
    "Unclassified",
    "Unspecified",
}


def fasta_size(path: str | Path) -> int:
    size = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                size += len("".join(line.split()))
    if size <= 0:
        raise ValueError(f"Genome FASTA contains no sequence: {path}")
    return size


def library_family_names(path: str | Path) -> set[str]:
    names: set[str] = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split(maxsplit=1)[0]
                names.add(identifier.split("#", 1)[0])
    return names


def family_source(
    family: str, curated_families: set[str], denovo_families: set[str]
) -> str:
    in_curated = family in curated_families
    in_denovo = family in denovo_families
    if in_curated and in_denovo:
        return "ambiguous_library"
    if in_denovo:
        return "de_novo"
    if in_curated:
        return "curated"
    return "repeatmasker_internal"


def is_interspersed(repeat_class: str) -> bool:
    return repeat_class.split("/", 1)[0] in INTERSPERSED_PREFIXES


def read_intervals(
    path: str | Path,
    curated_families: set[str],
    denovo_families: set[str],
) -> dict[str, list[tuple[str, int, int, str]]]:
    """Read non-suppressed RepeatMasker hits as zero-based half-open intervals."""
    intervals: dict[str, list[tuple[str, int, int, str]]] = {
        "all": [],
        "interspersed": [],
    }
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 15 or not fields[0].isdigit():
                continue
            if fields[-1] == "*":
                continue
            try:
                begin = int(fields[5])
                end = int(fields[6])
            except ValueError as error:
                raise ValueError(
                    f"Invalid RepeatMasker coordinates on line {line_number}"
                ) from error
            if begin > end:
                begin, end = end, begin
            family = fields[9]
            repeat_class = fields[10] or "Unknown"
            source = family_source(family, curated_families, denovo_families)
            interval = (fields[4], begin - 1, end, source)
            intervals["all"].append(interval)
            if is_interspersed(repeat_class):
                intervals["interspersed"].append(interval)
    return intervals


def source_sort_key(source: str) -> tuple[int, str]:
    return SOURCE_ORDER.get(source, len(SOURCE_ORDER)), source


def coverage_by_active_sources(
    intervals: list[tuple[str, int, int, str]],
) -> dict[str, int]:
    """Return disjoint coverage for each active-source combination."""
    events: dict[str, dict[int, list[tuple[str, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for sequence, start, end, source in intervals:
        if start < 0 or end <= start:
            raise ValueError(
                f"Invalid half-open interval: {sequence}:{start}-{end}"
            )
        events[sequence][start].append((source, 1))
        events[sequence][end].append((source, -1))

    coverage: dict[str, int] = defaultdict(int)
    for sequence_events in events.values():
        active: dict[str, int] = defaultdict(int)
        previous: int | None = None
        for position in sorted(sequence_events):
            if previous is not None and position > previous:
                sources = sorted(
                    (source for source, count in active.items() if count > 0),
                    key=source_sort_key,
                )
                if sources:
                    coverage["+".join(sources)] += position - previous
            for source, delta in sequence_events[position]:
                active[source] += delta
                if active[source] == 0:
                    del active[source]
                elif active[source] < 0:
                    raise ValueError(
                        f"Unbalanced interval events for source {source!r}"
                    )
            previous = position
        if active:
            raise ValueError("Unbalanced RepeatMasker interval events")
    return dict(coverage)


def write_report(
    intervals_by_scope: dict[str, list[tuple[str, int, int, str]]],
    genome_size: int,
    output_path: str | Path,
) -> None:
    rows: list[tuple[object, ...]] = []
    for scope in ("interspersed", "all"):
        coverage = coverage_by_active_sources(intervals_by_scope[scope])
        total = sum(coverage.values())
        for category in sorted(coverage):
            bp = coverage[category]
            rows.append(
                (
                    scope,
                    category,
                    bp,
                    f"{bp / genome_size * 100:.6f}",
                    f"{bp / total * 100:.6f}" if total else "0.000000",
                )
            )
        rows.append(
            (
                scope,
                "TOTAL_UNIQUE",
                total,
                f"{total / genome_size * 100:.6f}",
                "100.000000" if total else "0.000000",
            )
        )

    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "scope",
                "source_category",
                "masked_bp",
                "percent_genome",
                "percent_of_scope",
            )
        )
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeatmasker-out", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--curated-library", required=True)
    parser.add_argument("--denovo-library", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    curated = library_family_names(args.curated_library)
    denovo = library_family_names(args.denovo_library)
    intervals = read_intervals(args.repeatmasker_out, curated, denovo)
    write_report(intervals, fasta_size(args.genome), args.output)


if __name__ == "__main__":
    main()
