#!/usr/bin/env python3
"""Summarize RepeatMasker coverage per family and flag dominant families."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


HEADER = (
    "family",
    "class",
    "source",
    "masked_bp",
    "hit_count",
    "percent_genome",
    "percent_of_class",
    "status",
)


@dataclass
class FamilyStats:
    repeat_class: str
    source: str
    masked_bp: int = 0
    hit_count: int = 0


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


def repeatmasker_stats(
    path: str | Path, denovo_families: set[str]
) -> dict[str, FamilyStats]:
    stats: dict[str, FamilyStats] = {}
    current_sequence: str | None = None
    completed_sequences: set[str] = set()
    last_end_by_family: dict[str, int] = {}

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 15 or not fields[0].isdigit():
                continue
            if fields[-1] == "*":
                continue

            sequence = fields[4]
            try:
                begin = int(fields[5])
                end = int(fields[6])
            except ValueError as error:
                raise ValueError(
                    f"Invalid RepeatMasker coordinates on line {line_number}"
                ) from error
            if begin > end:
                begin, end = end, begin

            if sequence != current_sequence:
                if sequence in completed_sequences:
                    raise ValueError(
                        "RepeatMasker output is not grouped by query sequence; "
                        f"sequence {sequence!r} appears in multiple blocks"
                    )
                if current_sequence is not None:
                    completed_sequences.add(current_sequence)
                current_sequence = sequence
                last_end_by_family = {}

            family = fields[9]
            repeat_class = fields[10] or "Unknown"
            source = "de_novo" if family in denovo_families else "curated"
            record = stats.setdefault(
                family, FamilyStats(repeat_class=repeat_class, source=source)
            )
            if record.repeat_class != repeat_class:
                record.repeat_class = "Multiple"
            if source == "de_novo":
                record.source = source
            record.hit_count += 1

            previous_end = last_end_by_family.get(family, 0)
            if end > previous_end:
                record.masked_bp += end - max(begin, previous_end + 1) + 1
                last_end_by_family[family] = end
    return stats


def write_reports(
    stats: dict[str, FamilyStats],
    genome_size: int,
    summary_path: str | Path,
    flagged_path: str | Path,
    max_genome_percent: float,
    allowlist: set[str],
) -> list[str]:
    class_totals: dict[str, int] = {}
    for record in stats.values():
        class_totals[record.repeat_class] = (
            class_totals.get(record.repeat_class, 0) + record.masked_bp
        )

    rows: list[tuple[object, ...]] = []
    flagged_rows: list[tuple[object, ...]] = []
    flagged_names: list[str] = []
    for family, record in stats.items():
        percent_genome = record.masked_bp / genome_size * 100
        percent_class = (
            record.masked_bp / class_totals[record.repeat_class] * 100
            if class_totals[record.repeat_class]
            else 0.0
        )
        if family in allowlist:
            status = "allowed"
        elif record.source == "de_novo" and percent_genome > max_genome_percent:
            status = "flagged"
            flagged_names.append(family)
        else:
            status = "ok"
        row = (
            family,
            record.repeat_class,
            record.source,
            record.masked_bp,
            record.hit_count,
            f"{percent_genome:.6f}",
            f"{percent_class:.6f}",
            status,
        )
        rows.append(row)
        if status == "flagged":
            flagged_rows.append(row)

    rows.sort(key=lambda row: (-int(row[3]), str(row[0])))
    flagged_rows.sort(key=lambda row: (-int(row[3]), str(row[0])))
    for path, selected_rows in (
        (summary_path, rows),
        (flagged_path, flagged_rows),
    ):
        with Path(path).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(selected_rows)
    return flagged_names


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeatmasker-out", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--denovo-library", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--flagged", required=True)
    parser.add_argument("--max-genome-percent", type=float, default=5.0)
    parser.add_argument("--allow-family", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 < args.max_genome_percent <= 100:
        raise ValueError("--max-genome-percent must be in (0, 100]")
    genome_size = fasta_size(args.genome)
    denovo = library_family_names(args.denovo_library)
    stats = repeatmasker_stats(args.repeatmasker_out, denovo)
    flagged = write_reports(
        stats,
        genome_size,
        args.summary,
        args.flagged,
        args.max_genome_percent,
        set(args.allow_family),
    )
    print(
        f"[family_qc] Families: {len(stats)}, Flagged: {len(flagged)}, "
        f"Genome bp: {genome_size}"
    )


if __name__ == "__main__":
    main()
