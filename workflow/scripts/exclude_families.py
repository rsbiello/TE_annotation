#!/usr/bin/env python3
"""Remove explicitly reviewed repeat families from a FASTA library."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator, TextIO


def fasta_records(path: str | Path) -> Iterator[tuple[str, str, str]]:
    header: str | None = None
    sequence: list[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield family_name(header), header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError("FASTA sequence encountered before the first header")
                sequence.append(line)
    if header is not None:
        yield family_name(header), header, "".join(sequence)


def family_name(header: str) -> str:
    """Return the RepeatMasker family name without the ``#class`` suffix."""

    return header.split(maxsplit=1)[0].split("#", 1)[0]


def write_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def exclude_families(
    input_fasta: str | Path,
    output_fasta: str | Path,
    report_tsv: str | Path,
    exclusions: list[str],
) -> tuple[int, int]:
    requested = set(exclusions)
    if len(requested) != len(exclusions):
        raise ValueError("Duplicate names were supplied in the family exclusion list")

    seen: set[str] = set()
    kept = 0
    removed = 0
    with Path(output_fasta).open("w", encoding="utf-8") as output_handle, Path(
        report_tsv
    ).open("w", encoding="utf-8") as report_handle:
        report_handle.write("family\theader\n")
        for name, header, sequence in fasta_records(input_fasta):
            if name in requested:
                report_handle.write(f"{name}\t{header}\n")
                seen.add(name)
                removed += 1
                continue
            write_record(output_handle, header, sequence)
            kept += 1

    missing = requested - seen
    if missing:
        Path(output_fasta).unlink(missing_ok=True)
        Path(report_tsv).unlink(missing_ok=True)
        names = ", ".join(sorted(missing))
        raise ValueError(f"Excluded family name(s) not found in the input library: {names}")
    if kept == 0:
        Path(output_fasta).unlink(missing_ok=True)
        Path(report_tsv).unlink(missing_ok=True)
        raise ValueError("Family exclusions removed every sequence from the library")
    return kept, removed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input RepeatModeler FASTA")
    parser.add_argument("--output", required=True, help="Filtered output FASTA")
    parser.add_argument("--report", required=True, help="Removed-family TSV report")
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help="Family name before the #class suffix; may be repeated",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    kept, removed = exclude_families(
        args.input, args.output, args.report, args.exclude
    )
    print(f"[exclude_families] Kept: {kept}, Removed: {removed}")


if __name__ == "__main__":
    main()
