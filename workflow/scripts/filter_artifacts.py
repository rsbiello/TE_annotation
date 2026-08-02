#!/usr/bin/env python3
"""Filter RepeatModeler families with likely host-protein BLASTx hits.

Families without a BLASTx hit are retained. Families whose best hit contains a
recognized transposable-element term are also retained. Other protein-matching
families are written to the removal report.
"""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path
from typing import Iterator, TextIO


TE_TERMS = (
    "transposase",
    "reverse transcriptase",
    "retrotransposon",
    "transposable element",
    "integrase",
    "gag",
    "pol",
    "env",
    "line",
    "sine",
    "ltr",
    "gypsy",
    "copia",
    "helitron",
    "mariner",
    "tc1",
    "mutator",
    "harbinger",
    "cacta",
    "hat",
    "rte",
    "penelope",
    "piggybac",
)
TE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:" + "|".join(re.escape(term) for term in TE_TERMS) + r")(?![A-Za-z0-9])",
    re.IGNORECASE,
)


def parse_blast(blast_file: str | Path) -> dict[str, tuple[float, str]]:
    """Return the lowest-e-value BLASTx hit for each query.

    The expected input is the 13-column tabular format emitted by the
    workflow. Blank and malformed rows are ignored.
    """

    best: dict[str, tuple[float, str]] = {}
    with Path(blast_file).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t", 12)
            if len(columns) != 13:
                continue

            query_id = columns[0]
            try:
                evalue = float(columns[10])
            except ValueError:
                continue
            subject_title = columns[12].replace("\t", " ")

            if query_id not in best or evalue < best[query_id][0]:
                best[query_id] = (evalue, subject_title)
    return best


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(sequence_id, full_header, sequence)`` FASTA records."""

    header: str | None = None
    sequence: list[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError("FASTA sequence encountered before the first header")
                sequence.append(line)
    if header is not None:
        yield header.split()[0], header, "".join(sequence)


def write_fasta(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for line in textwrap.wrap(sequence, width=80) or [""]:
        handle.write(f"{line}\n")


def is_te_protein(subject_title: str) -> bool:
    return TE_PATTERN.search(subject_title) is not None


def filter_families(
    input_fasta: str | Path,
    blast_file: str | Path,
    output_fasta: str | Path,
    report_file: str | Path,
) -> tuple[int, int]:
    best_hits = parse_blast(blast_file)
    kept = 0
    removed = 0

    with Path(output_fasta).open("w", encoding="utf-8") as output_handle, Path(
        report_file
    ).open("w", encoding="utf-8") as report_handle:
        report_handle.write("family\tevalue\tprotein_hit\n")
        for sequence_id, header, sequence in parse_fasta(input_fasta):
            best_hit = best_hits.get(sequence_id)
            if best_hit is None or is_te_protein(best_hit[1]):
                write_fasta(output_handle, header, sequence)
                kept += 1
                continue

            evalue, subject_title = best_hit
            report_handle.write(f"{sequence_id}\t{evalue:g}\t{subject_title}\n")
            removed += 1

    return kept, removed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Classified family FASTA")
    parser.add_argument("--blast", required=True, help="13-column BLASTx TSV")
    parser.add_argument("--output", required=True, help="Filtered family FASTA")
    parser.add_argument("--report", required=True, help="Removed-family TSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    kept, removed = filter_families(args.input, args.blast, args.output, args.report)
    print(f"[filter_artifacts] Kept: {kept}, Removed: {removed}")


if __name__ == "__main__":
    main()
