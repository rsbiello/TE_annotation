#!/usr/bin/env python3
"""Split a FASTA into balanced chunks without splitting sequence records."""

from __future__ import annotations

import argparse
import gzip
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class FastaRecord:
    index: int
    sequence_id: str
    length: int


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def scan_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    seen_ids: set[str] = set()
    header: str | None = None
    sequence_id: str | None = None
    length = 0

    def finish_record() -> None:
        nonlocal header, sequence_id, length
        if header is None or sequence_id is None:
            return
        if length == 0:
            raise ValueError(f"FASTA record {sequence_id!r} has no sequence")
        records.append(FastaRecord(len(records), sequence_id, length))

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if raw_line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                sequence_id = header.split()[0]
                if sequence_id in seen_ids:
                    raise ValueError(
                        f"Duplicate FASTA sequence ID {sequence_id!r}; "
                        "TRF BED coordinates require unique IDs"
                    )
                seen_ids.add(sequence_id)
                length = 0
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence data occurs before the first FASTA header "
                        f"at line {line_number}"
                    )
                length += len("".join(line.split()))

    finish_record()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def assign_balanced_chunks(
    records: list[FastaRecord], requested_chunks: int
) -> tuple[list[int], list[int], list[int]]:
    if requested_chunks < 1:
        raise ValueError("The number of TRF chunks must be at least 1")
    if requested_chunks > 256:
        raise ValueError("The number of TRF chunks must not exceed 256")

    chunk_count = min(requested_chunks, len(records))
    loads = [0] * chunk_count
    counts = [0] * chunk_count
    assignment = [0] * len(records)

    # Longest-processing-time bin packing keeps complete FASTA records intact
    # while balancing total sequence length among jobs.
    for record in sorted(records, key=lambda item: (-item.length, item.index)):
        chunk_index = min(range(chunk_count), key=lambda i: (loads[i], i))
        assignment[record.index] = chunk_index
        loads[chunk_index] += record.length
        counts[chunk_index] += 1

    return assignment, loads, counts


def write_chunks(
    input_fasta: Path,
    output_dir: Path,
    manifest: Path,
    requested_chunks: int,
) -> None:
    records = scan_fasta(input_fasta)
    assignment, loads, counts = assign_balanced_chunks(records, requested_chunks)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    chunk_names = [f"chunk_{i + 1:04d}" for i in range(len(loads))]
    chunk_paths = [output_dir / f"{name}.fa" for name in chunk_names]

    with ExitStack() as stack:
        outputs = [
            stack.enter_context(path.open("w", encoding="utf-8"))
            for path in chunk_paths
        ]
        record_index = -1
        destination: TextIO | None = None
        with open_text(input_fasta) as source:
            for raw_line in source:
                if raw_line.startswith(">"):
                    record_index += 1
                    destination = outputs[assignment[record_index]]
                if destination is None:
                    if raw_line.strip():
                        raise ValueError("Invalid FASTA encountered during second pass")
                    continue
                destination.write(raw_line)
                if not raw_line.endswith("\n"):
                    destination.write("\n")

    with manifest.open("w", encoding="utf-8") as handle:
        handle.write("chunk\tsequence_count\ttotal_bases\tfasta\n")
        for name, count, bases, path in zip(chunk_names, counts, loads, chunk_paths):
            handle.write(f"{name}\t{count}\t{bases}\t{path.resolve()}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--chunks", required=True, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    write_chunks(args.input, args.output_dir, args.manifest, args.chunks)


if "snakemake" in globals():
    write_chunks(
        Path(str(snakemake.input.genome)),
        Path(str(snakemake.output.chunks)),
        Path(str(snakemake.output.manifest)),
        int(snakemake.params.chunk_count),
    )
elif __name__ == "__main__":
    main()
