#!/usr/bin/env python3
"""Patch RepeatMasker's divergence utility for reproducible landscapes.

RepeatMasker 4.2.4 reuses Kimura values embedded in modern .align files even
when -noCpGMod is requested. It also writes negative divergence values to bins
that createRepeatLandscape.pl never emits. This script patches a job-local copy
of the installed utility; it never changes the Conda environment.
"""

from __future__ import annotations

import argparse
from pathlib import Path


ORIGINAL_BLOCK = """    # Obtain divergence from modern *.align files directly
    $div = $result->getPctKimuraDiverge();

    if ( $div eq "" ) {

      # Calculate divergence on the fly
      ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
          $result->calcKimuraDivergence( divCpGMod => 1 );
      $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    }
"""

PATCHED_BLOCK = """    # Force a genuine standard-Kimura recalculation when requested.
    # Modern .align files already contain a CpG-adjusted value, so merely
    # changing the report heading is insufficient.
    if ( $options{'noCpGMod'} ) {
      ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
          $result->calcKimuraDivergence( divCpGMod => 0 );
      $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    }
    else {
      # Reuse the embedded CpG-adjusted value when it is available.
      $div = $result->getPctKimuraDiverge();
      if ( $div eq "" ) {
        ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
            $result->calcKimuraDivergence( divCpGMod => 1 );
        $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
      }
    }
"""

OVERLAP_MARKER = """  if ( $prevQueryName eq $seqName ) {
"""

CLAMP_BLOCK = """  # A negative corrected divergence has no biological age meaning and
  # would otherwise be stored in a bin omitted from the 0..70 landscape.
  # Retain its coverage in the zero-divergence bin.
  $div = 0 if ( $div < 0 );

  if ( $prevQueryName eq $seqName ) {
"""


def patch_source(source: str) -> str:
    if source.count(ORIGINAL_BLOCK) != 1:
        raise ValueError(
            "Unsupported calcDivergenceFromAlign.pl layout: expected the "
            "RepeatMasker 4.2.4 divergence block exactly once"
        )
    patched = source.replace(ORIGINAL_BLOCK, PATCHED_BLOCK, 1)
    if patched.count(OVERLAP_MARKER) != 1:
        raise ValueError(
            "Unsupported calcDivergenceFromAlign.pl layout: expected the "
            "overlap-processing marker exactly once"
        )
    return patched.replace(OVERLAP_MARKER, CLAMP_BLOCK, 1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = Path(args.input).read_text(encoding="utf-8")
    Path(args.output).write_text(patch_source(source), encoding="utf-8")


if __name__ == "__main__":
    main()
