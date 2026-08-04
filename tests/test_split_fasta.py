import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[1] / "workflow" / "scripts" / "split_fasta.py"
)
SPEC = importlib.util.spec_from_file_location("split_fasta", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class SplitFastaTests(unittest.TestCase):
    def test_records_are_preserved_and_balanced(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "genome.fa"
            fasta.write_text(
                ">chr1 description\n" + "A" * 100 + "\n"
                ">chr2\n" + "C" * 60 + "\n"
                ">scaffold_1\n" + "G" * 40 + "\n"
                ">scaffold_2\n" + "T" * 20 + "\n",
                encoding="utf-8",
            )

            output_dir = root / "chunks"
            manifest = root / "manifest.tsv"
            MODULE.write_chunks(fasta, output_dir, manifest, 2)

            chunk_fastas = sorted(output_dir.glob("chunk_*.fa"))
            self.assertEqual(len(chunk_fastas), 2)
            observed = "".join(path.read_text() for path in chunk_fastas)
            for sequence_id in ("chr1", "chr2", "scaffold_1", "scaffold_2"):
                self.assertEqual(observed.count(f">{sequence_id}"), 1)

            rows = manifest.read_text().strip().splitlines()
            self.assertEqual(rows[0], "chunk\tsequence_count\ttotal_bases\tfasta")
            loads = [int(row.split("\t")[2]) for row in rows[1:]]
            self.assertEqual(sorted(loads), [100, 120])

    def test_duplicate_sequence_ids_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "duplicate.fa"
            fasta.write_text(">chr1 first\nAAAA\n>chr1 second\nTTTT\n")
            with self.assertRaisesRegex(ValueError, "Duplicate FASTA sequence ID"):
                MODULE.scan_fasta(fasta)


if __name__ == "__main__":
    unittest.main()
