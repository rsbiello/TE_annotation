import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "family_qc.py"
SPEC = importlib.util.spec_from_file_location("family_qc", SCRIPT)
QC = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = QC
SPEC.loader.exec_module(QC)


class FamilyQCTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.genome = self.root / "genome.fa"
        self.genome.write_text(">chr1\n" + "A" * 200 + "\n")
        self.library = self.root / "denovo.fa"
        self.library.write_text(">famA#LINE/CR1\nACGT\n")
        self.rm_out = self.root / "genome.fa.out"
        self.rm_out.write_text(
            "header\nheader\nheader\n"
            " 100 1.0 0.0 0.0 chr1 1 60 (140) + famA LINE/CR1 1 60 (0) 1\n"
            " 100 1.0 0.0 0.0 chr1 50 100 (100) + famA LINE/CR1 1 51 (9) 2\n"
            " 100 1.0 0.0 0.0 chr1 20 40 (160) + ignored DNA/hAT 1 21 (0) 3 *\n"
            " 100 1.0 0.0 0.0 chr1 150 169 (31) + curated DNA/hAT 1 20 (0) 4\n"
        )

    def tearDown(self):
        self.temporary.cleanup()

    def test_union_coverage_source_and_flagging(self):
        stats = QC.repeatmasker_stats(
            self.rm_out, QC.library_family_names(self.library)
        )
        self.assertEqual(stats["famA"].masked_bp, 100)
        self.assertEqual(stats["famA"].hit_count, 2)
        self.assertEqual(stats["famA"].source, "de_novo")
        self.assertNotIn("ignored", stats)

        summary = self.root / "summary.tsv"
        flagged = self.root / "flagged.tsv"
        names = QC.write_reports(
            stats,
            QC.fasta_size(self.genome),
            summary,
            flagged,
            30.0,
            set(),
        )
        self.assertEqual(names, ["famA"])
        with flagged.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["family"], "famA")
        self.assertEqual(rows[0]["percent_genome"], "50.000000")

    def test_allowlist_preserves_validated_expansion(self):
        stats = QC.repeatmasker_stats(
            self.rm_out, QC.library_family_names(self.library)
        )
        summary = self.root / "summary.tsv"
        flagged = self.root / "flagged.tsv"
        names = QC.write_reports(
            stats,
            QC.fasta_size(self.genome),
            summary,
            flagged,
            30.0,
            {"famA"},
        )
        self.assertEqual(names, [])
        self.assertEqual(len(flagged.read_text().splitlines()), 1)
        self.assertIn("allowed", summary.read_text())

    def test_large_curated_family_is_reported_but_not_flagged(self):
        stats = QC.repeatmasker_stats(
            self.rm_out, QC.library_family_names(self.library)
        )
        summary = self.root / "summary.tsv"
        flagged = self.root / "flagged.tsv"
        names = QC.write_reports(
            {"curated": stats["curated"]},
            QC.fasta_size(self.genome),
            summary,
            flagged,
            5.0,
            set(),
        )
        self.assertEqual(names, [])
        self.assertIn("curated", summary.read_text())
        self.assertEqual(len(flagged.read_text().splitlines()), 1)


if __name__ == "__main__":
    unittest.main()
