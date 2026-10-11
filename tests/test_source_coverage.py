import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "source_coverage.py"
SPEC = importlib.util.spec_from_file_location("source_coverage", SCRIPT)
SOURCE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = SOURCE
SPEC.loader.exec_module(SOURCE)


class SourceCoverageTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.genome = self.root / "genome.fa"
        self.genome.write_text(">chr1\n" + "A" * 200 + "\n")
        self.curated = self.root / "curated.fa"
        self.curated.write_text(">curated#DNA/hAT\nACGT\n")
        self.denovo = self.root / "denovo.fa"
        self.denovo.write_text(">denovo#LINE/CR1\nACGT\n")
        self.rm_out = self.root / "genome.fa.out"
        self.rm_out.write_text(
            "header\nheader\nheader\n"
            " 100 1.0 0.0 0.0 chr1 1 100 (100) + denovo LINE/CR1 1 100 (0) 1\n"
            " 100 1.0 0.0 0.0 chr1 81 140 (60) + curated DNA/hAT 1 60 (0) 2\n"
            " 100 1.0 0.0 0.0 chr1 151 160 (40) + (AT)n Simple_repeat 1 10 (0) 3\n"
            " 100 1.0 0.0 0.0 chr1 170 180 (20) + ignored DNA/hAT 1 11 (0) 4 *\n"
        )

    def tearDown(self):
        self.temporary.cleanup()

    def test_disjoint_source_partition(self):
        intervals = SOURCE.read_intervals(
            self.rm_out,
            SOURCE.library_family_names(self.curated),
            SOURCE.library_family_names(self.denovo),
        )
        interspersed = SOURCE.coverage_by_active_sources(intervals["interspersed"])
        self.assertEqual(
            interspersed,
            {
                "de_novo": 80,
                "curated+de_novo": 20,
                "curated": 40,
            },
        )
        all_coverage = SOURCE.coverage_by_active_sources(intervals["all"])
        self.assertEqual(all_coverage["repeatmasker_internal"], 10)
        self.assertEqual(sum(all_coverage.values()), 150)

    def test_report_totals_are_unique(self):
        intervals = SOURCE.read_intervals(
            self.rm_out,
            SOURCE.library_family_names(self.curated),
            SOURCE.library_family_names(self.denovo),
        )
        report = self.root / "source_coverage.tsv"
        SOURCE.write_report(intervals, SOURCE.fasta_size(self.genome), report)
        with report.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        totals = {
            row["scope"]: row
            for row in rows
            if row["source_category"] == "TOTAL_UNIQUE"
        }
        self.assertEqual(totals["interspersed"]["masked_bp"], "140")
        self.assertEqual(totals["interspersed"]["percent_genome"], "70.000000")
        self.assertEqual(totals["all"]["masked_bp"], "150")


if __name__ == "__main__":
    unittest.main()
