import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "filter_artifacts.py"
SPEC = importlib.util.spec_from_file_location("filter_artifacts", SCRIPT)
FILTER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(FILTER)


def blast_row(query, evalue, title):
    columns = [
        query,
        "subject",
        "90.0",
        "100",
        "0",
        "0",
        "1",
        "100",
        "1",
        "100",
        str(evalue),
        "200",
        title,
    ]
    return "\t".join(columns) + "\n"


class FilterArtifactsTests(unittest.TestCase):
    def test_filtering_and_best_hit_selection(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            input_fasta = temp / "families.fa"
            blast_tsv = temp / "hits.tsv"
            output_fasta = temp / "filtered.fa"
            report_tsv = temp / "removed.tsv"

            input_fasta.write_text(
                ">no_hit#Unknown description retained\nACGT\n"
                ">te_hit#LTR/Gypsy\nTGCA\n"
                ">host_hit#Unknown\nAAAA\n"
                ">best_hit#Unknown\nCCCC\n",
                encoding="utf-8",
            )
            blast_tsv.write_text(
                blast_row("te_hit#LTR/Gypsy", "1e-40", "Gypsy retrotransposon protein")
                + blast_row("host_hit#Unknown", "1e-50", "conserved metabolic enzyme")
                + blast_row("best_hit#Unknown", "1e-5", "transposase")
                + blast_row("best_hit#Unknown", "1e-60", "ribosomal protein"),
                encoding="utf-8",
            )

            kept, removed = FILTER.filter_families(
                input_fasta, blast_tsv, output_fasta, report_tsv
            )

            self.assertEqual((kept, removed), (2, 2))
            output = output_fasta.read_text(encoding="utf-8")
            report = report_tsv.read_text(encoding="utf-8")
            self.assertIn(">no_hit#Unknown description retained", output)
            self.assertIn(">te_hit#LTR/Gypsy", output)
            self.assertNotIn(">host_hit#Unknown", output)
            self.assertIn("host_hit#Unknown", report)
            self.assertIn("best_hit#Unknown\t1e-60\tribosomal protein", report)

    def test_keyword_boundaries_avoid_substring_false_positives(self):
        self.assertTrue(FILTER.is_te_protein("LINE retrotransposon polyprotein"))
        self.assertFalse(FILTER.is_te_protein("lineage-specific host protein"))
        self.assertFalse(FILTER.is_te_protein("environmental response protein"))

    def test_malformed_blast_rows_are_ignored(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            blast_tsv = Path(temp_dir) / "hits.tsv"
            blast_tsv.write_text("bad row\n" + blast_row("q1", "not-a-number", "protein"))
            self.assertEqual(FILTER.parse_blast(blast_tsv), {})


if __name__ == "__main__":
    unittest.main()
