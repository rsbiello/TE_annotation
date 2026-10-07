import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "exclude_families.py"
SPEC = importlib.util.spec_from_file_location("exclude_families", SCRIPT)
EXCLUDE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(EXCLUDE)


class ExcludeFamiliesTests(unittest.TestCase):
    def test_excludes_base_name_before_class_suffix(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            source = temp / "families.fa"
            output = temp / "filtered.fa"
            report = temp / "report.tsv"
            source.write_text(
                ">rnd-1_family-1#LINE/CR1 description\nACGT\n"
                ">rnd-2_family-2#Unknown\nTGCA\n",
                encoding="utf-8",
            )

            result = EXCLUDE.exclude_families(
                source, output, report, ["rnd-1_family-1"]
            )

            self.assertEqual(result, (1, 1))
            self.assertNotIn("rnd-1_family-1", output.read_text())
            self.assertIn(">rnd-2_family-2#Unknown", output.read_text())
            self.assertIn("rnd-1_family-1\t", report.read_text())

    def test_no_exclusions_copies_every_record(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            source = temp / "families.fa"
            output = temp / "filtered.fa"
            report = temp / "report.tsv"
            source.write_text(">a#Unknown\nAAAA\n>b#DNA\nCCCC\n")

            self.assertEqual(
                EXCLUDE.exclude_families(source, output, report, []), (2, 0)
            )
            self.assertEqual(output.read_text(), source.read_text())

    def test_missing_requested_family_is_an_error(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            source = temp / "families.fa"
            source.write_text(">a#Unknown\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "not found"):
                EXCLUDE.exclude_families(
                    source,
                    temp / "filtered.fa",
                    temp / "report.tsv",
                    ["missing"],
                )


if __name__ == "__main__":
    unittest.main()
