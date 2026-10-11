import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[1]
    / "workflow"
    / "scripts"
    / "patch_repeatmasker_divergence.py"
)
SPEC = importlib.util.spec_from_file_location("patch_repeatmasker_divergence", SCRIPT)
PATCHER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = PATCHER
SPEC.loader.exec_module(PATCHER)


class DivergencePatcherTests(unittest.TestCase):
    def test_forces_no_cpg_recalculation_and_clamps_negative_values(self):
        source = (
            "prefix\n"
            + PATCHER.ORIGINAL_BLOCK
            + "middle\n"
            + PATCHER.OVERLAP_MARKER
            + "suffix\n"
        )
        patched = PATCHER.patch_source(source)
        self.assertIn("divCpGMod => 0", patched)
        self.assertIn("$div = 0 if ( $div < 0 );", patched)
        self.assertNotIn(PATCHER.ORIGINAL_BLOCK, patched)

    def test_rejects_an_unknown_repeatmasker_layout(self):
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            PATCHER.patch_source("unexpected source")


if __name__ == "__main__":
    unittest.main()
