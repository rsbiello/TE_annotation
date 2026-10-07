import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def dependencies(environment_name):
    environment = yaml.safe_load(
        (ROOT / "workflow" / "envs" / environment_name).read_text()
    )
    return set(environment["dependencies"])


class CondaConfigurationTests(unittest.TestCase):
    def test_repeatmodeler_versions_are_pinned(self):
        deps = dependencies("repeatmodeler.yaml")
        self.assertIn("repeatmodeler =2.0.9", deps)
        self.assertIn("repeatmasker =4.2.4", deps)
        self.assertIn("famdb =3.0.0", deps)
        self.assertIn("ltr_retriever =2.9.0", deps)

    def test_repeatmasker_version_is_pinned(self):
        deps = dependencies("repeatmasker.yaml")
        self.assertIn("repeatmasker =4.2.4", deps)
        self.assertIn("famdb =3.0.0", deps)

    def test_tetrimmer_version_is_pinned(self):
        deps = dependencies("tetrimmer.yaml")
        self.assertIn("tetrimmer =1.7.2", deps)
        self.assertIn("repeatmodeler =2.0.9", deps)
        self.assertIn("repeatmasker =4.2.4", deps)

    def test_workflow_has_no_container_runtime(self):
        snakefile = (ROOT / "workflow" / "Snakefile").read_text().lower()
        for term in ("docker", "podman", "tetools", "container_runtime"):
            self.assertNotIn(term, snakefile)

    def test_famdb_repeat_peps_subcommand(self):
        snakefile = (ROOT / "workflow" / "Snakefile").read_text()
        self.assertIn(" repeat_peps |", snakefile)
        self.assertNotIn(" repeatpeps |", snakefile)

    def test_snakemake_script_has_no_future_import(self):
        script = (ROOT / "workflow" / "scripts" / "split_fasta.py").read_text()
        self.assertNotIn("from __future__ import", script)

    def test_repeatmasker_uses_famdb_executable_wrapper(self):
        snakefile = (ROOT / "workflow" / "Snakefile").read_text()
        repeatmasker_rule = snakefile.split("rule repeatmasker:", 1)[1].split(
            "# 04A TANDEM REPEATS FINDER", 1
        )[0]
        self.assertIn("-lib {input.library:q}", repeatmasker_rule)
        self.assertIn('-famdb_dir "$FAMDB_WRAPPER_DIR"', repeatmasker_rule)
        self.assertIn('exec %q -i %q "$@"', repeatmasker_rule)
        self.assertNotIn("-famdb_dir {params.famdb_data}", repeatmasker_rule)

    def test_landscape_adds_repeatmasker_modules_to_perl_path(self):
        snakefile = (ROOT / "workflow" / "Snakefile").read_text()
        landscape_rule = snakefile.split("rule landscape:", 1)[1]
        self.assertIn("RepeatMaskerConfig.pm", landscape_rule)
        self.assertIn('export PERL5LIB="$RM_LIB_DIR', landscape_rule)


if __name__ == "__main__":
    unittest.main()
