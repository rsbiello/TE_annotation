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

    def test_workflow_has_no_container_runtime(self):
        snakefile = (ROOT / "workflow" / "Snakefile").read_text().lower()
        for term in ("docker", "podman", "tetools", "container_runtime"):
            self.assertNotIn(term, snakefile)


if __name__ == "__main__":
    unittest.main()
