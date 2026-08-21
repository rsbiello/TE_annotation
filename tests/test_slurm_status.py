import importlib.util
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "slurm_status.py"
SPEC = importlib.util.spec_from_file_location("slurm_status", SCRIPT)
STATUS = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(STATUS)


class SlurmStatusTests(unittest.TestCase):
    def test_terminal_states(self):
        self.assertEqual(STATUS.classify_state("COMPLETED", "0:0"), "success")
        self.assertEqual(STATUS.classify_state("FAILED", "1:0"), "failed")
        self.assertEqual(STATUS.classify_state("TIMEOUT", "0:0"), "failed")
        self.assertEqual(STATUS.classify_state("OUT_OF_MEMORY", "0:9"), "failed")
        self.assertEqual(STATUS.classify_state("CANCELLED by 123", "0:15"), "failed")

    def test_active_and_unknown_states_are_polled_again(self):
        self.assertEqual(STATUS.classify_state("RUNNING", "0:0"), "running")
        self.assertEqual(STATUS.classify_state("PENDING", "0:0"), "running")
        self.assertEqual(STATUS.classify_state(""), "running")

    @patch.object(STATUS, "current_queue_states", return_value=[])
    @patch.object(STATUS, "accounting_record", return_value=("TIMEOUT", "0:0"))
    def test_completed_job_is_read_from_accounting(self, accounting, queue):
        self.assertEqual(STATUS.job_status("641448"), "failed")
        queue.assert_called_once_with("641448")
        accounting.assert_called_once_with("641448")

    @patch.object(STATUS, "current_queue_states", return_value=[])
    @patch.object(STATUS, "accounting_record", return_value=None)
    def test_accounting_delay_does_not_cause_false_failure(self, accounting, queue):
        self.assertEqual(STATUS.job_status("123;cluster"), "running")
        queue.assert_called_once_with("123")
        accounting.assert_called_once_with("123")


if __name__ == "__main__":
    unittest.main()
