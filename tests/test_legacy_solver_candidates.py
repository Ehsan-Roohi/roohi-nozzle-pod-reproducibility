import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "legacy_solver_candidates"


class LegacySolverCandidateTests(unittest.TestCase):
    def test_retained_file_hashes(self):
        manifest = json.loads((CANDIDATES / "manifest.json").read_text())
        for relative, expected in manifest["retained_files_sha256"].items():
            actual = hashlib.sha256((CANDIDATES / relative).read_bytes()).hexdigest()
            self.assertEqual(actual, expected, relative)

    def test_p25_candidate_input(self):
        text = (CANDIDATES / "nozzle_2009" / "Inputdata.txt").read_text()
        self.assertIn("500.\t!FTMP", text)
        self.assertIn("1.E5\t!PIN", text)
        self.assertIn("25.E3\t!POUT", text)
        self.assertIn("100\t!NCX", text)
        self.assertIn("30\t!NCY", text)

    def test_candidate_is_not_claimed_as_exact_producer(self):
        report = (ROOT / "docs" / "LEGACY_SOLVER_RECOVERY.md").read_text()
        self.assertIn("candidate rather than", report)
        self.assertIn("QX", report)
        self.assertIn("QY", report)
        self.assertIn("Txy", report)


if __name__ == "__main__":
    unittest.main()
