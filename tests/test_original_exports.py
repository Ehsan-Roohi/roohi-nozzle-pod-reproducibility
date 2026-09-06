"""Small parser regression tests plus a full fifteen-case lineage check."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from verify_original_exports import VARIABLES, parse_zones, verify


class OriginalExportTests(unittest.TestCase):
    def fixture(self, row=None):
        return ('TITLE="test"\nVariables=' + ','.join(VARIABLES)
                + '\nZONE T="MAIN", I=1,J=1\nDATAPACKING=POINT\n'
                + (row or ' '.join(str(i) for i in range(12))) + '\n')

    def test_point_format(self):
        self.assertEqual(parse_zones(self.fixture())[0].shape, (1, 12))

    def test_quoted_variables(self):
        text = self.fixture().replace(','.join(VARIABLES), '\n'.join('"'+v+'"' for v in VARIABLES))
        self.assertEqual(parse_zones(text)[0][0, 7], 7)

    def test_incomplete_zone(self):
        with self.assertRaises(ValueError):
            parse_zones(self.fixture().replace('I=1', 'I=2'))

    def test_wrong_variable_order(self):
        with self.assertRaises(ValueError):
            parse_zones(self.fixture().replace('U,V', 'V,U'))

    def test_nonfinite(self):
        with self.assertRaises(ValueError):
            parse_zones(self.fixture('0 1 2 3 4 5 6 7 8 9 10 nan'))

    def test_wrong_columns(self):
        with self.assertRaises(ValueError):
            parse_zones(self.fixture('0 1 2'))

    def test_block_rejected(self):
        with self.assertRaises(ValueError):
            parse_zones(self.fixture().replace('POINT', 'BLOCK'))

    def test_no_zone(self):
        with self.assertRaises(ValueError):
            parse_zones('Variables=' + ','.join(VARIABLES))

    def test_all_fifteen_original_exports(self):
        report = verify()
        self.assertEqual(report['verified_cases'], 15)
        self.assertEqual(report['verified_values'], 792360)
        self.assertTrue(report['all_match_at_published_precision'])


if __name__ == '__main__':
    unittest.main()
