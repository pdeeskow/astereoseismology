import unittest

import numpy as np

from atlas_metrics import crossing_candidates, measure_d02


class AtlasMetricTests(unittest.TestCase):
    def test_measures_d02_from_matching_radial_orders(self):
        modes = {
            "freq": np.array([100.0, 96.0, 110.0, 105.5, 130.0]),
            "freq_err": np.array([0.2, 0.3, 0.1, 0.2, 0.1]),
            "l": np.array([0, 2, 0, 2, 1]),
            "n": np.array([10, 9, 11, 10, -1]),
        }

        result = measure_d02(modes)

        self.assertEqual(result.pairs, 2)
        np.testing.assert_allclose(result.values, [4.0, 4.5])
        self.assertGreater(result.uncertainty, 0.0)
        self.assertGreater(result.value, 4.0)
        self.assertLess(result.value, 4.5)

    def test_returns_no_d02_for_ambiguous_order_pair(self):
        modes = {
            "freq": np.array([100.0, 96.0, 96.2]),
            "freq_err": np.array([0.2, 0.3, 0.3]),
            "l": np.array([0, 2, 2]),
            "n": np.array([10, 9, 9]),
        }

        result = measure_d02(modes)

        self.assertIsNone(result.value)
        self.assertEqual(result.pairs, 0)

    def test_d02_uses_its_own_loose_quality_selection(self):
        modes = {
            "freq": np.array([100.0, 98.0, 110.0, 107.5, 120.0, 117.0]),
            "freq_err": np.full(6, 0.1),
            "l": np.array([0, 2, 0, 2, 0, 2]),
            "n": np.array([10, 9, 11, 10, 12, 11]),
            "quality": np.array([0.65, 0.7, 2.0, 0.64, 0.9, 0.8]),
        }

        result = measure_d02(modes)

        self.assertEqual(result.pairs, 2)
        np.testing.assert_allclose(result.values, [2.0, 3.0])
        self.assertEqual(result.value, 2.5)
        self.assertEqual(result.uncertainty, 0.5)

    def test_reports_crossing_candidate_evidence(self):
        candidates = crossing_candidates(np.array([90.0, 100.0, 116.0, 126.0]), 10.0)

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["midpoint"], 108.0)
        self.assertAlmostEqual(candidates[0]["relative_deviation"], 0.6)


if __name__ == "__main__":
    unittest.main()