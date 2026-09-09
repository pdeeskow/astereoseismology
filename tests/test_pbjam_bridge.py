import unittest

import numpy as np

from pbjam_bridge import (
    _apply_jax_compatibility,
    _apply_numpy_compatibility,
    _pbjam_version,
    extract_modes,
    filter_reliable_modes,
    find_avoided_crossings,
)


class PBjamBridgeTests(unittest.TestCase):
    def setUp(self):
        self.modeid = {
            "ell": np.array([0, 1, 1]),
            "enn": np.array([10, -1, -1]),
            "summary": {"freq": np.array([[100.0, 105.0, 112.0], [0.2, 0.3, 0.3]])},
        }
        self.peakbag = {
            "ell": np.array([[0, 1, 1]]),
            "summary": {
                "freq": np.array([[100.1, 105.2, 112.1], [0.1, 1.0, 0.05]]),
                "height": np.array([[4.0, 3.0, 2.0], [0.4, 0.3, 0.2]]),
                "width": np.array([[1.0, 0.5, 0.4], [0.1, 0.05, 0.04]]),
            },
        }

    def test_extracts_pbjam_2_result(self):
        modes = extract_modes(self.modeid, self.peakbag, deltanu=10.0)

        np.testing.assert_array_equal(modes["l"], [0, 1, 1])
        np.testing.assert_array_equal(modes["n"], [10, -1, -1])
        np.testing.assert_allclose(modes["quality"], [3.0, 0.3, 6.0])

    def test_filters_on_prior_posterior_ratio(self):
        modes = extract_modes(self.modeid, self.peakbag, deltanu=10.0)

        reliable = filter_reliable_modes(modes, quality_min=2.0)

        np.testing.assert_allclose(reliable["freq"], [100.1, 112.1])

    def test_finds_disturbed_l1_spacing(self):
        result = find_avoided_crossings(np.array([90.0, 100.0, 116.0, 126.0]), 10.0)

        np.testing.assert_allclose(result["crossings"], [108.0])

    def test_restores_numpy_trapz_for_pbjam_2(self):
        original = getattr(np, "trapz", None)
        if original is not None:
            delattr(np, "trapz")
        try:
            self.assertTrue(_apply_numpy_compatibility())
            np.testing.assert_allclose(np.trapz([0.0, 1.0], [0.0, 1.0]), 0.5)
            self.assertFalse(_apply_numpy_compatibility())
        finally:
            if original is None:
                delattr(np, "trapz")
            else:
                setattr(np, "trapz", original)

    def test_restores_old_jax_clip_keywords_for_pbjam_2(self):
        import jax.numpy as jnp

        original = jnp.clip
        try:
            self.assertTrue(_apply_jax_compatibility())
            values = jnp.array([-1.0, 0.5, 2.0])
            np.testing.assert_allclose(jnp.clip(values, a_min=0.0, a_max=1.0), [0.0, 0.5, 1.0])
            self.assertFalse(_apply_jax_compatibility())
        finally:
            setattr(jnp, "clip", original)

    def test_installed_pbjam_runtime_is_compatible(self):
        self.assertRegex(_pbjam_version(), r"^2\.")


if __name__ == "__main__":
    unittest.main()