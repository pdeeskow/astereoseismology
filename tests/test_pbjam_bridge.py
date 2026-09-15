import unittest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np

from pbjam_bridge import (
    adaptive_quality_threshold,
    _apply_jax_compatibility,
    _apply_numpy_compatibility,
    _pbjam_version,
    extract_modes,
    filter_reliable_modes,
    find_avoided_crossings,
    run_pbjam_modeid,
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

    def test_adapts_quality_threshold_and_filters_low_height(self):
        modes = extract_modes(self.modeid, self.peakbag, deltanu=10.0)
        modes["quality"] = np.array([0.7, 0.74, 0.9])
        modes["height"] = np.array([4.0, 3.0, 0.5])

        self.assertAlmostEqual(adaptive_quality_threshold(modes), 0.814)
        reliable = filter_reliable_modes(modes)

        self.assertEqual(len(reliable["freq"]), 0)

    def test_finds_disturbed_l1_spacing(self):
        result = find_avoided_crossings(np.array([90.0, 100.0, 116.0, 126.0]), 10.0)

        np.testing.assert_allclose(result["crossings"], [108.0])

    def test_explicitly_reuses_complete_cache_with_mismatched_signature(self):
        modes = extract_modes(self.modeid, self.peakbag, deltanu=10.0)
        with tempfile.TemporaryDirectory() as directory:
            cache_dir = Path(directory) / "cache"
            cache_dir.mkdir()
            (cache_dir / "KIC_1_modes.json").write_text(
                json.dumps({"schema": 1, "signature": "old", "modes": {key: value.tolist() for key, value in modes.items()}}),
                encoding="utf-8",
            )
            with patch("pbjam_bridge._pbjam_version") as version_check:
                loaded = run_pbjam_modeid(
                    np.array([1.0, 2.0, 3.0]),
                    np.array([1.0, 1.0, 1.0]),
                    numax=(100.0, 1.0),
                    deltanu=(10.0, 0.1),
                    teff=(5800.0, 100.0),
                    star_id="KIC 1",
                    outdir=directory,
                    reuse_existing_cache=True,
                )

        version_check.assert_not_called()
        np.testing.assert_allclose(loaded["freq"], modes["freq"])

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