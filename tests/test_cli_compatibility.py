import sys
import unittest
from unittest.mock import patch

import asteroseismologie as pipeline


class CLICompatibilityTests(unittest.TestCase):
    def test_existing_options_map_to_target_config(self):
        argv = [
            "asteroseismologie.py", "--tic", "KIC 10963065", "--name", "Test",
            "--teff", "6140", "--teff-sigma", "77", "--exptime", "60",
            "--sectors", "8", "--oversample", "2", "--pbjam", "--bp-rp", "0.70026",
            "--pbjam-quality-min", "0", "--echelle-replicas", "3",
            "--reuse-existing-pbjam",
        ]
        with patch.object(sys, "argv", argv):
            config = pipeline._target_config_from_args(pipeline._parse_args())

        self.assertEqual(config.id, "KIC 10963065")
        self.assertEqual(config.teff_sigma, 77.0)
        self.assertTrue(config.pbjam)
        self.assertEqual(config.pbjam_quality_min, 0.0)
        self.assertEqual(config.echelle_replicas, 3)
        self.assertTrue(config.pbjam_reuse_existing)

    def test_main_keeps_single_target_adapter(self):
        with patch.object(pipeline, "_parse_args") as parse_args, patch.object(pipeline, "analyze_target") as analyze:
            parse_args.return_value = pipeline.argparse.Namespace(
                tic="KIC 1", name="Test", teff=5800.0, teff_sigma=100.0,
                author=None, exptime=60, sectors=1, fmin=1.0, fmax=300.0,
                oversample=2, gauss_smooth=False, echelle_replicas=2, pbjam=False,
                pbjam_numax_sigma=None, pbjam_deltanu_sigma=None, pbjam_orders=7,
                pbjam_quality_min=2.0, pbjam_refresh=False, reuse_existing_pbjam=False,
                bp_rp=None, bp_rp_sigma=0.05,
            )

            pipeline.main()

        config = analyze.call_args.args[0]
        self.assertEqual(config.id, "KIC 1")
        self.assertFalse(analyze.call_args.kwargs.get("write_atlas", False))


if __name__ == "__main__":
    unittest.main()