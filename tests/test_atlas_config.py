import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from atlas_config import load_targets
from atlas.run_all import _cached_summary, _literature_summary, _stored_summaries
from atlas_models import TargetConfig


class AtlasConfigTests(unittest.TestCase):
    def test_target_uses_shared_adaptive_quality_settings(self):
        target = TargetConfig("KIC 1", "Test", "subgiant", 5800.0)

        self.assertIsNone(target.pbjam_quality_min)
        self.assertEqual(target.pbjam_quality_factor, 1.10)
        self.assertEqual(target.pbjam_quality_floor, 0.75)
        self.assertEqual(target.pbjam_quality_ceiling, 2.0)
        self.assertEqual(target.pbjam_height_min, 1.0)
        self.assertEqual(target.pbjam_d02_quality_min, 0.65)

    def test_loads_nested_analysis_separately_from_reference(self):
        content = """
targets:
  - id: KIC 1
    label: Test
    stage: subgiant
    status: analyzed
    analysis:
      teff: 5800
      deltanu_override:
        value: 50.0
        uncertainty: 0.2
        reason: vertical ridges
    reference:
      deltanu: 49.8
      source: Example 2020
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "targets.yaml"
            path.write_text(content, encoding="utf-8")
            targets = load_targets(path)

        status, target = targets[0]
        self.assertEqual(status, "analyzed")
        self.assertEqual(target.select_deltanu(49.0, 0.5), (50.0, 0.2, "override"))
        self.assertEqual(target.reference["deltanu"], 49.8)

    def test_literature_target_becomes_reference_only_summary(self):
        content = """
targets:
  - id: Sun
    label: Sonne
    stage: Hauptreihe (Referenz)
    literature_only: true
    analysis:
      teff: 5778
    reference:
      numax: 3090.0
      deltanu: 135.1
      source: Example
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "targets.yaml"
            path.write_text(content, encoding="utf-8")
            _, target = load_targets(path)[0]

        summary = _literature_summary(target)
        self.assertEqual(summary["qc_status"], "reference")
        self.assertEqual(summary["deltanu_source"], "literature")
        self.assertEqual(summary["numax"], 3090.0)

    def test_cache_comparison_ignores_one_time_pbjam_flags(self):
        target = TargetConfig("KIC 1", "Test", "subgiant", 5800.0)
        cached_target = TargetConfig(
            "KIC 1", "Test", "subgiant", 5800.0,
            pbjam_refresh=True,
            pbjam_reuse_existing=True,
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary_dir = root / "results" / "atlas" / "KIC_1"
            summary_dir.mkdir(parents=True)
            (summary_dir / "diagnostics.npz").write_bytes(b"data")
            (summary_dir / "summary.json").write_text(
                json.dumps({"schema": 1, "target": cached_target.__dict__}),
                encoding="utf-8",
            )
            with patch("atlas.run_all._atlas_output_dir", return_value=summary_dir):
                loaded = _cached_summary(target)

        self.assertIsNotNone(loaded)

    def test_stored_summaries_include_all_targets(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for target_id in ("KIC_1", "KIC_2"):
                target_dir = root / target_id
                target_dir.mkdir()
                (target_dir / "summary.json").write_text(
                    json.dumps({"target": {"id": target_id}}), encoding="utf-8"
                )

            summaries = _stored_summaries(root)

        self.assertEqual([item["target"]["id"] for item in summaries], ["KIC_1", "KIC_2"])


if __name__ == "__main__":
    unittest.main()