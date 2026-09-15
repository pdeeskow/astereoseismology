import json
import tempfile
import unittest
from pathlib import Path

from atlas_models import AtlasResult, DeltanuOverride, QCFlag, TargetConfig


class AtlasModelTests(unittest.TestCase):
    def test_selects_automatic_deltanu_without_override(self):
        target = TargetConfig("KIC 1", "Test", "subgiant", 5800.0)

        self.assertEqual(target.select_deltanu(50.0, 0.5), (50.0, 0.5, "auto"))

    def test_selects_documented_deltanu_override(self):
        target = TargetConfig(
            "KIC 1",
            "Test",
            "subgiant",
            5800.0,
            deltanu_override=DeltanuOverride(49.5, 0.2, "Ridges visually vertical"),
        )

        self.assertEqual(target.select_deltanu(50.0, 0.5), (49.5, 0.2, "override"))

    def test_rejects_undocumented_override(self):
        with self.assertRaises(ValueError):
            DeltanuOverride(49.5, 0.2, " ")

    def test_writes_versioned_result_with_aggregate_qc_status(self):
        target = TargetConfig("KIC 1", "Test", "subgiant", 5800.0)
        result = AtlasResult(
            target=target,
            numax=900.0,
            numax_sigma=20.0,
            deltanu_auto=50.0,
            deltanu_auto_sigma=0.5,
            deltanu_used=50.0,
            deltanu_used_sigma=0.5,
            deltanu_source="auto",
            stellar_parameters={"mass": 1.2},
            qc=[QCFlag("few_modes", "warning", "Only four reliable modes", 4)],
        )

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "summary.json"
            result.write_json(output)
            saved = json.loads(output.read_text(encoding="utf-8"))

        self.assertEqual(saved["schema"], 1)
        self.assertEqual(saved["qc_status"], "warning")
        self.assertEqual(saved["target"]["id"], "KIC 1")


if __name__ == "__main__":
    unittest.main()
