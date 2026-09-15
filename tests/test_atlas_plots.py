import csv
import tempfile
import unittest
from pathlib import Path

from atlas.plot_cd_diagram import plot_cd_diagram
from atlas.plot_echelle_grid import plot_echelle_grid
from atlas.plot_mixed_census import plot_mixed_census
from atlas.plot_scaling import plot_scaling


class AtlasPlotTests(unittest.TestCase):
    def test_all_plots_write_nonempty_pdfs_from_cached_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            modes_path = root / "modes.csv"
            with modes_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(("freq", "freq_err", "l", "n", "height", "height_err", "width", "width_err", "quality"))
                writer.writerows([
                    (96.0, 0.2, 2, 9, 2.0, 0.2, 0.5, 0.1, 3.0),
                    (100.0, 0.1, 0, 10, 4.0, 0.2, 0.5, 0.1, 4.0),
                    (105.0, 0.2, 1, -1, 3.0, 0.2, 0.5, 0.1, 3.0),
                    (110.0, 0.1, 0, 11, 4.0, 0.2, 0.5, 0.1, 4.0),
                    (116.0, 0.2, 1, -1, 3.0, 0.2, 0.5, 0.1, 3.0),
                ])
            summary = {
                "target": {
                    "id": "KIC 1", "label": "Test", "stage": "früher Unterriese",
                    "pbjam_quality_min": 2.0,
                    "reference": {"numax": 105.0, "deltanu": 10.0, "mass": 1.2, "radius": 2.0},
                },
                "numax": 105.0, "numax_sigma": 1.0,
                "deltanu_used": 10.0, "deltanu_used_sigma": 0.1, "deltanu_source": "auto",
                "d02": 4.0, "d02_sigma": 0.2,
                "stellar_parameters": {"mass": 1.2, "radius": 2.0, "logg": 3.9},
                "crossing_candidates": [{"midpoint": 110.5}],
                "mode_counts": {"l1": 2},
                "artifacts": {"pbjam_modes": str(modes_path)},
            }
            outputs = [
                plot_echelle_grid([summary], root / "echelle.pdf"),
                plot_cd_diagram([summary], root / "cd.pdf"),
                plot_scaling([summary], root / "scaling.pdf"),
                plot_mixed_census([summary], root / "mixed.pdf"),
            ]

            self.assertTrue(all(path.stat().st_size > 1000 for path in outputs))


if __name__ == "__main__":
    unittest.main()