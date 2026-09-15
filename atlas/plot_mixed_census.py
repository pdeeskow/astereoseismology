"""QC-gekennzeichneter Zensus von Avoided-Crossing-Kandidaten."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from atlas.plot_common import save_figure, stage_color


def plot_mixed_census(summaries: list[dict[str, Any]], output: str | Path) -> Path:
    usable = [summary for summary in summaries if summary.get("stellar_parameters", {}).get("logg") is not None]
    if not usable:
        raise ValueError("Keine log-g-Werte für den Mixed-Mode-Zensus gefunden.")
    fig, axis = plt.subplots(figsize=(7.2, 5.2))
    for summary in usable:
        count = len(summary.get("crossing_candidates", []))
        logg = summary["stellar_parameters"]["logg"]
        mode_count = summary.get("mode_counts", {}).get("l1", 0)
        axis.scatter(logg, count, s=40 + 8 * mode_count, color=stage_color(summary), alpha=0.85)
        axis.annotate(summary["target"]["label"], (logg, count), xytext=(5, 5), textcoords="offset points")
    axis.set(xlabel="log g", ylabel="Anzahl Crossing-Kandidaten", title="Heuristischer Mixed-Mode-Kandidatenzensus")
    axis.invert_xaxis()
    axis.grid(alpha=0.2)
    axis.text(0.01, 0.01, "Punktgröße: Anzahl zuverlässiger l=1-Moden", transform=axis.transAxes, fontsize=9)
    fig.tight_layout()
    path = save_figure(fig, output)
    plt.close(fig)
    return path