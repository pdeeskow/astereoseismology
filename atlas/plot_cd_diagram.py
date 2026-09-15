"""Christensen-Dalsgaard-Diagramm aus gespeicherten Atlaswerten."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from atlas.plot_common import save_figure, stage_color


def plot_cd_diagram(summaries: list[dict[str, Any]], output: str | Path) -> Path:
    usable = [summary for summary in summaries if summary.get("d02") and summary.get("deltanu_used")]
    if not usable:
        raise ValueError("Keine belastbaren delta-nu-02-Messungen für das C-D-Diagramm gefunden.")
    fig, axis = plt.subplots(figsize=(7.2, 5.4))
    for summary in usable:
        axis.errorbar(
            summary["deltanu_used"], summary["d02"],
            xerr=summary.get("deltanu_used_sigma"), yerr=summary.get("d02_sigma"),
            fmt="o", color=stage_color(summary), capsize=3,
        )
        axis.annotate(summary["target"]["label"], (summary["deltanu_used"], summary["d02"]), xytext=(5, 5), textcoords="offset points")
    axis.set(xscale="log", yscale="log", xlabel="Δν (μHz)", ylabel="δν₀₂ (μHz)", title="C-D-Diagramm")
    axis.grid(which="both", alpha=0.2)
    fig.tight_layout()
    path = save_figure(fig, output)
    plt.close(fig)
    return path