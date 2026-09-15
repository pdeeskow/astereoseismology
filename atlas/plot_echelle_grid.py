"""Vergleichbares Echelle-Raster aller atlasfähigen Sterne."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from atlas.plot_common import DEGREE_COLORS, load_modes, save_figure
from pbjam_bridge import adaptive_quality_threshold


def plot_echelle_grid(summaries: list[dict[str, Any]], output: str | Path) -> Path:
    usable = [(summary, load_modes(summary)) for summary in summaries]
    usable = [(summary, modes) for summary, modes in usable if modes is not None]
    if not usable:
        raise ValueError("Keine gespeicherten PBjam-Moden für das Echelle-Raster gefunden.")
    usable.sort(key=lambda item: item[0]["numax"], reverse=True)
    columns = min(4, len(usable))
    rows = math.ceil(len(usable) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(3.5 * columns, 4.0 * rows), squeeze=False)
    for axis, (summary, modes) in zip(axes.flat, usable):
        dnu = summary["deltanu_used"]
        numax = summary["numax"]
        target = summary["target"]
        quality_min = target.get("pbjam_quality_min")
        if quality_min is None:
            quality_min = adaptive_quality_threshold(
                modes,
                factor=target.get("pbjam_quality_factor", 1.10),
                floor=target.get("pbjam_quality_floor", 0.75),
                ceiling=target.get("pbjam_quality_ceiling", 2.0),
            )
        reliable = (
            (modes["quality"] >= quality_min)
            & (modes["height"] >= target.get("pbjam_height_min", 1.0))
        )
        for degree in (0, 1, 2):
            selected = reliable & (modes["l"] == degree)
            axis.scatter(
                np.mod(modes["freq"][selected], dnu),
                (modes["freq"][selected] - numax) / dnu,
                s=28,
                color=DEGREE_COLORS[degree],
                label=f"l={degree}",
                alpha=0.85,
            )
        source = "override" if summary["deltanu_source"] == "override" else "auto"
        axis.set_title(f"{summary['target']['label']}\nΔν={dnu:.2f} μHz ({source})", fontsize=10)
        axis.set_xlim(0, dnu)
        axis.set_xlabel("ν mod Δν (μHz)")
        axis.set_ylabel("(ν − νmax) / Δν")
        axis.grid(alpha=0.18)
    for axis in axes.flat[len(usable):]:
        axis.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Seismischer Atlas: Echelle-Diagramme", y=1.01)
    fig.tight_layout()
    path = save_figure(fig, output)
    plt.close(fig)
    return path