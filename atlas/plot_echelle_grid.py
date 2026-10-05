"""Vergleichbares Echelle-Raster aller atlasfähigen Sterne."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from atlas.plot_common import DEGREE_COLORS, load_modes, save_figure
from pbjam_bridge import adaptive_quality_threshold, tier_modes


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
        classified = tier_modes(
            modes,
            fap_gold=target.get("pbjam_fap_gold", 0.01),
            fap_silver=target.get("pbjam_fap_silver", 0.1),
            quality_min=quality_min,
            dnu=dnu,
            ridge_tol_uHz=target.get("pbjam_ridge_tol_uHz", 1.5),
            d02_fraction=target.get("pbjam_d02_fraction"),
            sequence_tolerance=target.get("pbjam_sequence_tolerance", 0.10),
            sequence_minimum=target.get("pbjam_sequence_minimum", 3),
        )
        for degree in (0, 1, 2):
            for confidence in ("silver", "gold"):
                selected = (classified["tier"] == confidence) & (classified["l"] == degree)
                axis.scatter(
                    np.mod(classified["freq"][selected], dnu),
                    (classified["freq"][selected] - numax) / dnu,
                    s=32,
                    facecolors=DEGREE_COLORS[degree] if confidence == "gold" else "none",
                    edgecolors=DEGREE_COLORS[degree],
                    linewidths=1.2,
                    alpha=0.9,
                )
        source = "override" if summary["deltanu_source"] == "override" else "auto"
        axis.set_title(f"{summary['target']['label']}\nΔν={dnu:.2f} μHz ({source})", fontsize=10)
        axis.set_xlim(0, dnu)
        axis.set_xlabel("ν mod Δν (μHz)")
        axis.set_ylabel("(ν − νmax) / Δν")
        axis.grid(alpha=0.18)
    for axis in axes.flat[len(usable):]:
        axis.set_visible(False)
    handles = [
        Line2D([], [], marker="o", linestyle="none", color=DEGREE_COLORS[degree], label=f"l={degree}")
        for degree in (0, 1, 2)
    ]
    handles.extend([
        Line2D([], [], marker="o", linestyle="none", color="#333333", label="Gold"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none", color="#333333", label="Silber"),
    ])
    fig.suptitle("Seismischer Atlas: Echelle-Diagramme", y=0.98)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=5, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    path = save_figure(fig, output)
    plt.close(fig)
    return path