"""Skalenrelation und Vergleich mit publizierten Werten."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from atlas.plot_common import save_figure, stage_color


def plot_scaling(summaries: list[dict[str, Any]], output: str | Path) -> Path:
    if not summaries:
        raise ValueError("Keine Atlasresultate für den Skalenrelationstest gefunden.")
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    relation_axis = axes[0, 0]
    numax_values = np.logspace(1, 3.6, 300)
    relation_axis.plot(numax_values, 0.263 * numax_values**0.772, color="#333333", label="Stello et al. 2009")
    for summary in summaries:
        if summary.get("numax") is None or summary.get("deltanu_used") is None:
            continue
        marker = "*" if summary["target"].get("literature_only") else "o"
        relation_axis.errorbar(summary["numax"], summary["deltanu_used"], xerr=summary.get("numax_sigma"), yerr=summary.get("deltanu_used_sigma"), fmt=marker, color=stage_color(summary))
    relation_axis.set(xscale="log", yscale="log", xlabel="νmax (μHz)", ylabel="Δν (μHz)", title="Skalenrelation")
    relation_axis.legend(frameon=False)
    relation_axis.grid(which="both", alpha=0.2)

    metrics = [
        ("numax", "numax", "νmax (μHz)"),
        ("deltanu_used", "deltanu", "Δν (μHz)"),
        ("mass", "mass", "M / M☉"),
        ("radius", "radius", "R / R☉"),
    ]
    for axis, (measured_key, reference_key, label) in zip(axes.flat[1:], metrics):
        pairs = []
        for summary in summaries:
            if summary["target"].get("literature_only"):
                continue
            reference = summary["target"].get("reference", {}).get(reference_key)
            measured = summary.get(measured_key)
            if measured is None:
                measured = summary.get("stellar_parameters", {}).get(measured_key)
            if reference is not None and measured is not None:
                pairs.append((float(reference), float(measured), summary))
        if not pairs:
            axis.text(0.5, 0.5, "Keine Referenzwerte", ha="center", va="center", transform=axis.transAxes)
            axis.set_axis_off()
            continue
        limits = [value for reference, measured, _ in pairs for value in (reference, measured)]
        low, high = min(limits), max(limits)
        padding = 0.08 * (high - low or high or 1.0)
        axis.plot([low - padding, high + padding], [low - padding, high + padding], color="#777777", linestyle="--")
        for reference, measured, summary in pairs:
            axis.scatter(reference, measured, color=stage_color(summary))
        axis.set(xlabel=f"publiziert: {label}", ylabel=f"gemessen: {label}", title=label)
        axis.grid(alpha=0.2)
    axes.flat[-1].set_visible(False)
    fig.suptitle("Skalenrelationen und Literaturvergleich")
    fig.tight_layout()
    path = save_figure(fig, output)
    plt.close(fig)
    return path