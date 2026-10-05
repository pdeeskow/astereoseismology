"""Gemeinsames Laden und Darstellen gespeicherter Atlasresultate."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


STAGE_COLORS = {
    "main_sequence": "#2878B5",
    "subgiant": "#D07C27",
    "red_giant": "#3A8D5D",
    "reference": "#555555",
}
DEGREE_COLORS = {0: "#185FA5", 1: "#0F9E66", 2: "#D85A30"}


def load_summaries(root: str | Path = "results/atlas") -> list[dict[str, Any]]:
    summaries = []
    for path in sorted(Path(root).glob("*/summary.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        data["_summary_path"] = str(path)
        summaries.append(data)
    return summaries


def stage_group(stage: str) -> str:
    normalized = stage.casefold()
    if "hauptreihe" in normalized or "main" in normalized:
        return "main_sequence"
    if "unterriese" in normalized or "subgiant" in normalized:
        return "subgiant"
    if "rgb" in normalized or "riese" in normalized or "giant" in normalized:
        return "red_giant"
    return "reference"


def stage_color(summary: dict[str, Any]) -> str:
    return STAGE_COLORS[stage_group(summary["target"]["stage"])]


def load_modes(summary: dict[str, Any]) -> dict[str, np.ndarray] | None:
    path_text = summary.get("artifacts", {}).get("pbjam_modes")
    if path_text is None:
        stem = "_".join(summary["target"]["id"].split())
        path = Path("results/tables") / f"{stem}_pbjam_modes.csv"
    else:
        path = Path(path_text)
    if not path.exists():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    modes = {
        "freq": np.asarray([float(row["freq"]) for row in rows]),
        "l": np.asarray([int(float(row["l"])) for row in rows]),
        "quality": np.asarray([float(row["quality"]) for row in rows]),
        "height": np.asarray([float(row["height"]) for row in rows]),
    }
    for key in ("height_snr_proxy", "fap_proxy"):
        if key in rows[0]:
            modes[key] = np.asarray([float(row[key]) for row in rows])
    return modes


def save_figure(fig: Any, path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(output, bbox_inches="tight")
    except PermissionError:
        output = output.with_suffix(".png")
        fig.savefig(output, bbox_inches="tight")
    return output