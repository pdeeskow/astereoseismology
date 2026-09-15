"""Laden der deklarativen Targetliste für den seismischen Atlas."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from atlas_models import TargetConfig, target_from_dict


def load_targets(path: str | Path) -> list[tuple[str, TargetConfig]]:
    source = Path(path)
    document = yaml.safe_load(source.read_text(encoding="utf-8"))
    entries = document.get("targets") if isinstance(document, dict) else None
    if not isinstance(entries, list):
        raise ValueError("Die Atlas-Konfiguration benötigt eine Liste 'targets'.")

    targets: list[tuple[str, TargetConfig]] = []
    seen: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Jeder Target-Eintrag muss ein Mapping sein.")
        status = str(entry.get("status", "pending"))
        analysis = entry.get("analysis", {})
        if not isinstance(analysis, dict):
            raise ValueError("Der Block 'analysis' muss ein Mapping sein.")
        values: dict[str, Any] = {
            "id": entry.get("id"),
            "label": entry.get("label"),
            "stage": entry.get("stage"),
            "reference": entry.get("reference", {}),
            "literature_only": entry.get("literature_only", False),
            **analysis,
        }
        target = target_from_dict(values)
        normalized_id = " ".join(target.id.upper().split())
        if normalized_id in seen:
            raise ValueError(f"Doppelte Target-ID: {target.id}")
        seen.add(normalized_id)
        targets.append((status, target))
    return targets