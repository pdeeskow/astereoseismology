"""Datenvertrag für reproduzierbare Einzelstern- und Atlasanalysen."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal


ATLAS_SCHEMA_VERSION = 1
QCStatus = Literal["pass", "warning", "failed"]


@dataclass(frozen=True)
class DeltanuOverride:
    value: float
    uncertainty: float
    reason: str

    def __post_init__(self) -> None:
        if self.value <= 0 or self.uncertainty <= 0:
            raise ValueError("Delta-nu-Override und Unsicherheit müssen positiv sein.")
        if not self.reason.strip():
            raise ValueError("Ein Delta-nu-Override benötigt eine Begründung.")


@dataclass(frozen=True)
class TargetConfig:
    id: str
    label: str
    stage: str
    teff: float
    teff_sigma: float = 100.0
    author: str | None = None
    exptime: int | None = None
    sectors: int = 1
    fmin: float = 1.0
    fmax: float = 300.0
    oversample: int = 5
    gauss_smooth: bool = False
    echelle_replicas: int = 2
    pbjam: bool = True
    pbjam_numax_sigma: float | None = None
    pbjam_deltanu_sigma: float | None = None
    pbjam_orders: int = 7
    pbjam_quality_min: float | None = None
    pbjam_quality_factor: float = 1.10
    pbjam_quality_floor: float = 0.75
    pbjam_quality_ceiling: float = 2.0
    pbjam_height_min: float = 1.0
    pbjam_d02_quality_min: float = 0.65
    pbjam_refresh: bool = False
    pbjam_reuse_existing: bool = False
    bp_rp: float | None = None
    bp_rp_sigma: float = 0.05
    literature_only: bool = False
    deltanu_override: DeltanuOverride | None = None
    reference: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.id.strip() or not self.label.strip() or not self.stage.strip():
            raise ValueError("Target-ID, Label und Entwicklungsstadium sind erforderlich.")
        if self.teff <= 0 or self.teff_sigma <= 0:
            raise ValueError("Teff und Teff-Unsicherheit müssen positiv sein.")
        if self.sectors < 1 or self.oversample < 1 or self.pbjam_orders < 1:
            raise ValueError("Sektoren, Oversampling und PBjam-Ordnungen müssen mindestens 1 sein.")
        if self.echelle_replicas not in (2, 3):
            raise ValueError("Es werden zwei oder drei Echelle-Replikate unterstützt.")
        if self.fmin <= 0 or self.fmax <= self.fmin:
            raise ValueError("Der Frequenzbereich ist ungültig.")
        optional_uncertainties = (self.pbjam_numax_sigma, self.pbjam_deltanu_sigma)
        if any(value is not None and value <= 0 for value in optional_uncertainties):
            raise ValueError("Optionale PBjam-Unsicherheiten müssen positiv sein.")
        if self.pbjam_quality_min is not None and self.pbjam_quality_min < 0:
            raise ValueError("Der explizite PBjam-Quality-Cut darf nicht negativ sein.")
        if self.pbjam_quality_factor <= 0 or self.pbjam_quality_floor < 0:
            raise ValueError("PBjam-Quality-Faktor und Untergrenze müssen gültig sein.")
        if self.pbjam_quality_ceiling < self.pbjam_quality_floor:
            raise ValueError("Die PBjam-Quality-Obergrenze muss mindestens der Untergrenze entsprechen.")
        if self.pbjam_height_min < 0 or self.pbjam_d02_quality_min < 0:
            raise ValueError("PBjam-Höhen- und Paargrenzen dürfen nicht negativ sein.")

    def select_deltanu(self, automatic: float, automatic_sigma: float) -> tuple[float, float, str]:
        if automatic <= 0 or automatic_sigma <= 0:
            raise ValueError("Automatisches Delta-nu und Unsicherheit müssen positiv sein.")
        if self.deltanu_override is None:
            return automatic, automatic_sigma, "auto"
        return self.deltanu_override.value, self.deltanu_override.uncertainty, "override"


@dataclass(frozen=True)
class QCFlag:
    code: str
    status: QCStatus
    message: str
    value: float | int | str | None = None

    def __post_init__(self) -> None:
        if self.status not in ("pass", "warning", "failed"):
            raise ValueError(f"Unbekannter QC-Status: {self.status}")


@dataclass
class AtlasResult:
    target: TargetConfig
    numax: float
    numax_sigma: float
    deltanu_auto: float
    deltanu_auto_sigma: float
    deltanu_used: float
    deltanu_used_sigma: float
    deltanu_source: Literal["auto", "override"]
    stellar_parameters: dict[str, float]
    d02: float | None = None
    d02_sigma: float | None = None
    d02_pairs: int = 0
    crossing_candidates: list[dict[str, Any]] = field(default_factory=list)
    mode_counts: dict[str, int | float | None] = field(default_factory=dict)
    qc: list[QCFlag] = field(default_factory=list)
    artifacts: dict[str, str] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    schema: int = ATLAS_SCHEMA_VERSION

    @property
    def qc_status(self) -> QCStatus:
        statuses = {flag.status for flag in self.qc}
        if "failed" in statuses:
            return "failed"
        if "warning" in statuses:
            return "warning"
        return "pass"

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["qc_status"] = self.qc_status
        return data

    def write_json(self, path: str | Path) -> None:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")


def target_from_dict(data: dict[str, Any]) -> TargetConfig:
    """Erzeugt eine Target-Konfiguration aus einem JSON-/YAML-kompatiblen Mapping."""
    values = dict(data)
    override = values.get("deltanu_override")
    if override is not None and not isinstance(override, DeltanuOverride):
        values["deltanu_override"] = DeltanuOverride(**override)
    return TargetConfig(**values)
