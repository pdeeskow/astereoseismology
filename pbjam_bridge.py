"""Optionale Schnittstelle zwischen der Pipeline und PBjam 2.x."""

from __future__ import annotations

import csv
import hashlib
import inspect
import json
import re
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np


PBJAM_MIN_VERSION = (2, 0)
BASE_MODE_KEYS = ("freq", "freq_err", "l", "n", "height", "height_err", "width", "width_err", "quality")
MODE_KEYS = BASE_MODE_KEYS + ("height_snr_proxy", "fap_proxy")
CACHE_SCHEMA_VERSION = 4


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "star"


def _pbjam_version() -> str:
    try:
        installed = version("pbjam")
    except PackageNotFoundError as exc:
        raise RuntimeError(
            "PBjam ist nicht installiert. Installation: uv sync --extra pbjam"
        ) from exc

    numbers = tuple(int(part) for part in re.findall(r"\d+", installed)[:2])
    if numbers < PBJAM_MIN_VERSION:
        raise RuntimeError(f"PBjam >= 2.0 ist erforderlich, installiert ist {installed}.")

    dynesty_version = version("dynesty")
    dynesty_major = int(re.match(r"\d+", dynesty_version).group())
    if dynesty_major >= 3:
        raise RuntimeError(
            f"PBjam {installed} ist nicht mit Dynesty {dynesty_version} kompatibel. "
            "Bitte `uv sync --extra pbjam` ausführen (benötigt Dynesty < 3)."
        )
    return installed


def _apply_numpy_compatibility() -> bool:
    """Stellt den von PBjam 2.0 verwendeten NumPy-Alias bereit."""
    if not hasattr(np, "trapz") and hasattr(np, "trapezoid"):
        setattr(np, "trapz", np.trapezoid)
        return True
    return False


def _apply_jax_compatibility() -> bool:
    """Übersetzt die von PBjam 2.0 verwendeten jnp.clip-Schlüssel."""
    import jax.numpy as jnp

    if "a_min" in inspect.signature(jnp.clip).parameters:
        return False

    original_clip = jnp.clip

    def compatible_clip(
        array: Any,
        minimum: Any = None,
        maximum: Any = None,
        *,
        a_min: Any = None,
        a_max: Any = None,
    ) -> Any:
        if minimum is not None and a_min is not None:
            raise TypeError("clip() erhielt sowohl minimum als auch a_min")
        if maximum is not None and a_max is not None:
            raise TypeError("clip() erhielt sowohl maximum als auch a_max")
        lower = minimum if a_min is None else a_min
        upper = maximum if a_max is None else a_max
        return original_clip(array, lower, upper)

    setattr(jnp, "clip", compatible_clip)
    return True


def _validate_inputs(freq_uHz: np.ndarray, power_psd: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    freq = np.asarray(freq_uHz, dtype=float)
    power = np.asarray(power_psd, dtype=float)
    if freq.ndim != 1 or power.ndim != 1 or freq.shape != power.shape:
        raise ValueError("Frequenz und PSD müssen gleich lange eindimensionale Arrays sein.")
    if len(freq) < 3 or not np.all(np.isfinite(freq)) or not np.all(np.isfinite(power)):
        raise ValueError("Frequenz und PSD müssen mindestens drei endliche Werte enthalten.")
    if not np.all(np.diff(freq) > 0) or np.any(power <= 0):
        raise ValueError("PBjam benötigt streng steigende Frequenzen und eine positive PSD.")
    return freq, power


def _cache_signature(
    freq: np.ndarray,
    power: np.ndarray,
    obs: dict[str, tuple[float, float]],
    n_orders: int,
) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(freq).tobytes())
    digest.update(np.ascontiguousarray(power).tobytes())
    digest.update(json.dumps(obs, sort_keys=True).encode("utf-8"))
    digest.update(str(n_orders).encode("ascii"))
    return digest.hexdigest()


def _cache_metadata(
    freq: np.ndarray,
    power: np.ndarray,
    obs: dict[str, tuple[float, float]],
    n_orders: int,
) -> dict[str, Any]:
    return {
        "frequency_sha256": hashlib.sha256(np.ascontiguousarray(freq).tobytes()).hexdigest(),
        "power_sha256": hashlib.sha256(np.ascontiguousarray(power).tobytes()).hexdigest(),
        "frequency_points": len(freq),
        "frequency_min": float(freq[0]),
        "frequency_max": float(freq[-1]),
        "obs": {key: list(value) for key, value in obs.items()},
        "n_orders": n_orders,
    }


def _cached_modes(cached: dict[str, Any]) -> dict[str, np.ndarray]:
    modes = cached.get("modes")
    if not isinstance(modes, dict) or any(key not in modes for key in BASE_MODE_KEYS):
        raise ValueError("Der PBjam-Cache enthält keine vollständigen Modenparameter.")
    arrays = {key: np.asarray(modes[key]) for key in BASE_MODE_KEYS}
    if len({len(value) for value in arrays.values()}) != 1:
        raise ValueError("Der PBjam-Cache enthält inkonsistente Modenarrays.")
    return detection_significance(arrays)


def _mode_orders(modeid_result: dict[str, Any], peakbag_result: dict[str, Any]) -> np.ndarray:
    peak_freq = np.asarray(peakbag_result["summary"]["freq"])[0]
    peak_l = np.asarray(peakbag_result["ell"], dtype=int).ravel()
    source_freq = np.asarray(modeid_result["summary"]["freq"])[0]
    source_l = np.asarray(modeid_result["ell"], dtype=int).ravel()
    source_n = np.asarray(modeid_result.get("enn", []), dtype=float).ravel()
    orders = np.full(len(peak_freq), -1, dtype=int)
    if len(source_n) != len(source_freq):
        return orders

    for index, (frequency, degree) in enumerate(zip(peak_freq, peak_l)):
        candidates = np.flatnonzero(source_l == degree)
        if len(candidates):
            nearest = candidates[np.argmin(np.abs(source_freq[candidates] - frequency))]
            if np.isfinite(source_n[nearest]):
                orders[index] = int(source_n[nearest])
    return orders


def extract_modes(
    modeid_result: dict[str, Any],
    peakbag_result: dict[str, Any],
    deltanu: float,
) -> dict[str, np.ndarray]:
    """Überführt PBjam-2.x-Ergebnisse in flache, pipelinefreundliche Arrays."""
    summary = peakbag_result["summary"]

    def values(key: str, row: int) -> np.ndarray:
        data = np.asarray(summary[key], dtype=float)
        if data.ndim != 2 or data.shape[0] != 2:
            raise ValueError(f"Unerwartete PBjam-Ausgabe für summary['{key}']: {data.shape}")
        return data[row]

    freq_err = values("freq", 1)
    modes = {
        "freq": values("freq", 0),
        "freq_err": freq_err,
        "l": np.asarray(peakbag_result["ell"], dtype=int).ravel(),
        "n": _mode_orders(modeid_result, peakbag_result),
        "height": values("height", 0),
        "height_err": values("height", 1),
        "width": values("width", 0),
        "width_err": values("width", 1),
        "quality": np.divide(
            0.03 * deltanu,
            freq_err,
            out=np.full_like(freq_err, np.nan),
            where=freq_err > 0,
        ),
    }
    lengths = {len(value) for value in modes.values()}
    if len(lengths) != 1:
        raise ValueError("PBjam lieferte Modenparameter mit inkonsistenten Längen.")
    return detection_significance(modes)


def detection_significance(
    modes: dict[str, np.ndarray],
    n_independent: int | None = None,
) -> dict[str, np.ndarray]:
    """Berechnet einen idealisierten FAP-Proxy aus PBjams Modenhöhe."""
    if n_independent is not None and n_independent < 1:
        raise ValueError("Die Zahl unabhängiger Frequenzbins muss positiv sein.")
    result = {key: np.asarray(value).copy() for key, value in modes.items()}
    lengths = {len(value) for value in result.values()}
    if len(lengths) != 1:
        raise ValueError("Modenarrays müssen gleich lang sein.")

    snr = np.asarray(result["height"], dtype=float)
    snr = np.where(np.isfinite(snr) & (snr >= 0), snr, 0.0)
    p_single = np.exp(-snr)
    if n_independent is None:
        fap = p_single
    else:
        fap = -np.expm1(n_independent * np.log1p(-p_single))
    result["height_snr_proxy"] = snr
    result["fap_proxy"] = np.clip(fap, 0.0, 1.0)
    return result


def global_false_alarm_probability(fap_proxy: float, n_modes_tested: int) -> float:
    """Nominale globale Rate für n unabhängige Tests desselben FAP-Proxys."""
    if not 0 <= fap_proxy <= 1 or n_modes_tested < 0:
        raise ValueError("FAP-Proxy muss in [0, 1] liegen und die Testzahl nichtnegativ sein.")
    return float(-np.expm1(n_modes_tested * np.log1p(-fap_proxy)))


def estimate_d02_fraction(
    modes: dict[str, np.ndarray],
    dnu: float,
    anchor_mask: np.ndarray,
) -> tuple[float, str]:
    """Bestimmt δν₀₂/Δν aus Gold-Paaren oder einer empirischen Δν-Relation."""
    if dnu <= 0:
        raise ValueError("Δν muss positiv sein.")
    frequency = np.asarray(modes["freq"], dtype=float)
    degree = np.asarray(modes["l"], dtype=int)
    order = np.asarray(modes.get("n", np.full(len(frequency), -1)), dtype=int)
    anchor = np.asarray(anchor_mask, dtype=bool)
    if len(anchor) != len(frequency):
        raise ValueError("Der Anker-Mask muss dieselbe Länge wie die Modenarrays haben.")

    fractions: list[float] = []
    for radial_index in np.flatnonzero(anchor & (degree == 0) & (order >= 0)):
        quadrupole = np.flatnonzero(
            anchor & (degree == 2) & (order == order[radial_index] - 1)
        )
        if len(quadrupole) == 1:
            separation = frequency[radial_index] - frequency[int(quadrupole[0])]
            if 0 < separation < 0.2 * dnu:
                fractions.append(float(separation / dnu))
    if fractions:
        return float(np.median(fractions)), "gold_pairs"
    return float(-0.0324 * np.log10(dnu) + 0.1388), "dnu_relation"


def _sequence_members(
    frequency: np.ndarray,
    candidate_mask: np.ndarray,
    dnu: float,
    tolerance_fraction: float,
    minimum_length: int,
) -> np.ndarray:
    """Markiert Kandidaten in Läufen konsekutiver radialer Ordnungen."""
    members = np.zeros(len(frequency), dtype=bool)
    indices = np.flatnonzero(candidate_mask)
    if len(indices) < minimum_length:
        return members
    indices = indices[np.argsort(frequency[indices])]
    spacings_ok = np.abs(np.diff(frequency[indices]) - dnu) <= tolerance_fraction * dnu
    start = 0
    for boundary in np.flatnonzero(~spacings_ok):
        stop = int(boundary) + 1
        if stop - start >= minimum_length:
            members[indices[start:stop]] = True
        start = stop
    if len(indices) - start >= minimum_length:
        members[indices[start:]] = True
    return members


def tier_modes(
    modes: dict[str, np.ndarray],
    fap_gold: float = 0.01,
    fap_silver: float = 0.1,
    quality_min: float | None = None,
    dnu: float | None = None,
    ridge_tol_uHz: float = 1.5,
    d02_fraction: float | None = None,
    sequence_tolerance: float = 0.10,
    sequence_minimum: int = 3,
) -> dict[str, np.ndarray]:
    """Kennzeichnet Moden als Gold, Silber oder Rest."""
    if not 0 <= fap_gold <= fap_silver <= 1:
        raise ValueError("FAP-Schwellen müssen 0 ≤ Gold ≤ Silber ≤ 1 erfüllen.")
    if dnu is not None and (dnu <= 0 or ridge_tol_uHz <= 0):
        raise ValueError("Ridge-Test benötigt positive Δν und Toleranz.")
    if d02_fraction is not None and not 0 < d02_fraction < 0.5:
        raise ValueError("δν₀₂/Δν muss zwischen 0 und 0.5 liegen.")
    if not 0 < sequence_tolerance < 1 or sequence_minimum < 3:
        raise ValueError("Sequenztest benötigt 0 < Toleranz < 1 und mindestens drei Moden.")
    result = detection_significance(modes)
    fap = np.asarray(result["fap_proxy"], dtype=float)
    gold = fap <= fap_gold
    if quality_min is not None:
        quality = np.asarray(result["quality"], dtype=float)
        gold |= np.isfinite(quality) & (quality >= quality_min)
    tier = np.full(len(fap), "rest", dtype="<U6")
    tier[(fap <= fap_silver) & ~gold] = "silver"
    tier[gold] = "gold"
    ridge_deviation = np.full(len(fap), np.nan)
    sequence_promoted = np.zeros(len(fap), dtype=bool)
    d02_fraction_used = np.nan
    d02_fraction_source = "none"
    if dnu is not None:
        degree = np.asarray(result["l"], dtype=int)
        frequency = np.asarray(result["freq"], dtype=float)
        radial_phases = (frequency[gold & (degree == 0)] % dnu) / dnu
        if len(radial_phases) >= 2:
            radial_phase = float(np.median(radial_phases))
            if d02_fraction is None:
                d02_fraction_used, d02_fraction_source = estimate_d02_fraction(result, dnu, gold)
            else:
                d02_fraction_used, d02_fraction_source = d02_fraction, "configured"
            expected_phases = {0: radial_phase, 2: (radial_phase - d02_fraction_used) % 1.0}
            for mode_degree, expected_phase in expected_phases.items():
                selected = degree == mode_degree
                phase_offset = np.abs((frequency[selected] % dnu) / dnu - expected_phase)
                ridge_deviation[selected] = np.minimum(phase_offset, 1.0 - phase_offset) * dnu
            ridge_silver = (
                (tier == "rest")
                & np.isin(degree, (0, 2))
                & (ridge_deviation <= ridge_tol_uHz)
            )
            tier[ridge_silver] = "silver"
            for mode_degree in (0, 2):
                sequence_promoted |= _sequence_members(
                    frequency,
                    (tier != "rest") & (degree == mode_degree),
                    dnu,
                    sequence_tolerance,
                    sequence_minimum,
                ) & (tier == "silver")
            tier[sequence_promoted] = "gold"
    result["tier"] = tier
    result["ridge_dev"] = ridge_deviation
    result["sequence_promoted"] = sequence_promoted
    result["d02_fraction_used"] = np.full(len(fap), d02_fraction_used)
    result["d02_fraction_source"] = np.full(len(fap), d02_fraction_source, dtype="<U16")
    return result


def select_mode_tiers(
    modes: dict[str, np.ndarray], tiers: tuple[str, ...] = ("gold",)
) -> dict[str, np.ndarray]:
    """Wählt vollständige Modenarrays anhand ihrer Konfidenzstufe aus."""
    if "tier" not in modes:
        raise ValueError("Moden müssen vor der Auswahl mit tier_modes klassifiziert werden.")
    mask = np.isin(np.asarray(modes["tier"]), tiers)
    return {key: np.asarray(value)[mask] for key, value in modes.items()}


def adaptive_quality_threshold(
    modes: dict[str, np.ndarray],
    factor: float = 1.10,
    floor: float = 0.75,
    ceiling: float = 2.0,
) -> float:
    """Leitet den Quality-Cut aus dem Median der endlichen Modenwerte ab."""
    quality = np.asarray(modes["quality"], dtype=float)
    finite_quality = quality[np.isfinite(quality)]
    if not len(finite_quality):
        return floor
    return float(np.clip(np.median(finite_quality) * factor, floor, ceiling))


def filter_reliable_modes(
    modes: dict[str, np.ndarray],
    quality_min: float | None = None,
    height_min: float = 1.0,
) -> dict[str, np.ndarray]:
    """Behält Moden mit informativem Frequenzposterior und messbarer Höhe."""
    if quality_min is None:
        quality_min = adaptive_quality_threshold(modes)
    quality = np.asarray(modes["quality"], dtype=float)
    height = np.asarray(modes["height"], dtype=float)
    mask = (quality >= quality_min) & (height >= height_min)
    return {key: np.asarray(value)[mask] for key, value in modes.items()}


def find_avoided_crossings(
    freq_l1: np.ndarray, deltanu: float, threshold: float = 0.25
) -> dict[str, np.ndarray]:
    """Findet durch gestörte l=1-Frequenzabstände angezeigte avoided crossings."""
    frequencies = np.sort(np.asarray(freq_l1, dtype=float))
    if deltanu <= 0:
        raise ValueError("Δν muss positiv sein.")
    spacings = np.diff(frequencies)
    deviations = np.abs(spacings - deltanu) / deltanu
    indices = np.flatnonzero(deviations > threshold)
    return {
        "crossings": (frequencies[indices] + frequencies[indices + 1]) / 2.0,
        "spacings": spacings,
        "deviations": deviations,
        "freq_l1": frequencies,
    }


def _write_outputs(
    modes: dict[str, np.ndarray],
    cache_path: Path,
    table_path: Path,
    signature: str,
    metadata: dict[str, Any],
) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    table_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(
        json.dumps(
            {
                "schema": CACHE_SCHEMA_VERSION,
                "signature": signature,
                "inputs": metadata,
                "modes": {key: value.tolist() for key, value in modes.items()},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(MODE_KEYS)
        writer.writerows(zip(*(modes[key] for key in MODE_KEYS)))


def run_pbjam_modeid(
    freq_uHz: np.ndarray,
    power_psd: np.ndarray,
    numax: tuple[float, float],
    deltanu: tuple[float, float],
    teff: tuple[float, float],
    bp_rp: tuple[float, float] | None = None,
    star_id: str = "star",
    outdir: str | Path = "results/pbjam",
    n_orders: int = 7,
    use_cache: bool = True,
    force: bool = False,
    reuse_existing_cache: bool = False,
) -> dict[str, np.ndarray]:
    """Führt PBjam ModeID und Peakbagging auf einer Parseval-normierten PSD aus."""
    freq, power = _validate_inputs(freq_uHz, power_psd)
    obs = {"numax": tuple(numax), "dnu": tuple(deltanu), "teff": tuple(teff)}
    if bp_rp is not None:
        obs["bp_rp"] = tuple(bp_rp)
    for key, (value, uncertainty) in obs.items():
        if not np.isfinite(value) or not np.isfinite(uncertainty) or uncertainty <= 0:
            raise ValueError(f"{key} benötigt einen endlichen Wert und eine positive Unsicherheit.")

    output_dir = Path(outdir)
    stem = _safe_name(star_id)
    cache_path = output_dir / "cache" / f"{stem}_modes.json"
    table_path = Path("results/tables") / f"{stem}_pbjam_modes.csv"
    signature = _cache_signature(freq, power, obs, n_orders)
    if use_cache and not force and cache_path.exists():
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
        if cached.get("schema") in (1, 2, 3, CACHE_SCHEMA_VERSION) and cached.get("signature") == signature:
            print(f"  PBjam-Cache geladen: {cache_path}")
            return _cached_modes(cached)
        if reuse_existing_cache:
            print(
                f"  WARNUNG: Vorhandener PBjam-Cache trotz abweichender Signatur übernommen: {cache_path}\n"
                "           Nur für die kontrollierte Migration bereits geprüfter Ergebnisse verwenden."
            )
            return _cached_modes(cached)
        print(
            f"  PBjam-Cache nicht verwendet: Signatur stimmt nicht überein ({cache_path}).\n"
            "  Für die kontrollierte Übernahme eines geprüften Alt-Caches: --reuse-existing-pbjam"
        )

    installed_version = _pbjam_version()
    numpy_compatibility = _apply_numpy_compatibility()
    jax_compatibility = _apply_jax_compatibility()
    import pbjam

    if numpy_compatibility:
        print("  NumPy-Kompatibilität für PBjam 2.0 aktiviert (trapz → trapezoid).")
    if jax_compatibility:
        print("  JAX-Kompatibilität für PBjam 2.0 aktiviert (clip-Schlüssel).")
    print(f"  Starte PBjam {installed_version} (ModeID + Peakbagging) ...")
    target = pbjam.star(star_id, freq, power, obs, outpath=output_dir, N_p=n_orders)
    modeid_result, peakbag_result = target()
    modes = extract_modes(modeid_result, peakbag_result, deltanu=deltanu[0])
    _write_outputs(modes, cache_path, table_path, signature, _cache_metadata(freq, power, obs, n_orders))
    print(f"  PBjam-Tabelle gespeichert: {table_path}")
    return modes