"""Abgeleitete Messgrößen für den seismischen Atlas."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class D02Measurement:
    value: float | None
    uncertainty: float | None
    pairs: int
    values: np.ndarray
    uncertainties: np.ndarray


def measure_d02(
    modes: dict[str, np.ndarray], quality_min: float = 0.65
) -> D02Measurement:
    """Misst delta-nu-02 mit einer lockeren, von der Hauptauswahl getrennten Schwelle."""
    frequency = np.asarray(modes["freq"], dtype=float)
    frequency_error = np.asarray(modes["freq_err"], dtype=float)
    degree = np.asarray(modes["l"], dtype=int)
    order = np.asarray(modes["n"], dtype=int)
    if len({len(frequency), len(frequency_error), len(degree), len(order)}) != 1:
        raise ValueError("Modenarrays müssen gleich lang sein.")
    if "quality" in modes:
        quality = np.asarray(modes["quality"], dtype=float)
        if len(quality) != len(frequency):
            raise ValueError("Modenarrays müssen gleich lang sein.")
        selected = np.isfinite(quality) & (quality >= quality_min)
        frequency = frequency[selected]
        frequency_error = frequency_error[selected]
        degree = degree[selected]
        order = order[selected]

    values: list[float] = []
    uncertainties: list[float] = []
    radial_indices = np.flatnonzero((degree == 0) & (order >= 0))
    for radial_index in radial_indices:
        quadrupole_indices = np.flatnonzero((degree == 2) & (order == order[radial_index] - 1))
        if len(quadrupole_indices) != 1:
            continue
        quadrupole_index = int(quadrupole_indices[0])
        error = float(np.hypot(frequency_error[radial_index], frequency_error[quadrupole_index]))
        separation = float(frequency[radial_index] - frequency[quadrupole_index])
        if separation > 0 and np.isfinite(error) and error > 0:
            values.append(separation)
            uncertainties.append(error)

    value_array = np.asarray(values, dtype=float)
    uncertainty_array = np.asarray(uncertainties, dtype=float)
    if not len(value_array):
        return D02Measurement(None, None, 0, value_array, uncertainty_array)

    mean = float(np.mean(value_array))
    scatter = float(np.std(value_array)) if len(value_array) > 1 else float(uncertainty_array[0])
    return D02Measurement(
        mean,
        scatter,
        len(value_array),
        value_array,
        uncertainty_array,
    )


def crossing_candidates(
    freq_l1: np.ndarray, deltanu: float, threshold: float = 0.25
) -> list[dict[str, float | int]]:
    """Liefert QC-fähige Kandidaten aus gestörten benachbarten l=1-Abständen."""
    frequencies = np.sort(np.asarray(freq_l1, dtype=float))
    if deltanu <= 0 or threshold <= 0:
        raise ValueError("Delta-nu und Crossing-Schwelle müssen positiv sein.")
    if len(frequencies) < 2:
        return []

    spacings = np.diff(frequencies)
    deviations = np.abs(spacings - deltanu) / deltanu
    indices = np.flatnonzero(deviations > threshold)
    return [
        {
            "lower_freq": float(frequencies[index]),
            "upper_freq": float(frequencies[index + 1]),
            "midpoint": float((frequencies[index] + frequencies[index + 1]) / 2.0),
            "spacing": float(spacings[index]),
            "relative_deviation": float(deviations[index]),
            "interval_index": int(index),
        }
        for index in indices
    ]
