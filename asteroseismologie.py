#!/usr/bin/env python3
"""
Asteroseismologie-Workflow: Solar-like Oscillators mit TESS und Kepler
=======================================================================
Vollständiger Pipeline von der Lichtkurve bis zu abgeleiteten
Stellar-Parametern (M, R, log g, L) via Skalenrelationen.

Unterstützte Archive (auto-detektiert aus dem Ziel-Präfix):
  TIC …  → TESS SPOC  (2 min, Nyquist ≈ 4167 μHz)
  KIC …  → Kepler LC  (30 min, Nyquist ≈  278 μHz)
  EPIC … → K2 LC      (30 min, Nyquist ≈  278 μHz)

νmax wird vollautomatisch bestimmt — kein Schätzwert erforderlich:
  [1] Harvey-Hintergrundfit im Log-Raum
  [2] SNR-Spektrum (PSD / Hintergrund)
  [3] Erster νmax-Kandidat aus geglättetem SNR (oberhalb Harvey-Knie)
  [4] Gauß-Fit zur finalen νmax-Bestimmung

Verwendung
----------
    uv run asteroseismologie.py                                   # η Serpentis (Default)
    uv run asteroseismologie.py --tic "KIC 4351319" --teff 4800   # Kepler-Roter-Riese
    uv run asteroseismologie.py --tic "TIC 272821450" --teff 4900  # ε Ophiuchi
    uv run asteroseismologie.py --oversample 2                     # schneller Test
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

# UTF-8-Ausgabe erzwingen — verhindert UnicodeEncodeError auf Windows (cp1252)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

# ---------------------------------------------------------------------------
# Solarkonstanten (Referenzwerte)
# ---------------------------------------------------------------------------
SUN_NUMAX   = 3090.0   # μHz   (Chaplin & Miglio 2013)
SUN_DELTANU = 135.1    # μHz
SUN_TEFF    = 5778.0   # K

# ---------------------------------------------------------------------------
# Standardstern: η Serpentis
# ---------------------------------------------------------------------------
DEFAULT_TIC  = "TIC 234295610"
DEFAULT_NAME = "η Serpentis"
DEFAULT_TEFF = 4970.0            # K  (SIMBAD)

FMIN_UHZ   = 1.0
FMAX_UHZ   = 300.0
OVERSAMPLE = 5


# ===========================================================================
# Schritt 1 — Lichtkurve laden und bereinigen
# ===========================================================================

def _mission_defaults(target_id: str) -> tuple[str, int, float]:
    """
    Liefert (author, exptime_s, min_span_days) passend zum Zielidentifikator.

    KIC  → Kepler Long Cadence (1800 s); Q0 (~9d) und Q1 (~33d) werden mit
           min_span=60d automatisch übersprungen, Q2+ (~90d) werden genutzt.
    EPIC → K2 Long Cadence    (1800 s); Kampagnen sind ~80 d, min_span=60d.
    sonst → TESS SPOC 2-min   (120 s);  Sektoren sind ~27 d, min_span=20d.
    """
    tid = target_id.strip().upper()
    if tid.startswith("KIC"):
        return "Kepler", 1800, 60.0
    if tid.startswith("EPIC"):
        return "K2", 1800, 60.0
    return "SPOC", 120, 20.0


def load_lightcurve(
    target_id: str,
    max_sectors: int = 1,
    author: str | None = None,
    exptime: int | None = None,
) -> lk.LightCurve:
    """
    Lädt Lichtkurven vom MAST-Archiv (TESS, Kepler oder K2).

    Mission wird automatisch aus dem Ziel-Präfix (TIC/KIC/EPIC) ermittelt;
    --author und --exptime überschreiben die Auto-Erkennung.
    Gibt eine in ppm normierte, zusammengeführte Lichtkurve zurück.
    """
    _author, _exptime, min_span = _mission_defaults(target_id)
    author  = author  or _author
    exptime = exptime or _exptime

    search = lk.search_lightcurve(target_id, author=author, exptime=exptime)
    if len(search) == 0:
        raise ValueError(
            f"Keine Daten für '{target_id}' gefunden "
            f"(author='{author}', exptime={exptime} s).\n"
            f"  Tipp: --author Kepler --exptime 1800  für Kepler-Sterne\n"
            f"        --author SPOC   --exptime 120    für TESS-Sterne"
        )

    # Lade etwas mehr als max_sectors, um sehr kurze Quartale (z. B. Kepler Q0
    # mit nur ~9 d) herausfiltern zu können, ohne zu wenig Daten zu erhalten.
    n_try = min(len(search), max_sectors + 3)
    print(f"  Gefunden: {len(search)} Sektor(en)/Quartale, lade Kandidaten ...")
    lc_coll = search[:n_try].download_all(quality_bitmask="hardest")

    lcs = []
    skipped = 0
    for lc in lc_coll:
        lc_clean = lc.remove_nans().remove_outliers(sigma=4.0)
        if len(lc_clean) == 0:
            skipped += 1
            continue
        span = float(lc_clean.time.value[-1] - lc_clean.time.value[0])
        if span < min_span and len(lcs) == 0 and (n_try - skipped) > max_sectors:
            # Zu kurz und wir haben noch Reservekandidaten — überspringen
            skipped += 1
            print(f"  Überspringe kurzes Segment ({span:.1f} d < {min_span:.0f} d)")
            continue
        median   = float(lc_clean.flux.value.mean())
        flux_ppm = (lc_clean.flux.value / median - 1.0) * 1e6
        lcs.append(lk.LightCurve(time=lc_clean.time, flux=flux_ppm))
        if len(lcs) >= max_sectors:
            break

    if not lcs:
        raise RuntimeError("Alle heruntergeladenen Sektoren sind leer oder zu kurz.")

    print(f"  Verwende {len(lcs)} Segment(e)")
    if len(lcs) > 1:
        return lk.LightCurveCollection(lcs).stitch()
    return lcs[0]


# ===========================================================================
# Schritt 2 — Powerspektrum (Lomb-Scargle, einseitig, Parseval-normiert)
# ===========================================================================

def compute_power_spectrum(
    lc: lk.LightCurve,
    fmin_uHz: float = FMIN_UHZ,
    fmax_uHz: float = FMAX_UHZ,
    oversample: int = OVERSAMPLE,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Berechnet die einseitige PSD in ppm²/μHz mittels Lomb-Scargle.

    Frequenzauflösung: δf = 1 / (oversample · T_gesamt)
    fmax wird automatisch auf 90 % der Nyquist-Frequenz des Datensatzes
    gecappt — verhindert Aliasing bei Kepler Long Cadence (Nyquist ≈ 278 μHz).
    """
    time_s   = lc.time.value * 86400.0
    flux_ppm = lc.flux.value

    dt_s         = float(np.median(np.diff(time_s)))
    nyquist_uHz  = 1e6 / (2.0 * dt_s)
    fmax_uHz     = min(fmax_uHz, 0.9 * nyquist_uHz)

    T_s      = time_s[-1] - time_s[0]
    df_uHz   = 1e6 / (oversample * T_s)
    freq_uHz = np.arange(fmin_uHz, fmax_uHz, df_uHz)
    freq_Hz  = freq_uHz * 1e-6

    ls    = LombScargle(time_s, flux_ppm, normalization="psd")
    power = ls.power(freq_Hz)

    # Einseitige PSD, Parseval-normiert
    return freq_uHz, power * 2.0 / df_uHz


# ===========================================================================
# Schritt 3 — νmax ohne Schätzwert: 4-stufige Bootstrapping-Pipeline
# ===========================================================================

def _harvey_model(f: np.ndarray, a: float, b: float, c: float, w: float) -> np.ndarray:
    """Harvey-Hintergrundmodell: Granulation + Weißrauschen."""
    return a / (1.0 + (f / b) ** c) + w


def fit_harvey_background(
    freq: np.ndarray, power: np.ndarray
) -> tuple[np.ndarray, tuple[float, float, float, float] | None]:
    """
    Passt das Harvey-Modell im Log-Raum an das Powerspektrum.

    Fit erfolgt nur auf dem unteren 10 % des Frequenzbereichs (mindestens
    50 Bins), wo die Granulation dominiert und der Oszillations-Buckel
    noch nicht stört. Das gefittete Modell wird auf den vollen Bereich
    extrapoliert. Damit wird b nicht durch den Oszillationsexzess verzerrt.
    """
    log_power = np.log10(np.clip(power, 1e-10, None))

    def log_harvey(f, log_a, log_b, c, log_w):
        return np.log10(_harvey_model(f, 10**log_a, 10**log_b, c, 10**log_w))

    # Fit-Bereich: 2 %–15 % von fmax.
    # Untergrenze (2 %) überspringt DC-Spike und sehr niederfrequentes Rot-Rauschen.
    # Obergrenze (15 %) liegt für typische Rote Riesen noch unterhalb des
    # Oszillations-Buckels (νmax > 15 % × fmax bei fmax = 300 μHz → νmax > 45 μHz).
    i_lo = max(0, np.searchsorted(freq, freq[-1] * 0.02))
    i_hi = np.searchsorted(freq, freq[-1] * 0.15)
    if i_hi - i_lo < 30:           # Fallback wenn Bereich zu schmal
        i_lo = 0
        i_hi = max(50, int(len(freq) * 0.10))
    f_bg   = freq[i_lo:i_hi]
    lp_bg  = log_power[i_lo:i_hi]

    # Datengetriebene Bounds (über gesamtes Spektrum für Extrapolations-Stabilität)
    lp_lo = log_power.min() - 1.0
    lp_hi = log_power.max() + 1.0
    lf_lo = np.log10(freq[0])
    lf_hi = np.log10(freq[-1])
    bounds_lo = [lp_lo, lf_lo, 0.5, lp_lo]
    bounds_hi = [lp_hi, lf_hi, 6.0, lp_hi]

    # Weißrauschen aus dem hohen Frequenzende schätzen
    w_guess = float(np.median(power[int(len(power) * 0.8):]))

    p0_raw = [
        np.log10(max(power[0], 1e-10)),          # log_a: Granulationsamplitude
        np.log10(freq[len(freq) // 10]),          # log_b: Knie bei ~10 % von fmax
        2.5,                                       # c: Exponent
        np.log10(max(w_guess, 1e-10)),            # log_w: Weißrauschen
    ]
    eps = 1e-6
    p0 = [float(np.clip(v, bounds_lo[i] + eps, bounds_hi[i] - eps))
          for i, v in enumerate(p0_raw)]

    try:
        popt, _ = curve_fit(
            log_harvey, f_bg, lp_bg,
            p0=p0, bounds=(bounds_lo, bounds_hi), maxfev=10_000
        )
        a, b, c, w = 10**popt[0], 10**popt[1], popt[2], 10**popt[3]
        # Extrapolation auf den vollen Frequenzbereich
        return _harvey_model(freq, a, b, c, w), (a, b, c, w)
    except RuntimeError:
        sigma    = max(1, len(freq) // 10)
        fallback = gaussian_filter1d(power, sigma=sigma)
        return fallback, None


def _compute_snr_spectrum(power: np.ndarray, bg_model: np.ndarray) -> np.ndarray:
    bg_safe = np.clip(bg_model, np.percentile(bg_model, 1), None)
    return power / bg_safe


def _find_numax_candidate(
    freq: np.ndarray,
    snr: np.ndarray,
    smooth_width_uHz: float = 15.0,
    fmin_search: float = 5.0,
) -> tuple[float, np.ndarray]:
    """
    Erster νmax-Kandidat aus dem breitbandig geglätteten SNR-Spektrum.

    15 μHz Glättungsbreite ist frequenzunabhängig und funktioniert für
    Rote Riesen (νmax 5–300 μHz) ohne jeden Prior.
    """
    df         = freq[1] - freq[0]
    sigma_bins = max(1, round(smooth_width_uHz / df))
    snr_smooth = gaussian_filter1d(snr, sigma=sigma_bins)

    mask = freq >= fmin_search
    idx  = np.argmax(snr_smooth[mask])
    return float(freq[mask][idx]), snr_smooth


def _refine_numax(
    freq: np.ndarray,
    snr_smooth: np.ndarray,
    numax_cand: float,
    band_lo_factor: float = 0.55,
    band_hi_factor: float = 1.55,
) -> tuple[float, float]:
    """
    Verfeinert den Kandidaten durch Gauß-Fit auf einem automatischen Bandpass.

    Faktoren 0.55 / 1.55 sind asymmetrisch, weil die Oszillationshülle zur
    niederfrequenten Seite steiler abfällt.
    """
    band_lo = numax_cand * band_lo_factor
    band_hi = numax_cand * band_hi_factor
    mask    = (freq >= band_lo) & (freq <= band_hi)
    f_fit   = freq[mask]
    s_fit   = snr_smooth[mask]

    def gauss_plus_bg(f, A, mu, sigma, c0):
        return A * np.exp(-0.5 * ((f - mu) / sigma) ** 2) + c0

    p0     = [s_fit.max() - s_fit.min(), numax_cand, numax_cand * 0.25, s_fit.min()]
    bounds = ([0, band_lo, 1.0, 0], [1e6, band_hi, band_hi - band_lo, np.inf])

    try:
        popt, _ = curve_fit(gauss_plus_bg, f_fit, s_fit, p0=p0, bounds=bounds, maxfev=5_000)
        numax_fit = float(popt[1])
        sigma_fit = float(abs(popt[2]))
        # Fit am Bandpassrand → kein echtes Maximum gefunden → Kandidat beibehalten
        margin = 0.05 * (band_hi - band_lo)
        if numax_fit < band_lo + margin or numax_fit > band_hi - margin:
            return numax_cand, numax_cand * 0.25
        return numax_fit, sigma_fit
    except RuntimeError:
        return numax_cand, numax_cand * 0.25


def estimate_numax_auto(
    freq: np.ndarray, power: np.ndarray, verbose: bool = True
) -> dict:
    """
    Vollautomatische νmax-Bestimmung ohne Schätzwert.

    Gibt dict zurück mit: numax, numax_sigma, numax_cand,
    harvey_bg, snr, snr_smooth, harvey_params.
    """
    # 3a: Harvey-Hintergrundmodell
    bg_model, harvey_params = fit_harvey_background(freq, power)
    if verbose and harvey_params:
        a, b, c, w = harvey_params
        print(f"  Harvey-Fit : a={a:.0f} ppm²/μHz  b={b:.1f} μHz  "
              f"c={c:.2f}  w={w:.1f} ppm²/μHz")

    # 3b: SNR-Spektrum
    snr = _compute_snr_spectrum(power, bg_model)

    # 3c: Erster Kandidat — Suche startet oberhalb des Harvey-Knies,
    # damit die Granulations-dominierte Region nicht als Kandidat landet.
    if harvey_params is not None:
        b_harvey = harvey_params[1]
        fmin_candidate = float(np.clip(3.0 * b_harvey, 5.0, freq[-1] * 0.5))
    else:
        fmin_candidate = 5.0
    numax_cand, snr_smooth = _find_numax_candidate(freq, snr, fmin_search=fmin_candidate)
    if verbose:
        print(f"  νmax-Kandidat (SNR-Maximum): {numax_cand:.1f} μHz  "
              f"[Suche ab {fmin_candidate:.1f} μHz]")

    # 3d: Gauß-Fit
    numax, numax_sigma = _refine_numax(freq, snr_smooth, numax_cand)
    if verbose:
        print(f"  νmax final  (Gauß-Fit):      {numax:.1f} μHz  "
              f"(σ = {numax_sigma:.1f} μHz)")

    return {
        "numax":          numax,
        "numax_sigma":    numax_sigma,
        "numax_cand":     numax_cand,
        "harvey_bg":      bg_model,
        "snr":            snr,
        "snr_smooth":     snr_smooth,
        "harvey_params":  harvey_params,
    }


# ===========================================================================
# Schritt 4 — Δν via Autokorrelation
# ===========================================================================

def estimate_deltanu(
    freq: np.ndarray, power: np.ndarray, numax: float, verbose: bool = True
) -> float:
    """
    Schätzt Δν aus der Autokorrelation des bandgefilterten Powerspektrums.

    Suchfenster eng um den empirischen Prior zentriert (±20 %):
        Δν_prior ≈ 0.263 · νmax^0.772  (Stello et al. 2009)

    Enger Prior verhindert, dass die ACF Subharmonische oder Aliase
    (z. B. Δν/√2) als Δν zurückgibt — ein häufiger Fehler bei niedrigem SNR.
    Liegt kein ACF-Maximum im validierten Fenster, wird der Prior verwendet.
    """
    deltanu_prior = 0.263 * numax ** 0.772

    # Bandpass um νmax für die Autokorrelation
    width = max(4 * deltanu_prior, numax * 0.5)
    mask  = (freq > numax - width) & (freq < numax + width)
    p_bp  = power[mask] - power[mask].mean()

    # Autokorrelation
    acf      = np.correlate(p_bp, p_bp, mode="full")
    acf      = acf[len(acf) // 2:]
    lag_freq = np.arange(len(acf)) * (freq[1] - freq[0])

    # Enges Suchfenster: ±20 % um Skalenrelations-Prior
    lo = 0.8 * deltanu_prior
    hi = 1.2 * deltanu_prior
    mask_lag = (lag_freq >= lo) & (lag_freq <= hi)

    if verbose:
        print(f"  Δν (Stello-Prior):   {deltanu_prior:.2f} μHz  "
              f"→ Suchfenster [{lo:.2f}, {hi:.2f}] μHz")

    if not np.any(mask_lag):
        if verbose:
            print(f"  Δν (ACF):            kein Peak im Fenster — Prior verwendet")
        return deltanu_prior

    deltanu_acf = float(lag_freq[mask_lag][np.argmax(acf[mask_lag])])

    if verbose:
        abw = abs(deltanu_acf - deltanu_prior) / deltanu_prior * 100
        flag = "  ✓" if abw < 15 else f"  ⚠ {abw:.0f} % vom Prior"
        print(f"  Δν (ACF):            {deltanu_acf:.2f} μHz{flag}")

    # Échelle-Kohärenz-Verfeinerung: feines Δν-Gitter, maximiere Spaltenvarianz
    deltanu_final = _refine_deltanu_echelle(freq, power, numax, deltanu_acf,
                                            verbose=verbose)
    return deltanu_final


def _echelle_column_variance(
    freq: np.ndarray, power: np.ndarray, numax: float,
    deltanu: float, n_cols: int = 40
) -> float:
    """
    Varianz des gefalteten Spaltenprofils — Kohärenzmaß für das Échelle.

    Maximiert wenn Moden vertikal übereinander stehen (korrektes Δν).
    """
    flo  = max(numax - 5.0 * deltanu, freq[0])
    fhi  = min(numax + 5.0 * deltanu, freq[-1])
    mask = (freq >= flo) & (freq <= fhi)
    if mask.sum() < 10:
        return 0.0
    x         = freq[mask] % deltanu
    bins       = np.linspace(0, deltanu, n_cols + 1)
    wsum, _    = np.histogram(x, bins=bins, weights=power[mask])
    counts, _  = np.histogram(x, bins=bins)
    with np.errstate(divide="ignore", invalid="ignore"):
        profile = np.where(counts > 0, wsum / counts, 0.0)
    return float(np.var(profile))


def _refine_deltanu_echelle(
    freq: np.ndarray, power: np.ndarray, numax: float, deltanu_acf: float,
    n_trial: int = 300, search_width: float = 0.15, verbose: bool = True,
) -> float:
    """
    Verfeinert Δν durch Maximierung der Échelle-Kohärenz (Spaltenvarianz).

    Sucht auf einem feinen Gitter von ±15 % um den ACF-Wert. Das optimale
    Δν lässt die Oszillations-Moden in vertikalen Säulen stehen.
    """
    dn_lo   = deltanu_acf * (1.0 - search_width)
    dn_hi   = deltanu_acf * (1.0 + search_width)
    dn_grid = np.linspace(dn_lo, dn_hi, n_trial)

    coherence = np.array([
        _echelle_column_variance(freq, power, numax, dn) for dn in dn_grid
    ])
    deltanu_opt = float(dn_grid[np.argmax(coherence)])

    if verbose:
        print(f"  Δν (Échelle-Kohärenz): {deltanu_opt:.2f} μHz")

    return deltanu_opt


# ===========================================================================
# Schritt 5 — Skalenrelationen
# ===========================================================================

def scaling_relations(numax: float, deltanu: float, teff: float) -> dict:
    """
    Ableitung von Stellar-Parametern aus νmax, Δν und Teff.

    Chaplin & Miglio (2013, ARA&A 51, 353).
    Genauigkeit: ~5 % in R, ~10–15 % in M für Rote Riesen.
    """
    r_nu = numax   / SUN_NUMAX
    r_dn = deltanu / SUN_DELTANU
    r_T  = teff    / SUN_TEFF

    mass   = r_nu ** 3 * r_dn ** -4 * r_T **  1.5
    radius = r_nu      * r_dn ** -2 * r_T **  0.5
    logg   = 4.44 + np.log10(r_nu) + 0.5 * np.log10(r_T)
    lumin  = radius ** 2 * r_T ** 4

    return {"mass": mass, "radius": radius, "logg": logg, "lumin": lumin}


# ===========================================================================
# Schritt 6 — Abbildungen
# ===========================================================================

def make_figure(
    lc: lk.LightCurve,
    freq: np.ndarray,
    power: np.ndarray,
    numax_result: dict,
    deltanu: float,
    params: dict,
    star_name: str,
    outpath: Path,
) -> None:
    """
    Erstellt eine 4-Panel-Abbildung:
      oben links  — Lichtkurve
      oben rechts — PSD + Harvey-Hintergrund
      unten links — SNR-Spektrum
      unten rechts — Échelle-Diagramm mit Stellar-Parametern
    """
    numax      = numax_result["numax"]
    bg         = numax_result["harvey_bg"]
    snr_smooth = numax_result["snr_smooth"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Asteroseismologie: {star_name}", fontsize=13, fontweight="bold")

    # --- Panel 1: Lichtkurve ---
    ax = axes[0, 0]
    t0 = lc.time.value[0]
    ax.plot(lc.time.value - t0, lc.flux.value, "k-", lw=0.3, alpha=0.6, rasterized=True)
    ax.set_xlabel("Zeit (d)")
    ax.set_ylabel("Flux (ppm)")
    ax.set_title("Lichtkurve")

    # --- Panel 2: PSD + Harvey ---
    ax = axes[0, 1]
    ax.semilogy(freq, power, color="0.65", lw=0.5, label="PSD", rasterized=True)
    ax.semilogy(freq, bg,    color="C1",   lw=1.8, label="Harvey-Hintergrund")
    ax.axvline(numax, color="C0", lw=1.5, ls="--",
               label=f"νmax = {numax:.1f} μHz")
    ax.set_xlabel("Frequenz (μHz)")
    ax.set_ylabel("PSD (ppm²/μHz)")
    ax.set_title("Powerspektrum")
    ax.legend(fontsize=8, loc="upper right")

    # --- Panel 3: SNR-Spektrum ---
    ax = axes[1, 0]
    ax.plot(freq, numax_result["snr"], color="0.75", lw=0.5,
            label="SNR", rasterized=True)
    ax.plot(freq, snr_smooth, color="C2", lw=1.8, label="SNR geglättet")
    ax.axvline(numax, color="C0", lw=1.5, ls="--",
               label=f"νmax = {numax:.1f} μHz")
    ax.set_xlabel("Frequenz (μHz)")
    ax.set_ylabel("SNR")
    ax.set_title("SNR-Spektrum")
    ax.legend(fontsize=8, loc="upper right")

    # --- Panel 4: Échelle-Diagramm (2D-Bild, geglättete PSD) ---
    ax  = axes[1, 1]
    df  = freq[1] - freq[0]
    flo = max(numax - 5.5 * deltanu, freq[0])
    fhi = min(numax + 5.5 * deltanu, freq[-1])
    mask = (freq >= flo) & (freq <= fhi)

    if np.sum(mask) > 20:
        # 1-μHz-Glättung vor der Faltung: reduziert Rauschen, erhält Modenstruktur
        sigma_sm = max(1, int(round(1.0 / df)))
        p_smooth = gaussian_filter1d(power, sigma=sigma_sm)

        x = freq[mask] % deltanu
        y = freq[mask]

        # 2D-Gitter: n_x Spalten (x = freq mod Δν), n_y Zeilen (y = freq)
        n_x   = max(30, int(deltanu / df / 4))   # ~4 Bins pro μHz
        n_y   = max(30, int((fhi - flo) / df / 4))
        xb    = np.linspace(0, deltanu, n_x + 1)
        yb    = np.linspace(flo, fhi,  n_y + 1)
        H, _, _ = np.histogram2d(x, y, bins=[xb, yb],
                                  weights=p_smooth[mask])
        C, _, _ = np.histogram2d(x, y, bins=[xb, yb])
        with np.errstate(divide="ignore", invalid="ignore"):
            grid = np.where(C > 0, H / C, np.nan)

        # Leichte 2D-Glättung für saubere Ridges
        from scipy.ndimage import gaussian_filter
        grid_s = gaussian_filter(np.nan_to_num(grid), sigma=1.0)

        im = ax.imshow(
            np.log10(np.clip(grid_s, 1e-3 * np.nanmax(grid_s), None)).T,
            origin="lower", aspect="auto",
            extent=[0, deltanu, flo, fhi],
            cmap="viridis",
        )
        plt.colorbar(im, ax=ax, label="log PSD (ppm²/μHz)")
    ax.set_xlabel(f"Frequenz mod {deltanu:.2f} μHz")
    ax.set_ylabel("Frequenz (μHz)")
    ax.set_title("Échelle-Diagramm")

    param_text = (
        f"νmax  = {numax:.1f} μHz\n"
        f"Δν    = {deltanu:.2f} μHz\n"
        f"M     = {params['mass']:.2f} M☉\n"
        f"R     = {params['radius']:.2f} R☉\n"
        f"log g = {params['logg']:.2f}\n"
        f"L     = {params['lumin']:.1f} L☉"
    )
    ax.text(0.03, 0.97, param_text,
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85))

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
    except PermissionError:
        # PDF ist noch im Viewer geöffnet (Windows-Sperre) → als PNG speichern
        outpath = outpath.with_suffix(".png")
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → Gespeichert: {outpath}")


# ===========================================================================
# Hauptprogramm
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Asteroseismologie-Workflow für solar-like Oscillators (TESS/Kepler/K2)"
    )
    p.add_argument("--tic",        default=DEFAULT_TIC,   help="Ziel-ID: TIC …, KIC …, oder EPIC …")
    p.add_argument("--name",       default=DEFAULT_NAME,  help="Sternname (für Ausgabe/Dateiname)")
    p.add_argument("--teff",       default=DEFAULT_TEFF,  type=float, help="Teff in K")
    p.add_argument("--fmin",       default=FMIN_UHZ,      type=float, help="Untere Frequenzgrenze (μHz)")
    p.add_argument("--fmax",       default=FMAX_UHZ,      type=float, help="Obere Frequenzgrenze (μHz); wird automatisch auf Nyquist gecappt")
    p.add_argument("--oversample", default=OVERSAMPLE,    type=int,   help="Oversampling-Faktor (Rechenzeit)")
    p.add_argument("--sectors",    default=1,             type=int,   help="Max. Anzahl Sektoren/Quartale (Default 1)")
    p.add_argument("--author",     default=None,          help="Datenprovider (SPOC/Kepler/K2); auto-detektiert aus TIC/KIC/EPIC-Präfix")
    p.add_argument("--exptime",    default=None,          type=int,   help="Belichtungszeit in s (120 für TESS, 1800 für Kepler LC)")
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    sep = "=" * 60
    print(f"\n{sep}")
    print(f"  Astroseismologie-Workflow: {args.name}")
    print(f"{sep}\n")

    lk.conf.cache_path = "data/cache"
    safe_name = args.name.replace(" ", "_").replace("/", "-").replace("η", "eta").replace("ε", "eps")
    outpath   = Path("results/figures") / f"asteroseismologie_{safe_name}.pdf"

    # --- Schritt 1 ---
    print("[1] Lade Lichtkurve ...")
    lc = load_lightcurve(args.tic, max_sectors=args.sectors,
                         author=args.author, exptime=args.exptime)
    span = lc.time.value[-1] - lc.time.value[0]
    print(f"  Zeitspanne: {span:.1f} d, {len(lc.flux)} Datenpunkte")

    # --- Schritt 2 ---
    print("\n[2] Berechne Powerspektrum ...")
    freq, power = compute_power_spectrum(
        lc, fmin_uHz=args.fmin, fmax_uHz=args.fmax, oversample=args.oversample
    )
    df = freq[1] - freq[0]
    dt_med = float(np.median(np.diff(lc.time.value))) * 86400.0
    nyquist = 1e6 / (2.0 * dt_med)
    print(f"  Frequenzauflösung: {df:.4f} μHz")
    print(f"  Frequenzbereich  : {freq[0]:.1f} – {freq[-1]:.1f} μHz  (Nyquist: {nyquist:.0f} μHz)")

    # --- Schritt 3 ---
    print("\n[3] Bestimme νmax (vollautomatisch, ohne Schätzwert) ...")
    numax_result = estimate_numax_auto(freq, power, verbose=True)
    numax        = numax_result["numax"]

    # --- Schritt 4 ---
    print("\n[4] Schätze Δν (Autokorrelation mit Stello-Prior) ...")
    deltanu = estimate_deltanu(freq, power, numax, verbose=True)

    # --- Schritt 5 ---
    print("\n[5] Berechne Stellar-Parameter (Skalenrelationen) ...")
    params = scaling_relations(numax, deltanu, args.teff)
    print(f"  Masse      : {params['mass']:.2f}  M☉")
    print(f"  Radius     : {params['radius']:.2f}  R☉")
    print(f"  log g      : {params['logg']:.2f}")
    print(f"  Leuchtkraft: {params['lumin']:.1f}  L☉")

    # --- Schritt 6 ---
    print("\n[6] Erstelle Abbildung ...")
    make_figure(lc, freq, power, numax_result, deltanu, params, args.name, outpath)

    print(f"\n{sep}\n")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
