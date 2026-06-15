# Astroseismologie-Workflow: Solar-like Oscillators mit TESS

Dieses Dokument beschreibt den vollständigen Python-Workflow zur Analyse
pulsierender Roter Riesen mit öffentlichen TESS-Daten — von der
Lichtkurve bis zu abgeleiteten Stellar-Parametern.

---

## Projektstruktur

```
asteroseismologie/
├── .venv/                        # uv-verwaltete virtuelle Umgebung
├── data/
│   └── cache/                    # lightkurve-Cache (automatisch befüllt)
├── results/
│   ├── figures/                  # erzeugte PDFs / PNGs
│   └── tables/                   # CSV-Ausgaben
├── asteroseismologie.py          # Hauptskript (vollständiger Workflow)
├── pyproject.toml                # Projektdefinition für uv
└── README.md                     # dieses Dokument
```

---

## Voraussetzungen

- Python ≥ 3.11
- [uv](https://docs.astral.sh/uv/) als Paketmanager
- Internetzugang (für MAST-Download der TESS-Daten)

### Empfohlene VS Code-Erweiterungen

| Erweiterung | Zweck |
|---|---|
| `ms-python.python` | Python-Interpreter, Linting |
| `ms-python.vscode-pylance` | Typprüfung, Autovervollständigung |
| `ms-toolsai.jupyter` | Interaktive Zellausführung (optional) |
| `tamasfe.even-better-toml` | `pyproject.toml`-Syntax |

---

## Einrichtung

### 1 · Projekt anlegen und Umgebung erstellen

```bash
mkdir asteroseismologie && cd asteroseismologie
uv init --no-workspace
uv add lightkurve astropy scipy matplotlib numpy
```

Das erzeugt automatisch `pyproject.toml` und `.venv/`.

### 2 · VS Code öffnen und Interpreter wählen

```bash
code .
```

Dann: **Strg+Shift+P** → *Python: Select Interpreter* →
`.venv/bin/python` (Linux/macOS) bzw. `.venv\Scripts\python.exe` (Windows).

### 3 · lightkurve-Cache konfigurieren (optional)

Damit heruntergeladene Lichtkurven nicht bei jedem Lauf neu geholt werden,
trägt man in `~/.lightkurve/config` ein:

```ini
[cache]
cache_dir = /pfad/zum/projekt/data/cache
```

Oder per Code am Skriptanfang:

```python
import lightkurve as lk
lk.conf.cache_path = "data/cache"
```

---

## `pyproject.toml`

```toml
[project]
name = "asteroseismologie"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "lightkurve>=2.4",
    "astropy>=6.0",
    "scipy>=1.12",
    "matplotlib>=3.8",
    "numpy>=1.26",
]

[tool.uv]
dev-dependencies = [
    "ipykernel",      # Jupyter-Unterstützung in VS Code
    "jupyterlab",
]
```

---

## Workflow-Übersicht

Der Ablauf besteht aus sechs Schritten, die jeweils eine klar abgegrenzte
Aufgabe haben:

```
TESS-Daten
    │
    ▼
[1] Lichtkurve laden & bereinigen
    │  lightkurve.search_lightcurve()
    │  QUALITY-Flags entfernen, Median-Normierung, Flux → ppm
    ▼
[2] Powerspektrum berechnen
    │  astropy.timeseries.LombScargle
    │  einseitige PSD in ppm²/μHz
    ▼
[3] νmax bestimmen
    │  Gaussian-Glättung der Hüllkurve
    │  Gauß-Fit im erwarteten Frequenzbereich
    ▼
[4] Δν schätzen
    │  Bandpass-Filter um νmax
    │  Autokorrelation → erstes Nebenmaximum = Δν
    ▼
[5] Skalenrelationen
    │  M, R, log g, L aus νmax, Δν, Teff
    ▼
[6] Échelle-Diagramm + 4-Panel-Figure → PDF
```

---

## Schritt 1 · Lichtkurve laden

```python
import lightkurve as lk

def load_lightcurve(tic_id: str) -> lk.LightCurve:
    """
    Lädt TESS-2-Minuten-Lichtkurven vom MAST-Archiv.

    Parameter
    ---------
    tic_id : str
        TESS Input Catalog ID, z. B. "TIC 272821450"

    Rückgabe
    --------
    lk.LightCurve
        Bereinigte Lichtkurve in ppm, alle verfügbaren Sektoren
        zusammengeführt.
    """
    search  = lk.search_lightcurve(tic_id, author="SPOC", exptime=120)
    n       = min(len(search), 4)          # maximal 4 Sektoren
    lc_coll = search[:n].download_all(quality_bitmask="hardest")

    lcs = []
    for lc in lc_coll:
        lc_clean  = lc.remove_nans().remove_outliers(sigma=4.0)
        median    = float(lc_clean.flux.value.mean())
        flux_ppm  = (lc_clean.flux.value / median - 1.0) * 1e6
        lcs.append(lk.LightCurve(time=lc_clean.time, flux=flux_ppm))

    return lk.LightCurve.from_stacked(lcs) if len(lcs) > 1 else lcs[0]
```

**Wichtige Parameter:**

| Parameter | Empfehlung | Erklärung |
|---|---|---|
| `exptime=120` | 2-Minuten-Kadenz | Rote Riesen: νmax < 300 μHz, Nyquist bei 4167 μHz — passt |
| `quality_bitmask="hardest"` | strengste Filterung | entfernt Kosmische Strahlen, Satelliten-Manöver etc. |
| `sigma=4.0` | Ausreißer-Schwelle | konservativ; 3.0 entfernt zu viele echte Pulse |

---

## Schritt 2 · Powerspektrum

```python
import numpy as np
from astropy.timeseries import LombScargle

def compute_power_spectrum(
    lc          : lk.LightCurve,
    fmin_uHz    : float = 1.0,
    fmax_uHz    : float = 300.0,
    oversample  : int   = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Lomb-Scargle PSD in μHz / ppm²/μHz.

    Frequenzauflösung: δf = 1 / (oversample · T)
    """
    time_s   = lc.time.value * 86400.0
    flux_ppm = lc.flux.value

    T_s      = time_s[-1] - time_s[0]
    df_uHz   = 1e6 / (oversample * T_s)
    freq_uHz = np.arange(fmin_uHz, fmax_uHz, df_uHz)
    freq_Hz  = freq_uHz * 1e-6

    ls    = LombScargle(time_s, flux_ppm, normalization="psd")
    power = ls.power(freq_Hz)

    # einseitige PSD, Parseval-normiert
    return freq_uHz, power * 2.0 / df_uHz
```

> **Hinweis:** Der `oversample`-Faktor verlängert die Rechenzeit spürbar.
> Für erste Tests `oversample=2` verwenden; für Publikationen `oversample=5`.

---

## Schritt 3 · νmax bestimmen

νmax ist die Frequenz, bei der die Oszillations-Hüllkurve ihr Maximum hat.
Sie hängt über eine Skalenrelation direkt mit der Oberflächengravitation
und Teff zusammen:

```
νmax / νmax,☉  ≈  (g / g☉) · (Teff / Teff,☉)^(-0.5)
```

```python
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

def fit_numax(freq, power, numax_guess):
    sigma_sm = max(1, int(0.2 * numax_guess / (freq[1] - freq[0])))
    smooth   = gaussian_filter1d(power, sigma=sigma_sm)

    mask = (freq > numax_guess / 3) & (freq < numax_guess * 3)

    def gauss(f, A, mu, sigma, bg):
        return A * np.exp(-0.5 * ((f - mu) / sigma)**2) + bg

    p0   = [smooth[mask].max(), numax_guess, numax_guess * 0.3,
            smooth[mask].min()]
    popt, _ = curve_fit(gauss, freq[mask], smooth[mask], p0=p0)

    return abs(popt[1]), smooth
```

---

## Schritt 4 · Δν via Autokorrelation

Die große Separation Δν ist der mittlere Frequenzabstand zwischen
aufeinanderfolgenden radialen Ordnungen desselben `l`-Modus.
Die Autokorrelation des (bandgefilterten) Powerspektrums verstärkt
diese periodische Struktur:

```python
def estimate_deltanu(freq, power, numax, deltanu_guess):
    width  = max(4 * deltanu_guess, numax * 0.5)
    mask   = (freq > numax - width) & (freq < numax + width)
    p_bp   = power[mask] - power[mask].mean()

    acf      = np.correlate(p_bp, p_bp, mode="full")
    acf      = acf[len(acf) // 2:]
    lag_freq = np.arange(len(acf)) * (freq[1] - freq[0])

    mask_lag = (lag_freq > 0.5 * deltanu_guess) & \
               (lag_freq < 2.0 * deltanu_guess)
    return float(lag_freq[mask_lag][np.argmax(acf[mask_lag])])
```

---

## Schritt 5 · Skalenrelationen

```python
# Sonne als Referenz
SUN_NUMAX   = 3090.0  # μHz
SUN_DELTANU = 135.1   # μHz
SUN_TEFF    = 5778.0  # K

def scaling_relations(numax, deltanu, teff):
    """Kjeldsen & Bedding (1995); Chaplin & Miglio (2013)."""
    r_nu = numax   / SUN_NUMAX
    r_dn = deltanu / SUN_DELTANU
    r_T  = teff    / SUN_TEFF

    mass   = r_nu**3  * r_dn**-4 * r_T**1.5   # M/M☉
    radius = r_nu     * r_dn**-2 * r_T**0.5   # R/R☉
    logg   = 4.44 + np.log10(r_nu) + 0.5 * np.log10(r_T)
    lumin  = radius**2 * r_T**4               # L/L☉

    return {"mass": mass, "radius": radius,
            "logg": logg, "lumin": lumin}
```

**Genauigkeit:** Skalenrelationen liefern für Rote Riesen typisch
±5 % in Radius, ±10–15 % in Masse. Für bessere Ergebnisse:
Grid-Modellierung mit BASTA oder AMP (siehe Weiterführendes).

---

## Schritt 6 · Échelle-Diagramm

Das Échelle-Diagramm entsteht durch Auftragen von `ν mod Δν` gegen `ν`.
Moden gleichen Drehimpuls-Quantums `l` erscheinen als vertikale Reihen:

```python
import matplotlib.pyplot as plt

def plot_echelle(freq, power, numax, deltanu, outpath):
    flo  = numax - 5.5 * deltanu
    fhi  = numax + 5.5 * deltanu
    mask = (freq >= flo) & (freq <= fhi)

    x = freq[mask] % deltanu   # x-Achse: Frequenz mod Δν
    y = freq[mask]             # y-Achse: Frequenz

    fig, ax = plt.subplots(figsize=(5, 7))
    sc = ax.scatter(x, y,
                    c=np.log10(power[mask]),
                    cmap="viridis", s=2, rasterized=True)
    ax.set_xlabel(f"Frequenz mod {deltanu:.1f} μHz")
    ax.set_ylabel("Frequenz (μHz)")
    plt.colorbar(sc, ax=ax, label="log PSD")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
```

**Interpretation der Spalten:**

| Spalte (x-Position) | Modus `l` | Physik |
|---|---|---|
| rechts (nahe Δν) | l = 0 | radiale Schwingung, Knoten nur im Zentrum |
| links (nahe 0) | l = 1 | dipolare Schwingung, ein Knotenkreis |
| eng neben l = 0 | l = 2 | quadrupolare Schwingung |
| Abstand l=0 zu l=2 | δν (klein) | empfindlich auf Kernzustand, Altersindikator |

---

## Empfohlene Einstiegssterne

| Stern | TIC | V | νmax (μHz) | Δν (μHz) | Bemerkung |
|---|---|---|---|---|---|
| ε Ophiuchi | TIC 272821450 | 3.2 | 59 | 5.3 | Referenz Barban+2007 |
| η Serpentis | TIC 234295610 | 3.2 | 173 | 13.4 | heller, gut beobachtet |
| KIC 4351319 | — | ~11 | 47 | 4.7 | Kepler-Archiv, viele Referenzwerte |
| KIC 11026764 | — | ~11 | 1060 | 51 | Subgiant, asymptotische Näherung gut |

Für Kepler-Sterne: `lk.search_lightcurve("KIC 4351319", author="Kepler")`.

---

## Weiterführende Schritte

### Pre-whitening für δ Scuti-Sterne

Wenn das Ziel δ Scuti-Sterne sind, kommt nach dem Powerspektrum
ein iterativer Pre-whitening-Loop:

```python
def prewhiten(freq, power, n_modes=30):
    """Iterative Frequenzextraktion."""
    residual = power.copy()
    modes    = []
    for _ in range(n_modes):
        i_max = np.argmax(residual)
        f0    = freq[i_max]
        A0    = np.sqrt(residual[i_max])
        # Subtrahiere Sinusbeitrag vom Residuum
        sigma = 3 * (freq[1] - freq[0])
        residual -= A0**2 * np.exp(-0.5 * ((freq - f0) / sigma)**2)
        residual  = np.clip(residual, 0, None)
        modes.append((f0, A0))
    return modes
```

### Grid-basierte Modellierung

Für publizierbare Stellar-Parameter:
- **BASTA** (BAyesian STellar Algorithm): `pip install basta`
- **AMP** (Asteroseismic Modeling Portal): webbasiert, mpi.arizona.edu/amp
- Eigene MESA-Gitter: erfordert MESA-Installation, aber frei skalierbar

### Datenquellen

| Archiv | URL | Inhalt |
|---|---|---|
| MAST | mast.stsci.edu | TESS, Kepler, K2 — über lightkurve direkt zugänglich |
| ASAS-SN | asas-sn.osu.edu | Bodenbasierte V-Band-Photometrie, gut für Cepheiden |
| Gaia DR3 | gea.esac.esa.int | Variabilitätskatalog mit vorklassifizierten Typen |
| KASOC | kasoc.phys.au.dk | Aufbereitete Kepler-Lichtkurven für Astroseismologie |

### Literatur (Open Access)

- Aerts, Christensen-Dalsgaard & Kurtz (2010): *Asteroseismology* — Standardwerk
- Chaplin & Miglio (2013, ARA&A 51, 353): Solar-like oscillations, Review
- Mosser & Appourchaux (2009, A&A 508, 877): Autokorrelationsmethode für Δν
- Bowman (2021, Front. Astron. Space Sci.): δ Scuti-Sterne, aktueller Review

---

## Komplettes Ausführungsbeispiel

```bash
# Umgebung aktivieren (uv)
source .venv/bin/activate          # Linux/macOS
.venv\Scripts\activate             # Windows

# Einzellauf
python asteroseismologie.py

# Oder mit uv direkt
uv run asteroseismologie.py
```

Erwartete Konsolenausgabe für ε Ophiuchi:

```
============================================================
  Astroseismologie-Workflow: ε Ophiuchi
============================================================

[1] Lade TESS-Lichtkurve ...
  Gefunden: 3 Sektor(en)
  Lade 3 Sektor(en) ...
  Zeitspanne: 81.3 Tage, 56891 Datenpunkte

[2] Berechne Powerspektrum ...
  Frequenzauflösung: 0.0028 μHz
  Frequenzbereich  : 1.0 – 300.0 μHz

[3] Bestimme νmax ...
  νmax = 58.4 μHz  (Literatur: 59.0 μHz)

[4] Schätze Δν (Autokorrelation) ...
  Δν    = 5.31 μHz  (Literatur: 5.3 μHz)

[5] Berechne Stellar-Parameter (Skalenrelationen) ...
  Masse    : 1.87  M☉   (Lit. ~1.85 M☉)
  Radius   : 10.4  R☉   (Lit. ~10.6 R☉)
  log g    : 2.75       (Lit. ~2.73)
  Leuchtkraft: 52    L☉   (Lit. ~54 L☉)

[6] Erstelle Abbildung ...
  → Gespeichert: results/figures/asteroseismologie_eps_oph.pdf
```
