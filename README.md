# Asteroseismologie-Workflow

Python-Pipeline zur vollautomatischen Analyse solar-like Oscillators
(Rote Riesen, Unterriesen, Sonnenanaloga) mit öffentlichen Raumdaten
von TESS und Kepler/K2 — von der Lichtkurve bis zu abgeleiteten
Stellar-Parametern.

## Funktionsumfang

- **Multi-Mission**: TESS SC (2 min), Kepler LC (30 min), Kepler SC (1 min),
  K2 LC (30 min) — automatisch erkannt aus TIC/KIC/EPIC-Präfix
- **Lichtkurvenvorverarbeitung**: Qualitätsfilter, Savitzky-Golay-Detrending
  (7-Tage-Fenster, entfernt Rotation und Instrumentendrifts), RMS-Filterung
- **Lomb-Scargle-PSD** in ppm²/μHz (Parseval-normiert); fmax wird bei
  Kurzkadenz automatisch auf 0.9 × Nyquist angehoben
- **νmax vollautomatisch** (ohne jeden Schätzwert):
  1. Harvey-Hintergrundfit im Log-Raum (adaptiver Fit-Bereich LC/SC)
  2. SNR-Spektrum (PSD / Hintergrundmodell)
  3. Geglättetes SNR → erster Kandidat
  4. Gauß-Fit → finales νmax
- **Δν** via Autokorrelation (Stello-Prior ±20 %) + Échelle-Kohärenz-Verfeinerung
- **Skalenrelationen** M, R, log g, L (Chaplin & Miglio 2013)
- **4-Panel-PDF**: Lichtkurve · PSD + Harvey · SNR-Spektrum · Échelle-Diagramm

## Voraussetzungen

- Python ≥ 3.11
- [uv](https://docs.astral.sh/uv/) als Paketmanager
- Internetzugang (MAST-Download; beim ersten Lauf je Stern)

## Einrichtung

```bash
uv sync
```

## Verwendung

```bash
# Standardstern: η Serpentis (TIC 234295610, Roter Riese, TESS SC)
uv run asteroseismologie.py

# TESS-Riese: ε Ophiuchi, mehrere Sektoren für besseres SNR
uv run asteroseismologie.py --tic "TIC 272821450" --name "eps Oph" --teff 4700 --sectors 4

# Kepler Roter Riese (Langkadenz, Nyquist 278 μHz)
uv run asteroseismologie.py --tic "KIC 4351319" --name "KIC 4351319" --teff 4800

# Kepler Sonnenanaloge (Kurzkadenz, Nyquist 8333 μHz)
uv run asteroseismologie.py --tic "KIC 10963065" --name "KIC 10963065" --teff 5900 --exptime 60 --sectors 8

# Schneller Test (halbe Frequenzauflösung)
uv run asteroseismologie.py --oversample 2

# Alle Optionen
uv run asteroseismologie.py --help
```

## Empfohlene Einstiegssterne

| Stern | ID | Teff (K) | νmax (μHz) | Δν (μHz) | Mission |
|-------|----|---------|-----------|---------|--------|
| η Serpentis | TIC 234295610 | 4970 | ~173 | ~13.4 | TESS SC |
| ε Ophiuchi | TIC 272821450 | 4700 | ~59 | ~5.3 | TESS SC |
| KIC 4351319 | KIC 4351319 | 4800 | ~47 | ~4.7 | Kepler LC |
| KIC 10963065 | KIC 10963065 | 5900 | ~2000 | ~103 | Kepler SC (`--exptime 60`) |

## Projektstruktur

```
asteroseismologie/
├── asteroseismologie.py            # vollständige Pipeline
├── asteroseismologie_workflow.md   # ausführliche Dokumentation
├── numax_ohne_schaetzwert.md       # νmax-Pipeline im Detail
├── pyproject.toml
├── data/cache/                     # lightkurve-Cache (automatisch)
└── results/figures/                # erzeugte PDFs
```

## Literatur

- Chaplin & Miglio (2013, ARA&A 51, 353) — Solar-like oscillations, Review
- Stello et al. (2009, ApJ 700, 1589) — νmax–Δν Skalenrelation
- Mosser & Appourchaux (2009, A&A 508, 877) — Autokorrelationsmethode
- Harvey (1985, ESA SP-235) — Granulations-Hintergrundmodell


## Funktionsumfang

- TESS-Lichtkurven automatisch vom MAST-Archiv laden (lightkurve/SPOC)
- Lomb-Scargle-Powerspektrum (ppm²/μHz, Parseval-normiert)
- **νmax vollautomatisch**, ohne Schätzwert:
  1. Harvey-Hintergrundfit im Log-Raum
  2. SNR-Spektrum (PSD / Hintergrund)
  3. Breitbandig geglättetes SNR → erster Kandidat
  4. Gauß-Fit → finales νmax
- Δν via Autokorrelation (Suchbereich automatisch aus νmax abgeleitet)
- Stellar-Parameter M, R, log g, L aus Skalenrelationen
- 4-Panel-PDF: Lichtkurve, PSD, SNR-Spektrum, Échelle-Diagramm

## Voraussetzungen

- Python ≥ 3.11
- [uv](https://docs.astral.sh/uv/) als Paketmanager
- Internetzugang (MAST-Download)

## Einrichtung

```bash
uv sync
```

## Verwendung

```bash
# Standardstern: η Serpentis (TIC 234295610)
uv run asteroseismologie.py

# Anderer Stern
uv run asteroseismologie.py --tic "TIC 272821450" --name "eps Oph" --teff 4900

# Schneller Test (geringerer Oversample-Faktor)
uv run asteroseismologie.py --oversample 2

# Hilfe
uv run asteroseismologie.py --help
```

Ergebnis-PDFs werden unter `results/figures/` gespeichert.

## Empfohlene Einstiegssterne

| Stern        | TIC          | νmax (μHz) | Δν (μHz) |
|--------------|--------------|-----------|---------|
| η Serpentis  | TIC 234295610 | 173       | 13.4    |
| ε Ophiuchi   | TIC 272821450 | 59        | 5.3     |

## Projektstruktur

```
asteroseismologie/
├── asteroseismologie.py   # vollständige Pipeline
├── pyproject.toml
├── data/cache/            # lightkurve-Cache (automatisch befüllt, nicht eingecheckt)
└── results/figures/       # erzeugte PDFs (nicht eingecheckt)
```

## Literatur

- Chaplin & Miglio (2013, ARA&A 51, 353) — Solar-like oscillations, Review
- Stello et al. (2009, ApJ 700, 1589) — νmax–Δν Skalenrelation
- Mosser & Appourchaux (2009, A&A 508, 877) — Autokorrelationsmethode
