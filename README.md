# Asteroseismologie-Workflow

Python-Pipeline zur Analyse solar-like Oscillators (pulsierender Roter Riesen)
mit öffentlichen TESS-Daten — von der Lichtkurve bis zu abgeleiteten Stellar-Parametern.

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
