# Asteroseismologie-Workflow

Python-Pipeline zur vollautomatischen Analyse solar-like Oscillators
(Rote Riesen, Unterriesen, Sonnenanaloga) mit Daten aus TESS, Kepler und K2:
von der Lichtkurve bis zu abgeleiteten Stellar-Parametern.

## Funktionsumfang

- Multi-Mission-Unterstützung mit Auto-Erkennung aus Zielpräfix:
  TIC -> TESS (SPOC), KIC -> Kepler, EPIC -> K2
- Lichtkurvenaufbereitung: Quality-Filter, Outlier-Filter, 7-Tage-Detrending,
  Segmentselektion über Zeitspanne und RMS
- Lomb-Scargle-PSD in ppm²/μHz (einseitig, Parseval-normiert)
- Vollautomatische νmax-Bestimmung ohne Startschätzwert:
  1. Harvey-Hintergrundfit
  2. SNR-Spektrum (PSD / Hintergrund)
  3. Pre-Whitening + adaptive Glättung
  4. Gauß-Verfeinerung des Peaks
- Δν-Schätzung aus ACF (Stello-Prior) plus Échelle-Kohärenz/Ridge-Refinement
- Skalenrelationen für M, R, log g, L
- Ausgabe als 4-Panel-Figur (Lichtkurve, PSD+Harvey, SNR, Échelle)
- Zusätzliches repliziertes Échelle mit farbcodierten l=0/1/2-Modenkandidaten
- Optionale Bayesianische Modenidentifikation und Peakbagging mit PBjam >= 2.0
- Evidenzgestützte PBjam-Auswahl mit getrennten Gold- und Silber-Moden

## Voraussetzungen

- Python >= 3.11
- uv als Paketmanager: https://docs.astral.sh/uv/
- Internetzugang für den MAST-Download (beim ersten Lauf je Ziel)

## Installation

```bash
uv sync
```

Für die optionale PBjam-Integration:

```bash
uv sync --extra pbjam
```

## Verwendung

```bash
# Standardziel: eta Serpentis (TIC 234295610)
uv run asteroseismologie.py

# TESS-Roter-Riese
uv run asteroseismologie.py --tic "TIC 272821450" --name "eps Oph" --teff 4700 --sectors 4

# Kepler-Roter-Riese (LC)
uv run asteroseismologie.py --tic "KIC 4351319" --name "KIC 4351319" --teff 4800

# Kepler-Sonnenanaloge (SC)
uv run asteroseismologie.py --tic "KIC 10963065" --name "KIC 10963065" --teff 5900 --exptime 60 --sectors 8

# Alternativer SNR-Glätter (Gauß via FFT)
uv run asteroseismologie.py --gauss-smooth

# PBjam ModeID + Peakbagging mit Gaia-Farbe
uv run asteroseismologie.py --tic "KIC 10273246" --name "Mulder" --teff 6150 --exptime 60 --sectors 8 --pbjam --bp-rp 0.60

# Schneller Test
uv run asteroseismologie.py --oversample 2

# Alle Optionen
uv run asteroseismologie.py --help
```

Die Einzelstern-CLI bleibt unabhängig von der Atlas-Konfiguration vollständig
nutzbar. Für einen seriellen, wiederaufnehmbaren Atlaslauf:

```bash
# Bereits analysierte Targets aus config/targets.yaml
uv run python -m atlas.run_all

# Einzelnes ausstehendes Target erstmals analysieren
uv run python -m atlas.run_all --target "KIC 4351319"

# Bereits geprüften PBjam-Alt-Cache bei einer Migration explizit übernehmen
uv run python -m atlas.run_all --target "KIC 10963065" --reuse-existing-pbjam --no-legacy-figures

# Atlasplots ausschließlich aus vorhandenen Ergebnisartefakten erzeugen
uv run python -m atlas.generate_plots
```

Atlas-Ergebnisse liegen pro Stern unter `results/atlas/<target>/`. Ein
`summary.json` enthält Messwerte, Unsicherheiten, Qualitätsflags und Provenienz;
`diagnostics.npz` hält die Arrays für reproduzierbare Folgeauswertungen. Manuell
geprüfte Δν-Werte sind nur als begründeter `deltanu_override` in der
Target-Konfiguration zulässig. KIC 4351319 ist derzeit als noch nicht untersucht
(`pending`) markiert.

Wichtige CLI-Optionen:

- --tic: Ziel-ID (TIC..., KIC..., EPIC...)
- --name: Name für Ausgabe/Dateinamen
- --teff: Effektive Temperatur in K
- --sectors: Maximalzahl geladener Segmente
- --author: Datenquelle überschreiben (SPOC, Kepler, K2)
- --exptime: Belichtungszeit in s (z. B. 120 oder 1800, Kepler SC: 60)
- --fmin, --fmax: Frequenzbereich in μHz
- --oversample: Frequenzauflösung vs. Laufzeit
- --gauss-smooth: Gauß-Glättung statt Boxcar
- --echelle-replicas: Zwei oder drei horizontale Échelle-Kopien (Standard: 2)
- --pbjam: PBjam ModeID und Peakbagging aktivieren (sonst unveränderter schneller Modus)
- --pbjam-numax-sigma, --pbjam-deltanu-sigma, --teff-sigma: Unsicherheiten der PBjam-Priors
- --bp-rp, --bp-rp-sigma: Optionale Gaia-Farbe mit Unsicherheit
- --pbjam-orders: Anzahl modellierter radialer Ordnungen (Standard: 7)
- --pbjam-quality-min: Optionaler fester Quality-Cut; standardmäßig adaptiv (Median × 1.10, begrenzt auf 0.75 bis 2.0)
- --pbjam-fap-gold, --pbjam-fap-silver: Empirisch kalibrierte Schwellen des Höhen-FAP-Proxys
- --pbjam-ridge-tol: Maximale Abweichung schwacher l=0/2-Moden vom Anker-Ridge in μHz
- --pbjam-d02-fraction: Optionaler fester δν₀₂/Δν-Wert; sonst Gold-Paare oder Δν-Relation
- --pbjam-sequence-tolerance, --pbjam-sequence-minimum: Regularitätstest konsekutiver l=0/2-Moden
- --pbjam-refresh: Vorhandenen Ergebnis-Cache ignorieren

## Konfidenz der PBjam-Moden

PBjams `height` korreliert mit der Detektionsstärke, ist aber keine extern
kalibrierte SNR. Deshalb wird `exp(-height)` als **FAP-Proxy** und nicht als
kalibrierte Falschalarmwahrscheinlichkeit ausgewiesen. Die Standardgrenzen
`0.01` und `0.1` entsprechen den empirischen Höhenschnitten 4.61 und 2.30.
Die Summary speichert zusätzlich die nominale globale Rate
`fap_global_gold = 1 - (1 - fap_gold) ** n_modes_tested`. Diese Šidák-Größe
setzt unabhängige Tests und die ideale Höhenstatistik voraus; sie beschreibt
nicht die zusätzliche Quality-, Ridge- oder Sequenzevidenz.

Gold-Moden dürfen für quantitative Messungen verwendet werden. Silber-Moden
werden im Échelle offen dargestellt. Der l=0-Ridge wird ausschließlich aus
Gold-Ankern bestimmt. Der l=2-Versatz wird zuerst aus vorhandenen Gold-Paaren
geschätzt; ohne solche Paare gilt
`δν₀₂/Δν = -0.0324 log10(Δν) + 0.1388`. Bereits anderweitig bestätigte
l=0/2-Moden werden zu Gold hochgestuft, wenn mindestens drei konsekutive
Ordnungen Abstände innerhalb `Δν ± 10 %` bilden. Der Sequenztest wird niemals
auf l=1 angewandt, weil gemischte Moden und Avoided Crossings dessen
Regularitätsannahme verletzen.

## Ausgabe

- Cache: data/cache/
- Figuren: results/figures/
- PBjam-Ergebnisse und Cache: results/pbjam/
- PBjam-Modentabellen: results/tables/

Standardmäßig werden die 4-Panel-Übersicht und ein repliziertes Échelle als
separate PDF-Dateien gespeichert. Das replizierte Diagramm markiert automatisch
extrahierte Kandidaten für l=0 (blau), l=1 (grün), l=2 (orange) und unklare bzw.
gemischte Moden (grau). Die Zuordnung aus der Échelle-Position ist diagnostisch
und ersetzt kein wissenschaftliches Peak-Bagging. Falls eine PDF unter Windows
noch geöffnet ist (Dateisperre), schreibt das Skript automatisch eine PNG-Datei.
Mit `--pbjam` werden stattdessen die Bayesianischen PBjam-Labels und
Peakbagging-Höhen verwendet. Identische Eingaben laden automatisch den JSON-Cache;
`--pbjam-refresh` erzwingt einen neuen, potenziell lang laufenden Sampling-Lauf.

## Weiterführende Dokumentation

- Gesamtworkflow und methodische Details: [asteroseismologie_workflow.md](asteroseismologie_workflow.md)
- Detaillierte Beschreibung der νmax-Bestimmung ohne Schätzwert: [numax_ohne_schaetzwert.md](numax_ohne_schaetzwert.md)

## Empfohlene Einstiegsziele

| Stern | ID | Teff (K) | νmax (μHz) | Δν (μHz) | Mission |
|---|---|---:|---:|---:|---|
| eta Serpentis | TIC 234295610 | 4970 | ~173 | ~13.4 | TESS SC |
| eps Ophiuchi | TIC 272821450 | 4700 | ~59 | ~5.3 | TESS SC |
| KIC 4351319 | KIC 4351319 | 4800 | ~47 | ~4.7 | Kepler LC |
| KIC 10963065 | KIC 10963065 | 5900 | ~2000 | ~103 | Kepler SC |

## Projektstruktur

```text
asteroseismology/
|- asteroseismologie.py
|- pbjam_bridge.py
|- asteroseismologie_workflow.md
|- numax_ohne_schaetzwert.md
|- pyproject.toml
|- data/
|  |- cache/
|- results/
|  |- figures/
|  |- tables/
```

## Abhängigkeiten

- lightkurve
- astropy
- scipy
- matplotlib
- numpy

## Literatur

- Chaplin & Miglio (2013, ARA&A 51, 353)
- Stello et al. (2009, ApJ 700, 1589)
- Mosser & Appourchaux (2009, A&A 508, 877)
- Harvey (1985, ESA SP-235)
