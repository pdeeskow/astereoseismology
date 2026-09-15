# Konzept: PBjam-Integration in die Astroseismologie-Pipeline

Erweiterung der bestehenden Pipeline um Bayesianische Modenidentifikation
und Peakbagging — mit dem Ziel, gemischte l=1-Moden und avoided crossings
bei Unterriesen (KIC 10273246 "Mulder") korrekt zu identifizieren.

---

## 1  Was PBjam leistet

PBjam ist ein Bayesianisches Toolkit für automatisierte Modenidentifikation
und Peakbagging von solar-like oscillators (Nielsen et al. 2021, AJ 161, 62).

### 1.1  Die zwei Kernaufgaben

| Schritt | Aufgabe | Warum das schwierig ist |
|---|---|---|
| **Modenidentifikation** | Welcher Peak hat welchen Grad l? | Bei gemischten Moden sitzen l=1-Peaks nicht an der erwarteten asymptotischen Position |
| **Peakbagging** | Präzise Frequenz, Höhe, Breite jeder Mode | Erfordert Modellfit ans Spektrum mit korrekter l-Zuordnung als Vorwissen |

### 1.2  Der methodische Unterschied zur eigenen Pipeline

Die bisherige Pipeline ordnet Moden **nach ihrer Position** im Échelle zu
(Ansatz A aus dem Échelle-Konzept). Das versagt systematisch bei gemischten
Moden, weil diese definitionsgemäß von der erwarteten Position abweichen —
sie landen in der Kategorie "unklar".

PBjam geht anders vor:

1. **l=2,0-Paare zuerst**: Die asymptotische Relation für p-Moden wird an
   die l=2,0-Paare gefittet. Diese sind stabil und nahezu unbeeinflusst von
   Kopplung — sie bilden den Anker.

2. **l=1 mit Kopplungsmodell**: Anschließend wird eine Auswahl von Modellen
   für die l=1-Moden angewandt, wobei jedes Modell für ein anderes
   Entwicklungsstadium geeignet ist. Für Unterriesen im Kopplungsregime
   beschreibt das Modell die gemischten Moden explizit — inklusive der
   avoided crossings.

3. **Prior aus 13.000+ Sternen**: Der Prozess stützt sich auf einen großen
   Satz früherer Beobachtungen (Kepler, K2, TESS), aus denen eine
   Prior-Verteilung konstruiert wird. Das macht die Identifikation robust,
   auch bei mittelmäßigem SNR.

4. **Detailliertes Peakbagging**: Nach der Identifikation folgt ein Modell
   mit weniger Constraints, das die Frequenzen um einige μHz anpassen darf —
   für akustische Glitches, Rotation und Effekte, die das
   Identifikationsmodell nicht erfasst.

### 1.3  Version 2 — entscheidend für Unterriesen

> Die initiale Version von PBjam war darauf beschränkt, nur Moden mit
> ℓ = 0 und ℓ = 2 zu identifizieren, da diese am einfachsten konsistent
> über verschiedene Entwicklungsstadien zu behandeln sind.

PBjam **2.0** ergänzt die Messung von Dipol-Moden-Frequenzen in
Kopplungsregimen von der Hauptreihe bis zu niedrigleuchtkräftigen Roten
Riesen. Genau dieser Bereich umfasst Mulder (log g 3.90, entwickelter
Unterriese) und KIC 10963065 (log g 4.29, früher Unterriese).

**Konsequenz:** Für dieses Projekt ist zwingend PBjam ≥ 2.0 erforderlich.

---

## 2  Architektur — drei Nutzungsebenen

PBjam bietet drei Abstraktionsebenen. Die Wahl hängt davon ab, wie viel
Kontrolle die eigene Pipeline behalten soll.

```
┌──────────────────────────────────────────────────────┐
│  session   — höchste Ebene                           │
│  · lädt Daten selbst (via lightkurve)                │
│  · berechnet PSD selbst                              │
│  · macht ModeID + Peakbagging + Output               │
│  → für Stapelverarbeitung vieler Sterne              │
├──────────────────────────────────────────────────────┤
│  star      — mittlere Ebene                          │
│  · nimmt eigene Zeitreihe oder eigenes Spektrum      │
│  · macht ModeID + Peakbagging                        │
│  → EMPFOHLEN für Integration in eigene Pipeline      │
├──────────────────────────────────────────────────────┤
│  modeID / peakbag — niedrigste Ebene                 │
│  · Einzelschritte separat ansteuerbar                │
│  → für Debugging und methodische Experimente         │
└──────────────────────────────────────────────────────┘
```

Die `session`-Klasse ist im Grunde ein Wrapper um die `star`-Klasse.
Für die Integration in eine bestehende Pipeline ist **`star`** die richtige
Wahl: Die eigene Pipeline behält die Kontrolle über Datenbeschaffung,
Bereinigung und PSD-Berechnung; PBjam übernimmt nur die Modenidentifikation.

---

## 3  Integrationskonzept

### 3.1  Schnittstelle zur bestehenden Pipeline

Die bestehende Pipeline liefert bereits alles, was PBjam als Eingabe
braucht:

| Bestehende Pipeline liefert | PBjam braucht |
|---|---|
| `freq`, `power` (PSD aus Lomb-Scargle) | Frequenzachse + PSD |
| `numax` aus Harvey-/Gauß-Fit | νmax mit Unsicherheit |
| `deltanu` aus Autokorrelation | Δν mit Unsicherheit |
| Teff (aus Literatur/Gaia) | Teff mit Unsicherheit |
| — | Gaia-Farbe (bp_rp), optional |

**Wichtig:** PBjam erwartet νmax und Δν jeweils als `[Wert, Unsicherheit]`.
Die Unsicherheiten fließen in den Prior ein — zu enge Angaben können die
Identifikation verfälschen, zu weite verlangsamen die Konvergenz.

### 3.2  Datenfluss

```
   Bestehende Pipeline                    PBjam
   ───────────────────                    ─────
   TESS/Kepler laden
        ↓
   Lichtkurve bereinigen
        ↓
   Lomb-Scargle → PSD
        ↓
   Harvey-Fit → SNR
        ↓
   νmax bestimmen ──────────────────────►  Prior-Input
        ↓
   Δν (ACF + Stello-Prior) ─────────────►  Prior-Input
        ↓                                       ↓
   [freq, power] ───────────────────────►  star(...)
                                               ↓
                                          ModeID (l=0,2 → l=1)
                                               ↓
                                          Peakbagging
                                               ↓
   ◄──────────────────────────────────  {ν, l, n, Höhe, Breite, σ}
        ↓
   Repliziertes Échelle mit
   echten l-Labels
        ↓
   avoided crossings detektieren
```

### 3.3  Was zurückkommt

PBjam liefert für jede Mode:

- **Frequenz** ν mit Unsicherheit
- **Grad l** (0, 1, 2 — bei v2 auch gemischte l=1)
- **radiale Ordnung n**
- **Modenhöhe** und **-breite** (Lorentzprofil-Parameter)
- **Prior/Posterior-Breitenverhältnis** als Qualitätsmetrik

Zur Qualitätsmetrik: Moden mit einem Verhältnis von Prior- zu
Posterior-Breite größer als etwa 2 können verwendet werden, um Moden
auszuwählen, bei denen das Spektrum die Information dominiert (statt des
Priors). Das ist der entscheidende Filter für vertrauenswürdige Moden.

---

## 4  Implementierungsvorschlag

### 4.1  Neue Modulstruktur

```
asteroseismologie/
├── pipeline/
│   ├── data.py           # bestehend: Lichtkurve laden
│   ├── spectrum.py       # bestehend: PSD, Harvey, νmax, Δν
│   ├── echelle.py        # bestehend: repliziertes Échelle
│   ├── pbjam_bridge.py   # NEU: Schnittstelle zu PBjam
│   └── mixed_modes.py    # NEU: avoided-crossing-Analyse
├── results/
│   └── pbjam/            # PBjam-Ausgabeverzeichnis
└── main.py
```

### 4.2  Die Brücke — `pbjam_bridge.py`

```python
"""
Schnittstelle zwischen der eigenen Pipeline und PBjam.

Übergibt die bereits berechnete PSD und die seismischen Globalparameter
an PBjam, holt die Modenidentifikation zurück.
"""

import numpy as np
from pathlib import Path


def run_pbjam_modeid(
    freq_uHz    : np.ndarray,
    power_psd   : np.ndarray,
    numax       : tuple[float, float],   # (Wert, Unsicherheit) in μHz
    deltanu     : tuple[float, float],   # (Wert, Unsicherheit) in μHz
    teff        : tuple[float, float],   # (Wert, Unsicherheit) in K
    bp_rp       : tuple[float, float] | None = None,
    star_id     : str = "star",
    outdir      : str = "results/pbjam",
):
    """
    Führt PBjams Modenidentifikation und Peakbagging auf einem bereits
    berechneten Spektrum aus.

    Parameter
    ---------
    freq_uHz  : Frequenzachse der PSD in μHz
    power_psd : Leistungsdichte (ppm²/μHz)
    numax     : νmax mit 1σ-Unsicherheit
    deltanu   : Δν mit 1σ-Unsicherheit
    teff      : Effektivtemperatur mit Unsicherheit (aus Gaia/Literatur)
    bp_rp     : Gaia-Farbe BP−RP mit Unsicherheit (optional, verbessert Prior)
    star_id   : Bezeichner für Ausgabedateien
    outdir    : Ausgabeverzeichnis

    Rückgabe
    --------
    dict mit den identifizierten Moden:
        'freq'   : Frequenzen (μHz)
        'freq_err': Unsicherheiten
        'l'      : Grad
        'n'      : radiale Ordnung
        'height' : Modenhöhe
        'width'  : Modenbreite
        'quality': Prior/Posterior-Breitenverhältnis

    Hinweis
    -------
    Die genauen Argumentnamen der star-Klasse können sich zwischen
    PBjam-Versionen unterscheiden. Vor der Umsetzung die aktuellen
    Beispiel-Notebooks prüfen:
    https://pbjam.readthedocs.io/en/latest/examples.html
    """
    import pbjam

    Path(outdir).mkdir(parents=True, exist_ok=True)

    # PBjam erwartet die PSD als (frequenz, power)-Tupel
    # Die star-Klasse akzeptiert entweder eine Zeitreihe oder ein Spektrum
    st = pbjam.star(
        ID       = star_id,
        pg       = (freq_uHz, power_psd),   # bereits berechnete PSD
        numax    = list(numax),
        dnu      = list(deltanu),
        teff     = list(teff),
        bp_rp    = list(bp_rp) if bp_rp else None,
        path     = outdir,
    )

    # Vollständiger Durchlauf: ModeID + Peakbagging
    st()

    return _extract_modes(st)


def _extract_modes(st) -> dict:
    """
    Extrahiert die Modenparameter aus dem PBjam-star-Objekt in ein
    einfaches Dictionary für die weitere Verarbeitung.

    Die Struktur der PBjam-Ausgabe hängt von der Version ab — dieser
    Adapter kapselt diese Abhängigkeit an einer Stelle.
    """
    # PBjam legt die Ergebnisse als Summary-DataFrame ab
    summary = st.peakbag.summary   # exakter Attributname versionsabhängig

    modes = {
        "freq"     : summary["nu_mean"].values,
        "freq_err" : summary["nu_std"].values,
        "l"        : summary["l"].values.astype(int),
        "n"        : summary["n"].values.astype(int),
        "height"   : summary.get("height", np.full(len(summary), np.nan)),
        "width"    : summary.get("width",  np.full(len(summary), np.nan)),
    }

    # Qualitätsmetrik: Prior/Posterior-Breitenverhältnis
    if "prior_width" in summary and "nu_std" in summary:
        modes["quality"] = summary["prior_width"].values / summary["nu_std"].values
    else:
        modes["quality"] = np.full(len(summary), np.nan)

    return modes


def adaptive_quality_threshold(modes: dict, factor: float = 1.10,
                               floor: float = 0.75,
                               ceiling: float = 2.0) -> float:
    """Leitet den Quality-Cut aus dem Median der endlichen Werte ab."""
    quality = np.asarray(modes["quality"], dtype=float)
    quality = quality[np.isfinite(quality)]
    if not len(quality):
        return floor
    return float(np.clip(np.median(quality) * factor, floor, ceiling))


def filter_reliable_modes(modes: dict, quality_min: float | None = None,
                          height_min: float = 1.0) -> dict:
    """Filtert gemeinsam nach informativem Posterior und messbarer Höhe."""
    if quality_min is None:
        quality_min = adaptive_quality_threshold(modes)
    mask = ((modes["quality"] >= quality_min)
            & (modes["height"] >= height_min))
    return {k: (v[mask] if isinstance(v, np.ndarray) else v)
            for k, v in modes.items()}
```

### 4.3  Anschluss ans replizierte Échelle

Die bestehende Échelle-Funktion braucht nur die l-Labels aus PBjam statt
der eigenen Positionszuordnung:

```python
# ALT — eigene Positionszuordnung:
# labels = classify_l_by_pairs(freq_modes, amp_modes, dnu)

# NEU — PBjam-Identifikation:
modes  = run_pbjam_modeid(freq, power_psd,
                          numax=(856.2, 25.0),
                          deltanu=(48.60, 0.15),
                          teff=(6150, 100),
                          star_id="KIC10273246")
modes  = filter_reliable_modes(modes, quality_min=2.0)

plot_replicated_echelle(
    freq_modes = modes["freq"],
    amp_modes  = modes["height"],
    labels     = modes["l"],       # echte Bayesianische l-Zuordnung
    dnu        = 48.60,
    numax      = 856.2,
    n_replicas = 2,
)
```

---

## 5  Avoided-crossing-Analyse mit echten l-Labels

Mit korrekt identifizierten l=1-Moden wird die Detektion der avoided
crossings deutlich zuverlässiger als über die Positionsabweichung.

### 5.1  Prinzip

Bei reinen p-Moden ist der Frequenzabstand zwischen aufeinanderfolgenden
l=1-Moden konstant ≈ Δν. An einem avoided crossing wird dieser Abstand
**gestört** — eine zusätzliche gemischte Mode schiebt sich dazwischen, und
die lokale Abweichung von Δν springt.

```python
def find_avoided_crossings(freq_l1: np.ndarray, dnu: float,
                           threshold: float = 0.25) -> dict:
    """
    Detektiert avoided crossings über gestörte Frequenzabstände im
    l=1-Ridge.

    Parameter
    ---------
    freq_l1   : Frequenzen der l=1-Moden (μHz), aufsteigend sortiert
    dnu       : große Separation (μHz)
    threshold : relative Abweichung von Δν, ab der ein Crossing gilt

    Rückgabe
    --------
    dict mit 'crossings' (Frequenzen), 'spacings', 'deviations'
    """
    f = np.sort(freq_l1)
    spacings   = np.diff(f)
    deviations = np.abs(spacings - dnu) / dnu

    idx = np.where(deviations > threshold)[0]
    # Crossing liegt zwischen den beiden Moden
    crossings = (f[idx] + f[idx + 1]) / 2

    return {
        "crossings"  : crossings,
        "spacings"   : spacings,
        "deviations" : deviations,
        "freq_l1"    : f,
    }
```

### 5.2  Erwartungswerte zur Validierung

| Stern | erwartete avoided crossings | Quelle |
|---|---|---|
| KIC 10273246 (Mulder) | 2 | Campante et al. 2011 |
| KIC 10920273 (Scully) | 1 (+1 wahrscheinlich) | Campante et al. 2011 |
| KIC 10963065 | 0 bis wenige (kompakter, jünger) | — |

Die Reproduktion dieser Zahlen ist der Test, ob die PBjam-Integration
korrekt arbeitet.

---

## 6  Praktische Hinweise

### 6.1  Installation

```bash
uv add pbjam
# oder
pip install pbjam
```

PBjam bringt schwergewichtige Abhängigkeiten mit (PyMC / Sampling-Backend).
Bei Installationsproblemen die aktuelle Anleitung prüfen:
https://pbjam.readthedocs.io/en/latest/setup.html

### 6.2  Rechenzeit

PBjam führt MCMC-Sampling durch — das dauert deutlich länger als die
bisherige Pipeline. Grobe Größenordnung: Minuten bis Zehnminuten pro Stern,
abhängig von Spektrumsgröße und Modenzahl. Für Stapelverarbeitung ist die
`session`-Klasse mit Parallelisierung vorgesehen.

**Empfehlung:** Ergebnisse cachen. Die PBjam-Ausgabe pro Stern einmal
berechnen und als Datei ablegen, statt bei jedem Plot-Durchlauf neu zu
sampeln.

```python
import json
from pathlib import Path

def cached_pbjam(star_id, cache_dir="results/pbjam/cache", **kwargs):
    """Führt PBjam nur aus, wenn kein Cache existiert."""
    cache = Path(cache_dir) / f"{star_id}_modes.json"
    if cache.exists():
        with open(cache) as fh:
            data = json.load(fh)
        return {k: np.array(v) for k, v in data.items()}

    modes = run_pbjam_modeid(star_id=star_id, **kwargs)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, "w") as fh:
        json.dump({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                   for k, v in modes.items()}, fh, indent=2)
    return modes
```

### 6.3  API-Stabilität — wichtiger Vorbehalt

Die exakten Argumentnamen und Rückgabestrukturen der `star`-Klasse haben
sich zwischen PBjam 1.x und 2.x geändert und können sich weiter ändern.
Der Code oben zeigt die **Struktur** der Integration; die genauen
Signaturen sind vor der Umsetzung gegen die aktuellen Beispiel-Notebooks
zu prüfen:

- Session-Beispiel: https://pbjam.readthedocs.io/en/latest/Examples/example-session.html
- Star-Beispiel: https://pbjam.readthedocs.io/en/latest/Examples/example-star.html
- ModeID-Beispiel: https://pbjam.readthedocs.io/en/latest/Examples/example-modeID.html
- Peakbag-Beispiel: https://pbjam.readthedocs.io/en/latest/Examples/example-peakbag.html

Der Adapter `_extract_modes()` kapselt diese Versionsabhängigkeit bewusst
an einer einzigen Stelle — bei einem API-Wechsel muss nur diese Funktion
angepasst werden.

### 6.4  Eingabequalität

PBjam ist auf gute Startwerte angewiesen. Zwei Punkte aus den bisherigen
Pipeline-Erfahrungen:

- **Δν muss stimmen.** Ein Subharmonischer (Δν/2) verfälscht die gesamte
  Identifikation. Die Stello-Relation Δν ≈ 0.263·νmax^0.772 als
  Plausibilitätsprüfung vorschalten.
- **Kadenz beachten.** Für νmax > ~280 μHz zwingend Kurzkadenz — sonst
  liegen die Moden über Nyquist oder im Rauschen (Lehre aus dem
  gescheiterten TESS-Versuch bei KIC 10963065).

---

## 7  Vergleich: eigene Pipeline vs. PBjam

| Aspekt | Eigene Pipeline | PBjam |
|---|---|---|
| νmax, Δν | Harvey-Fit + ACF | eigener Fit (nutzt Eingabe als Prior) |
| l=0, l=2 | Positionszuordnung | asymptotische Relation, Bayesianisch |
| l=1 rein | Positionszuordnung | asymptotische Relation |
| **l=1 gemischt** | **scheitert → "unklar"** | **Kopplungsmodell, v2** |
| Frequenzgenauigkeit | Peak-Position | Lorentzfit mit Unsicherheiten |
| Rechenzeit | Sekunden | Minuten (MCMC) |
| Transparenz | vollständig einsehbar | Blackbox-artiger |
| Publikationsreife | Übungsniveau | ja |

**Empfohlene Arbeitsteilung:** Die eigene Pipeline für schnelle
Exploration, Datenbereinigung, νmax/Δν-Bestimmung und Visualisierung.
PBjam für die finale Modenidentifikation, wenn gemischte Moden
interessieren oder Ergebnisse belastbar sein sollen.

---

## 8  Umsetzungsreihenfolge

1. **PBjam installieren** und an einem bekannten Stern testen
   (z. B. dem Beispielstern aus der Dokumentation), um die API-Signaturen
   der installierten Version zu verifizieren.

2. **`pbjam_bridge.py` anpassen** — Argumentnamen und `_extract_modes()`
   an die tatsächliche Ausgabestruktur angleichen.

3. **KIC 10963065 als Kontrolle**: Der kompakte, junge Unterriese sollte
   ein sauberes l=0/1/2-Muster ohne oder mit wenigen gemischten Moden
   liefern. Vergleich mit den eigenen Pipeline-Ergebnissen
   (νmax 2175, Δν 103.34).

4. **Mulder als Zielobjekt**: Zwei avoided crossings erwartet. Wenn
   `find_avoided_crossings()` diese reproduziert, ist die Integration
   validiert.

5. **Repliziertes Échelle mit PBjam-Labels** erzeugen — jetzt sollten die
   gemischten l=1-Moden als zusammenhängende, verbogene Kurve erscheinen
   statt als graue "unklar"-Punkte.

---

## 9  Literatur

- Nielsen M.B. et al. (2021): *PBjam: A Python Package for Automating
  Asteroseismology of Solar-like Oscillators.* AJ 161, 62.
  — Grundlagenpaper, Architektur und Methode.

- Nielsen M.B. et al. (2023): *A&A 676, A117.* — Konstruktion der
  Prior-Wahrscheinlichkeitsdichten der aktuellen Version.

- Nielsen M.B., Ong J.M.J. et al. (2025): *Asteroseismology with PBjam 2.0:
  Measuring Dipole Mode Frequencies in Coupling Regimes from Main-sequence
  to Low-luminosity Red Giant Stars.* arXiv:2506.20382.
  — Die l=1-Erweiterung, zentral für gemischte Moden.

- Campante T.L. et al. (2011): *Asteroseismology of two Kepler subgiants:
  KIC 10273246 and KIC 10920273.* A&A 534, A6.
  — Referenzwerte für die Validierung.

- Dokumentation: https://pbjam.readthedocs.io/
- Quellcode: https://github.com/grd349/PBjam
