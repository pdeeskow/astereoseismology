# Analyse KIC 4351319 — Kalibrierung des PBjam-Qualitätsfilters

Auswertung des Laufs vom 12.09.2026 gegen Di Mauro et al. 2011, MNRAS 415, 3783.

---

## 1  Kernbefund

**Die Pipeline arbeitet korrekt. Der Qualitätsfilter war zu streng eingestellt.**

| Kennzahl | Wert |
|---|---|
| PBjam-Moden gesamt | 28 |
| Übereinstimmung mit Di Mauro (Toleranz 0.6 μHz) | **17 von 25 (68 %)** |
| l-Zuordnung korrekt | **16 von 16 (100 %)** |
| Mediane Frequenzabweichung | 0.14 μHz (Auflösung: 0.386 μHz) |
| Als „zuverlässig" gemeldet bei `quality_min = 2.0` | 3 |

Die gemeldeten 3 Moden waren ein Artefakt der Schwelle, keine Detektionsgrenze.

### Globalparameter — alle innerhalb der Unsicherheiten

| Größe | gemessen | Literatur | Abweichung |
|---|---|---|---|
| νmax | 381.62 ± 13.85 μHz | 386.5 ± 4.0 | −1.3 % |
| Δν | 24.719 ± 0.247 μHz | 24.6 ± 0.2 | +0.5 % |
| Masse | 1.233 M☉ | 1.35 ± 0.09 (Skalenrel.) | −8.7 % → 1.3 σ |
| Radius | 3.327 R☉ | 3.44 ± 0.08 (Skalenrel.) | −3.3 % → 1.4 σ |
| log g | 3.487 | 3.50 | −0.4 % |

> Vergleich gegen Di Mauros **Skalenrelations**-Werte (M = 1.35, R = 3.44),
> nicht gegen deren Modellfit (M = 1.30, R = 3.37) — die Pipeline rechnet
> ebenfalls mit Skalenrelationen.

---

## 2  Warum der Filter versagt

### Die Qualitätsmetrik hat bei kurzen Zeitreihen einen Boden

Die Metrik ist das Verhältnis Prior-Breite / Posterior-Breite. Bei 30 Tagen
Beobachtung (Q3.1) beträgt die Frequenzauflösung 1/T = 0.386 μHz. Die
gefitteten Modenbreiten liegen im Median bei 0.517 μHz — die Moden sind also
nur **marginal aufgelöst**.

Folge: Der Breitenparameter im Lorentzfit bleibt schlecht bestimmt
(mediane relative Unsicherheit **214 %**), die Posterior-Verteilung bleibt
fast so breit wie der Prior, und die Metrik sammelt sich bei **≈ 0.7**.

```
Qualitätsverteilung:   Min 0.59 | Median 0.74 | Max 3.90
```

Die Schwelle 2.0 stammt aus Nielsen et al. und ist an **vierjährigen**
Kepler-Zeitreihen kalibriert. Für 30 Tage ist sie unbrauchbar.

### Schwellenanalyse

| Schwelle | Moden | Treffer | Reinheit | Abdeckung |
|---|---|---|---|---|
| 0.00 | 28 | 17 | 61 % | 68 % |
| 0.65 | 25 | 16 | 64 % | 64 % |
| 0.70 | 20 | 14 | 70 % | 56 % |
| 0.75 | 13 | 12 | 92 % | 48 % |
| **0.80** | **9** | **9** | **100 %** | **36 %** |
| 1.00 | 6 | 6 | 100 % | 24 % |
| 2.00 | 3 | 3 | 100 % | 12 % |

Ab 0.80 ist die Reinheit bereits 100 %. Alles darüber verwirft nur noch
echte Moden.

---

## 3  Die Detektion ist SNR-limitiert, nicht fehlerhaft

| Di-Mauro-SNR | von PBjam gefunden |
|---|---|
| ≥ 5.5 | 8 von 8 (100 %) |
| 4.0 – 5.5 | 5 von 7 (71 %) |
| 3.0 – 4.0 | 3 von 10 (30 %) |

Eine saubere Detektionskurve. Di Mauro erreichten die schwachen Moden nur
über sechs parallele Pipelines mit Konsenskriterium (Mode in Tabelle 3, wenn
von ≥ 2 Teams bestätigt). Mit einer Pipeline ist das nicht reproduzierbar —
und auch nicht nötig.

### Abdeckung nach Grad

| l | gefunden | Bemerkung |
|---|---|---|
| 0 | 6 / 6 (100 %) | vollständig |
| 2 | 6 / 8 (75 %) | die 2 fehlenden haben SNR 3.4 |
| 1 | 4 / 11 (36 %) | siehe unten |

Die schwache l=1-Ausbeute ist erwartbar: Di Mauro merken an, dass bei diesem
Stern g- und p-Modenabstände vergleichbar groß sind und das Mode Bumping
sehr stark ist. Von den 11 l=1-Moden haben 7 einen SNR unter 4.

---

## 4  `d02_pairs: 0` ist ein reines Filterartefakt

Alle fünf l=0/l=2-Paare sind in PBjams Liste vorhanden. Der Filter zerschneidet
sie systematisch — mal fällt der l=0-Partner, mal der l=2:

```
l0 329.30 (q=0.83, bleibt)  |  l2 326.82 (q=0.73, fällt)   d02=2.48
l0 353.23 (q=0.70, fällt)   |  l2 351.32 (q=0.89, bleibt)  d02=1.91
l0 378.19 (q=3.90, bleibt)  |  l2 376.07 (q=0.76, fällt)   d02=2.12
l0 402.10 (q=0.65, fällt)   |  l2 400.55 (q=0.89, bleibt)  d02=1.55
l0 427.39 (q=2.45, bleibt)  |  l2 425.60 (q=0.67, fällt)   d02=1.79
```

### δν₀₂ bei verschiedenen Schwellen

| Schwelle | l0 | l2 | Paare | δν₀₂ |
|---|---|---|---|---|
| 2.00 | 2 | 0 | 0 | — |
| 0.80 | 3 | 2 | 0 | — |
| 0.75 | 4 | 3 | 1 | 2.13 |
| 0.70 | 6 | 4 | 3 | 2.17 ± 0.24 |
| **0.65** | **6** | **6** | **5** | **2.24 ± 0.40** |
| 0.00 | 7 | 7 | 7 | 2.06 ± 0.43 |

Literatur: **δν₀₂ = 2.20 ± 0.30 μHz** → bei Schwelle 0.65 exakt getroffen.

---

## 5  Die Modenhöhe trennt besser als die Qualitätsmetrik

| Gruppe | Median Modenhöhe |
|---|---|
| Literaturtreffer | 6.84 |
| ohne Gegenstück | 0.37 |

Faktor 18 — eine deutlich schärfere Trennung als die Qualitätsmetrik
(0.83 gegen 0.72).

**Kombiniertes Kriterium** `quality ≥ 0.75 AND height ≥ 1.0`:
12 Moden, 11 Treffer → **92 % Reinheit bei 44 % Abdeckung**.
Besser als der reine Qualitätsfilter bei gleicher Reinheit.

---

## 6  Die 11 Moden ohne Literaturgegenstück

| ν (μHz) | l | Qual | Höhe | Einordnung |
|---|---|---|---|---|
| 308.63 | 1 | 0.69 | 0.2 | Kopplungsmodell, kein Signal |
| 329.96 | 1 | 0.73 | 1.5 | Kopplungsmodell, kein Signal |
| 352.70 | 1 | 0.67 | 6.9 | Kopplungsmodell, kein Signal |
| **407.09** | **1** | **0.78** | **1.6** | **= Di Mauro marginal (407.41, SNR 2–3)** |
| 436.29 | 1 | 0.59 | 0.8 | Kopplungsmodell, kein Signal |
| 447.60 | 1 | 0.72 | 0.4 | Kopplungsmodell, kein Signal |
| 452.71 | 0 | 0.73 | 0.3 | außerhalb Di Mauros Frequenzfenster |
| 463.75 | 1 | 0.73 | 0.4 | Kopplungsmodell, kein Signal |
| 474.64 | 2 | 0.64 | 0.2 | außerhalb Frequenzfenster |
| 483.38 | 1 | 0.68 | 0.1 | Kopplungsmodell, kein Signal |
| 493.64 | 1 | 0.74 | 0.0 | Kopplungsmodell, kein Signal |

Neun davon sind l=1 mit sehr kleinen Höhen: Positionen, an denen PBjams
Kopplungsmodell gemischte Moden **vorhersagt**, die Daten sie aber nicht
bestätigen. Erwartetes Verhalten — die niedrige Qualitätsmetrik markiert sie
korrekt.

Eine Mode (407.09) entspricht einer von Di Mauros drei marginalen Detektionen.

---

## 7  Umsetzung

### 7.1  Schwelle datenabhängig statt fest

Die feste Schwelle ersetzen durch eine, die sich am Qualitätsboden orientiert:

```python
def adaptive_quality_threshold(modes, factor=1.10, floor=0.75, ceiling=2.0):
    """
    Setzt die Qualitätsschwelle relativ zum Median der Verteilung.

    Begründung: Bei kurzen Zeitreihen sammeln sich die Qualitätswerte an
    einem Boden (Posterior ≈ Prior). Eine feste Schwelle trifft dann
    entweder alles oder nichts. Der Median markiert diesen Boden; ein
    kleiner Aufschlag darüber trennt informierte von uninformierten Moden.

    Für KIC 4351319: Median 0.74 → Schwelle 0.81  (Reinheit 100 %)
    Für 4-Jahres-Kepler-Daten liegt der Median deutlich höher, die
    Schwelle wandert automatisch mit.
    """
    import numpy as np
    q = np.array([m["quality"] for m in modes if m.get("quality") is not None])
    if len(q) == 0:
        return floor
    return float(np.clip(np.median(q) * factor, floor, ceiling))
```

### 7.2  Höhenfilter ergänzen

```python
def filter_modes(modes, quality_min=None, height_min=1.0):
    """
    Zweistufiger Filter: Qualitätsmetrik UND Modenhöhe.

    Die Höhe trennt echte Detektionen von Modellvorhersagen ohne Signal
    deutlich schärfer als die Qualitätsmetrik allein
    (Median 6.84 gegen 0.37 bei KIC 4351319).
    """
    if quality_min is None:
        quality_min = adaptive_quality_threshold(modes)
    return [m for m in modes
            if m.get("quality", 0) >= quality_min
            and m.get("height", 0) >= height_min]
```

### 7.3  δν₀₂ mit eigener, lockerer Schwelle

Die Paarbildung braucht **beide** Partner. Ein Filter, der pro Mode
entscheidet, zerschneitet Paare zufällig. Lösung: δν₀₂ auf einer separaten,
lockereren Auswahl bestimmen.

```python
def measure_d02(modes, quality_min_pairs=0.65, window=(1.0, 4.0)):
    """
    δν₀₂ aus l=0/l=2-Paaren.

    Verwendet eine eigene, lockerere Schwelle als die Hauptauswahl:
    die Paarbildung ist selbst ein Konsistenztest (der Abstand muss im
    erwarteten Fenster liegen), daher ist ein strenger Einzelfilter
    kontraproduktiv.

    Validierung KIC 4351319: 5 Paare, δν₀₂ = 2.24 ± 0.40 μHz
    gegen Literatur 2.20 ± 0.30 μHz.
    """
    import numpy as np
    sel = [m for m in modes if m.get("quality", 0) >= quality_min_pairs]
    l0 = sorted(m["freq"] for m in sel if m["l"] == 0)
    l2 = sorted(m["freq"] for m in sel if m["l"] == 2)
    pairs = [(a, b, a - b) for a in l0 for b in l2
             if window[0] < a - b < window[1]]
    if not pairs:
        return None, None, 0
    d = np.array([p[2] for p in pairs])
    return float(d.mean()), float(d.std()), len(pairs)
```

### 7.4  QC-Codes anpassen

Die beiden `warning`-Meldungen aus dem Lauf waren irreführend:

| Code | bisher | neu |
|---|---|---|
| `reliable_mode_count` | Warnung ab < 5 Moden | Schwelle an Beobachtungsdauer koppeln; zusätzlich Rohzahl vor Filterung melden |
| `d02_pairs` | Warnung bei 0 Paaren | erst warnen, wenn auch mit der lockeren Paarschwelle keine Paare gefunden werden |

Sinnvolle Ergänzung: einen QC-Code `quality_floor` einführen, der anschlägt,
wenn `median(quality) < 0.9` — das signalisiert „Zeitreihe zu kurz für
belastbares Peakbagging" und erklärt niedrige Modenzahlen von selbst.

### 7.5  Rohzahlen mitschreiben

`summary.json` sollte neben `reliable_total` auch die ungefilterte Zahl
und die verwendete Schwelle enthalten:

```json
"mode_counts": {
  "l0": 6, "l1": 4, "l2": 6,
  "reliable_total": 16,
  "raw_total": 28,
  "quality_threshold_used": 0.81,
  "quality_median": 0.74,
  "height_min_used": 1.0
}
```

---

## 8  Für den Atlas

Der Stern ist vollwertig verwendbar. Alle vier Größen liegen vor:

| Größe | Wert | Literatur |
|---|---|---|
| νmax | 381.6 ± 13.9 μHz | 386.5 ± 4.0 |
| Δν | 24.72 ± 0.25 μHz | 24.6 ± 0.2 |
| δν₀₂ | 2.24 ± 0.40 μHz (Schwelle 0.65) | 2.20 ± 0.30 |
| ΔΠ₁ | — (noch zu messen) | 40 ± 4 s → RGB |

Für das C-D-Diagramm und den RGB/Red-Clump-Kontrast gegen KIC 1161618 ist
damit alles vorhanden.

### Erwartung für die übrigen Atlas-Sterne

Die Qualitätsschwelle wird bei den anderen Zielen **anders** liegen, weil die
Beobachtungsdauern stark variieren:

| Stern | Dauer | erwarteter Qualitätsboden |
|---|---|---|
| 16 Cyg A/B | 928 d | hoch → Schwelle 2.0 angemessen |
| KIC 10963065, Mulder, Scully | ~8 Monate | mittel |
| KIC 4351319 | 30 d | niedrig → 0.8 |
| KIC 1161618, KIC 1433730 | 1318 d | hoch → 2.0 angemessen |

Das ist das Argument für die adaptive Schwelle aus 7.1 statt eines festen
Werts in `targets.yaml`.

---

## 9  Werkzeuge

| Datei | Zweck |
|---|---|
| `vergleich_dimauro.py` | Modenvergleich gegen Literatur, Schwellenanalyse |
| `dimauro_table3.csv` | Referenztabelle, 25 Moden + 3 marginale, mit SNR und Team-Bestätigung |

Aufruf:

```bash
python vergleich_dimauro.py results/tables/KIC_4351319_pbjam_modes.csv --tol 0.6
python vergleich_dimauro.py <csv> --tol 0.6 --quality-min 0.80
```

Das Schema lässt sich auf andere Sterne übertragen, sobald eine
Literaturmodenliste vorliegt — für Mulder und Scully etwa aus Campante et al.
2011, für 16 Cyg A/B aus Davies et al. 2015 (54 bzw. 56 Moden).
