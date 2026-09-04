# Asteroseismologie-Workflow: Solar-like Oscillators mit TESS und Kepler

Dieses Dokument beschreibt den vollständigen Python-Workflow zur Analyse
pulsierender Sterne (Rote Riesen, Sonnenanaloga, Unterriesen) mit öffentlichen
Raumdaten von TESS und Kepler/K2 — von der Lichtkurve bis zu abgeleiteten
Stellar-Parametern.

---

## Projektstruktur

```
asteroseismologie/
├── .venv/                          # uv-verwaltete virtuelle Umgebung
├── data/
│   └── cache/                      # lightkurve-Cache (automatisch befüllt)
├── results/
│   ├── figures/                    # erzeugte PDFs / PNGs
│   └── tables/                     # CSV-Ausgaben (für künftige Erweiterungen)
├── asteroseismologie.py            # Hauptskript — vollständige Pipeline
├── asteroseismologie_workflow.md   # dieses Dokument
├── numax_ohne_schaetzwert.md       # Detail-Beschreibung der νmax-Pipeline
├── pyproject.toml                  # Projektdefinition für uv
└── README.md                       # Kurzreferenz
```

---

## Unterstützte Missionen

| Präfix | Mission | Kadenz | Nyquist | Geeignet für |
|--------|---------|--------|---------|------|
| `TIC …` | TESS SPOC | 2 min | ~4167 μHz | RGB (νmax < 300 μHz) und Sonnenanaloga (νmax > 1000 μHz) |
| `KIC …` | Kepler LC | 30 min | ~278 μHz | Rote Riesen (νmax < 250 μHz) |
| `KIC …` + `--exptime 60` | Kepler SC | 1 min | ~8333 μHz | Sonnenanaloga (νmax > 280 μHz) |
| `EPIC …` | K2 LC | 30 min | ~278 μHz | Rote Riesen (νmax < 250 μHz) |

Die Mission wird automatisch aus dem Ziel-Präfix erkannt. `--author` und
`--exptime` können sie überschreiben.

---

## Voraussetzungen

- Python ≥ 3.11
- [uv](https://docs.astral.sh/uv/) als Paketmanager
- Internetzugang (MAST-Download, beim ersten Lauf je Stern)

```bash
uv sync          # erzeugt .venv und installiert alle Abhängigkeiten
```

---

## Verwendung

```bash
# Standardstern: η Serpentis (TIC 234295610, Roter Riese)
uv run asteroseismologie.py

# TESS-Riese mit mehreren Sektoren
uv run asteroseismologie.py --tic "TIC 272821450" --name "eps Oph" --teff 4700 --sectors 4

# Kepler Roter Riese (Langkadenz)
uv run asteroseismologie.py --tic "KIC 4351319" --name "KIC 4351319" --teff 4800

# Kepler Sonnenanaloge (Kurzkadenz, Nyquist 8333 μHz)
uv run asteroseismologie.py --tic "KIC 10963065" --name "KIC 10963065" --teff 5900 --exptime 60 --sectors 8

# Mit Gauß-Glättung des SNR (glattere Flanken, etwas langsamer)
uv run asteroseismologie.py --tic "KIC 10963065" --name "KIC 10963065" --teff 5900 --exptime 60 --sectors 8 --gauss-smooth

# Schneller Test (halbe Frequenzauflösung)
uv run asteroseismologie.py --oversample 2

# Alle Optionen
uv run asteroseismologie.py --help
```

### Kommandozeilen-Argumente

| Argument | Standard | Beschreibung |
|---|---|---|
| `--tic` | `TIC 234295610` | Ziel-ID: `TIC …`, `KIC …` oder `EPIC …` |
| `--name` | `η Serpentis` | Sternname für Ausgabe und Dateinamen |
| `--teff` | `4970.0` | Effektive Temperatur in K (aus SIMBAD/TIC) |
| `--fmin` | `1.0` | Untere Frequenzgrenze in μHz |
| `--fmax` | `300.0` | Obere Frequenzgrenze in μHz (wird automatisch angehoben bei SC) |
| `--oversample` | `5` | Oversampling-Faktor der Frequenzauflösung |
| `--sectors` | `1` | Max. Sektoren / Quartale (mehr = besseres SNR, mehr Download) |
| `--author` | auto | Datenprovider (`SPOC`, `Kepler`, `K2`) |
| `--exptime` | auto | Belichtungszeit in s (z. B. `60` für Kepler SC, `120` für TESS SC) |
| `--gauss-smooth` | (aus) | SNR-Glättung per Gauß-FFT statt Boxcar (glattere Flanken, O(N log N)) |

---

## Pipeline-Übersicht

```
MAST-Archiv (TESS · Kepler · K2)
        │
        ▼
[1] Lichtkurve laden & bereinigen
    │  Mission-Auto-Erkennung (TIC/KIC/EPIC)
    │  Qualitätsfilter (quality_bitmask="hardest")
    │  Ausreißer-Entfernung (σ = 4)
    │  Kurzzeit-Segmente verwerfen
    │  Savitzky-Golay-Detrending (7-Tage-Fenster)
    │  RMS-basierte Rausch-Filterung (≤ 3× minimales RMS)
    │  Flux → ppm, Segmente zusammenfügen
        │
        ▼
[2] Powerspektrum berechnen (Lomb-Scargle)
    │  Automatische fmax-Anhebung bei Kurzkadenz
    │  Einseitige PSD in ppm²/μHz (Parseval-normiert)
        │
        ▼
[3] νmax vollautomatisch bestimmen (ohne Schätzwert)
    │  [3a] Harvey-Hintergrundfit im Log-Raum
    │  [3b] SNR-Spektrum: PSD / Hintergrundmodell
    │  [3c] Pre-Whitening: SNR / Boxcar-Trend (≥ 1500 μHz) → SNR-Hintergrund ≈ 1
    │  [3d] Erster Kandidat: adaptiv geglättetes SNR, Suche ab 3·b_Harvey
    │       Glättungsbreite ~ 1/40 des Suchbereichs, mind. 15 μHz (≈ 1–2 × Δν)
    │       Standard: Boxcar (O(N)); optional: Gauß via FFT (--gauss-smooth)
    │  [3e] Gauß-Fit auf adaptivem Bandpass → finales νmax
        │
        ▼
[4] Δν via Autokorrelation + Échelle-Verfeinerung
    │  Stello-Prior (±20 %)
    │  ACF-Peak im validierten Fenster
    │  Échellekolonnen-Varianz → optimales Δν
        │
        ▼
[5] Skalenrelationen
    │  Masse, Radius, log g, Leuchtkraft
    │  Referenz: Chaplin & Miglio (2013)
        │
        ▼
[6] 4-Panel-Abbildung → PDF / PNG
    │  Lichtkurve · PSD + Harvey · SNR-Spektrum · Échelle-Diagramm
```

---

## Schritt 1 · Lichtkurve laden und bereinigen

### Mission-Auto-Erkennung

Die Funktion `_mission_defaults(target_id)` liefert passende Voreinstellungen
aus dem Präfix des Ziel-IDs:

- **TIC** → TESS SPOC, 120 s Kadenz, min_span 20 d (1 Sektor ≈ 27 d)
- **KIC** → Kepler, 1800 s Kadenz, min_span 60 d (1 Quartal ≈ 90 d)
  Mit `--exptime 60`: Kepler SC, Kadenz 60 s, min_span 20 d
- **EPIC** → K2, 1800 s Kadenz, min_span 60 d (1 Kampagne ≈ 80 d)

Der min_span-Wert wird automatisch auf 20 d abgesenkt wenn `exptime ≤ 120 s`,
weil SC-Segmente deutlich kürzer als LC-Quartale sind.

### Download und Filterung

```python
search = lk.search_lightcurve(target_id, author=author, exptime=exptime)
lc_coll = search[:n_try].download_all(quality_bitmask="hardest")
```

Das `quality_bitmask="hardest"` entfernt Kadenz-Datenpunkte mit bekannten
Artifakten (kosmische Strahlen, Manöver-Events, Attitude-Tweaks u. a.).
Für Kepler SC werden typisch 20–28 % der Kadenz verworfen.

### Detrending (Savitzky-Golay)

Nach dem Qualitätsfilter wird jedes Segment einzeln mit einem
Savitzky-Golay-Filter (Polynomgrad 2) geglättet und durch seinen
Trend dividiert:

```python
dt_sec  = median(diff(lc_clean.time)) * 86400   # Kadenz in Sekunden
win_pts = 7 Tage / dt_sec                        # Fensterlänge (ungerade)
lc_clean = lc_clean.flatten(window_length=win_pts)
```

Ein 7-Tage-Fenster entfernt:
- Instrumentelle Drifts und Quarter-Anfangs-Systematics (Kepler)
- Stellare Rotationsmodulation (typisch 10–30 d für sonnenähnliche Sterne)

Das seismische Signal bleibt vollständig erhalten, weil die
Oszillationsperioden (Minuten bis wenige Stunden) weit unter dem
7-Tage-Cutoff liegen.

### Flux-Normierung und RMS-Filter

```python
flux_ppm = (flux / mean(flux) - 1.0) * 1e6
```

Anschließend werden Segmente verworfen, deren RMS mehr als das 3-fache
des niedrigsten RMS beträgt. Das entfernt Quartale mit anomal hohem
Instrumentenrauschen (z. B. nach Detektor-Wechsel) ohne manuellen Eingriff.

### Zusammenfügen (Stitch)

Mehrere Segmente werden mit `LightCurveCollection.stitch()` zusammengeführt.
Der Parameter `corrector_func=lambda lc: lc` verhindert, dass stitch()
intern erneut normiert — die Segmente sind bereits in ppm (Mittelwert ≈ 0),
sodass eine zweite Normierung ±∞-Artefakte erzeugen würde.

---

## Schritt 2 · Powerspektrum (Lomb-Scargle)

### Grundlagen

Der Lomb-Scargle-Periodogramm-Algorithmus berechnet die einseitige PSD
auch für lückenhafte Zeitreihen (Datenlücken zwischen Kepler-Quartalen,
Kosmische-Strahlen-Lücken usw.):

```
Frequenzauflösung: δf = 1 / (oversample · T_gesamt)
Frequenzbereich:   fmin … min(fmax, 0.9 · Nyquist)
```

Die PSD wird in ppm²/μHz zurückgegeben (Parseval-normiert, einseitig).

### Automatische fmax-Anpassung bei Kurzkadenz

Bei der Standard-Einstellung `fmax=300 μHz` würde das Spektrum für
TESS SC (Nyquist 4167 μHz) oder Kepler SC (Nyquist 8333 μHz) weit
unter dem Oszillationsbereich enden. Die Pipeline erkennt das automatisch:

```
Wenn Nyquist > 400 μHz  UND  fmax < 0.5 · Nyquist:
    fmax ← min(0.9 · Nyquist, 8000 μHz)
```

Für TESS SC wird fmax auf ≈ 3750 μHz, für Kepler SC auf ≈ 7500 μHz
angehoben — ohne manuellen Eingriff.

### Nyquist-Warnungen

- **fmax > Nyquist**: Aliasing-Warnung mit Empfehlung (Kepler SC)
- **νmax > 45 % Nyquist** (nach Detektion): Ergebnis-Warnung,
  Kadenz und Nyquist werden korrekt aus der Lichtkurve berechnet
  (nicht aus dem Spektralende abgeleitet, was bei manuellem `--fmax`
  zu falschen Kadenz-Angaben führen würde)

---

## Schritt 3 · νmax vollautomatisch bestimmen

Die Funktion `estimate_numax_auto()` läuft vollständig ohne externen
Schätzwert. Jeder der vier Sub-Schritte verwendet nur die Ausgabe
des vorherigen.

### 3a · Harvey-Hintergrundfit

**Physikalischer Hintergrund:**
Konvektive Granulation erzeugt ein Lorentzian-förmiges Rausch-Spektrum,
das Harvey-Modell genannt wird:

$$
B(f) = \sum_i \frac{a_i}{1 + (f/b_i)^{c_i}} + w
$$

- $a$ — Granulations-Amplitude (ppm²/μHz)
- $b$ — charakteristische Frequenz (Granulations-Knie, μHz)
- $c$ — Exponent (typisch 2–4; bei Roter Riese ≈ 2, bei Sonnenanaloga ≈ 4)
- $w$ — Weißrauschen-Plateau (Photonrauschen, Instrumentenrauschen)

Im LC-Regime genügt ein Harvey-Term. Im SC-Regime verwendet die Pipeline
zwei geordnete Terme für die langsame Aktivitäts-/Granulationskomponente
und die schnelle Granulation.

**Fit-Strategie:**
Der Fit erfolgt im Log-Raum (`log₁₀(PSD)` vs. `f`), weil das die
logarithmische Streuung des χ²-verteilten Spektrums berücksichtigt
und mehrere Dekaden Dynamikbereich gleichmäßig gewichtet.

Der Fit-Bereich ist adaptiv und wird durch die Nyquist-Frequenz gesteuert:

| Regime | Nyquist | Fit-Bereich | Begründung |
|--------|---------|-------------|------------|
| LC | ≤ 400 μHz | gesamter Frequenzbereich, logarithmisch gebinnt | Bestimmt Granulationsknie und Weißrauschplateau gemeinsam |
| SC | > 400 μHz | gesamter Frequenzbereich, logarithmisch gebinnt | Erfasst beide Knie auch bei Unterriesen mit νmax < 1500 μHz |

Im SC-Regime wird die stark streuende PSD zunächst in 100 logarithmischen
Frequenzbins durch den Median zusammengefasst. Die Medianwerte werden auf den
Erwartungswert einer exponentiell verteilten PSD korrigiert. Ein robuster
Least-Squares-Fit mit `soft_l1`-Verlust reduziert anschließend den Einfluss
des Oszillationsbuckels und einzelner Ausreißer. Die Knie sind durch
$b_1 < 200\,\mu\mathrm{Hz}$ und $b_2 > 200\,\mu\mathrm{Hz}$ geordnet;
für beide Exponenten gilt $1.5 \le c \le 6$.

Fällt der Fit fehl (z. B. zu wenig Frequenzpunkte), wird eine breit
geglättete PSD als Fallback-Hintergrund verwendet.

### 3b · SNR-Spektrum

```
SNR(f) = PSD(f) / B(f)
```

Das Harvey-Modell wird auf den gesamten Frequenzbereich extrapoliert.
Der Oszillationsexzess erscheint im SNR-Spektrum als lokaler Buckel über
dem Plateau bei SNR ≈ 1.

### 3c · Pre-Whitening

Das Harvey-Modell erfasst nicht perfekt alle Leistung bei tiefen Frequenzen
(Aktivitätsreste, SC-Granulationsübergang). Das SNR zeigt daher meist
einen breitbandigen Abfall von niedrig nach hoch, der die Lokalisierung
des Oszillationsbuckels erschwert.

Das Pre-Whitening entfernt diesen Trend:

```python
wide_bins = max(1, min(round(wide_phys / df), len(snr) // 2))
# wide_phys = max(1500 μHz, 10 × smooth_width)  —  größer als die Hüllkurve
snr_trend = uniform_filter1d(snr, size=wide_bins)  # O(N) Boxcar
snr_white = snr / snr_trend                        # SNR-Hintergrund ≈ 1
```

Nach dieser Normierung ist der Hintergrund überall ≈ 1.0; der
Oszillationsexzess tritt als relativer Buckel klar hervor.

### 3d · Erster νmax-Kandidat

Das vorverweißte SNR-Spektrum (`snr_white`) wird mit einer adaptiven
Breite geglättet, die automatisch aus dem Suchbereich bestimmt wird:

```
smooth_width = max(15 μHz, (fmax − fmin_search) / 40)
```

Diese Formel liefert für alle Regime eine Breite von etwa 1–2 × Δν:

| Regime | Beispiel | smooth_width | Faktor × Δν |
|--------|---------|-------------|-------------|
| Rote Riesen (LC) | Δν ≈ 5 μHz, fmax ≈ 250 μHz | ~15 μHz | ~3× |
| Unterriesen (SC) | Δν ≈ 50 μHz, fmax ≈ 5000 μHz | ~75 μHz | ~1.5× |
| Sonnenanaloga (SC) | Δν ≈ 100 μHz, fmax ≈ 8000 μHz | ~150 μHz | ~1.5× |

Die Glättung überbrückt die Einzelmoden und zeigt die Gauß-Hüllkurve
des Oszillationsexzesses. Das Maximum dieser geglätteten Kurve ist
der erste νmax-Kandidat.

**Glättungsfilter (wählbar):**

| Filter | Flag | Komplexität | Eigenschaft |
|--------|------|------------|-------------|
| Boxcar (`uniform_filter1d`) | Standard | O(N) | Schnell, für Peak-Finding gleichwertig |
| Gauß via FFT (`fftconvolve`) | `--gauss-smooth` | O(N log N) | Glattere Flanken, kein Gibbs |

Die Suche beginnt oberhalb des Harvey-Knies:

```
fmin_search = clip(3 · b_Harvey, 5 μHz, 0.5 · fmax)
```

Ein Plausibilitätscheck warnt, wenn der Kandidat direkt an der Untergrenze
liegt (mögliche Fehldetektion durch zu geringes SNR oder falsch kalibrierten
Harvey-b-Wert).

### 3e · Gauß-Fit (finales νmax)

Auf einem adaptiven Bandpass um den Kandidaten wird eine Gauß-Funktion
mit konstantem Hintergrund gefittet:

```
Bandpass: [0.55 · ν₀, 1.55 · ν₀]  (asymmetrisch: niedrig-freq. Seite steiler)
Modell:   A · exp(−½((f−μ)/σ)²) + c₀
```

Das Ergebnis ist $\nu_{\max}$ (= Gauß-Zentrum μ) und seine Unsicherheit
$\sigma_{\nu_{\max}}$ (= Gauß-Breite). Liegt das Gauß-Zentrum am Bandpassrand,
wird stattdessen der Kandidat aus 3c zurückgegeben.

---

## Schritt 4 · Δν via Autokorrelation und Échelle-Verfeinerung

### Stello-Prior

Die empirische Skalenrelation von Stello et al. (2009) liefert einen
a-priori-Schätzwert für Δν:

$$
\Delta\nu_{\text{prior}} \approx 0.263 \cdot \nu_{\max}^{0.772}
$$

Dieser Prior definiert das Suchfenster der Autokorrelation:
```
Suche Δν in [0.80 · Δν_prior, 1.20 · Δν_prior]
```

Ein enges Fenster (±20 %) verhindert, dass Subharmonische (Δν/2, Δν/3)
oder Aliase (Δν/√2) als Δν identifiziert werden — ein häufiger Fehler
bei niedrigem SNR.

### Autokorrelation

```
Bandpass: [νmax − max(4·Δν_prior, 0.5·νmax), νmax + max(4·Δν_prior, 0.5·νmax)]
ACF:      Autokorrelation des Bandpass-Spektrums (Mittelwert subtrahiert)
Lag-Gitter: δf · k  (k = 0, 1, 2, …)
```

Der erste ACF-Peak im validierten Fenster ist Δν_ACF.
Liegt kein Peak im Fenster, wird Δν_prior zurückgegeben.

### Échelle-Kohärenz-Verfeinerung

Auf einem feinen Δν-Gitter (±15 % um Δν_ACF, 300 Punkte) wird die
Spaltenvarianz des Échelle-Diagramms maximiert:

```
E(Δν) = Var[ ∑_bins PSD(f mod Δν)  für flo ≤ f ≤ fhi ]
```

Maximale Spaltenvarianz bedeutet, dass die Moden vertikal übereinander
stehen — das korrekte Δν. Diese Methode ist robuster als der ACF-Peak
allein, weil sie die gesamte Moden-Kohärenz über mehrere Radialordnungen
ausnutzt.

---

## Schritt 5 · Skalenrelationen

```python
SUN_NUMAX   = 3090.0   # μHz  (Chaplin & Miglio 2013)
SUN_DELTANU = 135.1    # μHz
SUN_TEFF    = 5778.0   # K

r_nu = numax   / SUN_NUMAX
r_dn = deltanu / SUN_DELTANU
r_T  = teff    / SUN_TEFF

mass   = r_nu**3  * r_dn**-4 * r_T**1.5    # M / M☉
radius = r_nu     * r_dn**-2 * r_T**0.5    # R / R☉
logg   = 4.44 + log10(r_nu) + 0.5·log10(r_T)
lumin  = radius**2 * r_T**4                # L / L☉
```

**Typische Genauigkeit:**

| Größe | Roter Riese | Sonnenanaloge |
|-------|-------------|---------------|
| Radius | ≈ ±5 % | ≈ ±3 % |
| Masse | ≈ ±10–15 % | ≈ ±5–10 % |
| log g | ≈ ±0.03 dex | ≈ ±0.02 dex |

Für Publikationen: grid-basierte Modellierung (BASTA, AMP) oder
korrigierte Skalenrelationen (Sharma et al. 2016).

---

## Schritt 6 · 4-Panel-Abbildung

Die Funktion `make_figure()` erzeugt ein 4-Panel-PDF (12 × 9 Zoll, 150 dpi):

### Panel 1: Lichtkurve (oben links)
Zeitreihe in ppm, rasterisiert für kompakte Datei-Größe.

### Panel 2: Powerspektrum (oben rechts)
- PSD in ppm²/μHz (logarithmische y-Achse, grau)
- Harvey-Hintergrundmodell B(f) (orange)
- νmax als vertikale gestrichelte Linie (blau)

### Panel 3: SNR-Spektrum (unten links)
- Vorverweißtes SNR = `snr_white` = SNR(f) / Boxcar-Trend
  (normiert, Hintergrund ≈ 1.0; hellgrau); Einzelmoden als schmale Spitzen
- Geglättetes vorverweißtes SNR (grün) — zeigt die Gauß-Hüllkurve
  mit dem Maximum bei νmax; Boxcar oder Gauß-FFT je nach `--gauss-smooth`
- νmax-Markierung (blau) liegt am Maximum der grünen Kurve

### Panel 4: Échelle-Diagramm (unten rechts)
Folgende Schritte erzeugen das Diagramm:
1. Fenster: νmax ± 5.5·Δν
2. 1D-Glättung σ = min(Δν/20, 1.0) μHz (Rauschunterdrückung, schmaler als Modenbreite)
3. SNR-Exzess = clip(SNR − 1, 0)
4. 2D-Histogramm: x = freq mod Δν, y = freq, gewichtet mit SNR-Exzess
5. 2D-Gaußglättung (σ = 1.2 Bins) für kontinuierliche Ridges
6. Stellar-Parameter als Text-Box

Ridges (vertikale helle Linien) entsprechen den Schwingungs-Moden:
- `l=0` (radiale Moden): dominante rechte Säule
- `l=1` (dipolare Moden): linke Säule, im RGB oft als gemischte Moden
- `l=2` (quadrupolare Moden): eng neben l=0

Bei Windows-Datei-Sperren (PDF geöffnet im Viewer) wird automatisch
als PNG gespeichert.

---

## Bekannte Einschränkungen und Fallstricke

### Harvey-Fit bei aktiven Sternen
Für magnetically aktive Sonnenanaloga (Rotationsperiode 15–30 d) können
die Rotationsmodulation und Aktivitätssignal das SNR-Plateau anheben
und den Harvey-Fit erschweren. Die 7-Tage-Detrending-Fenster hilft,
entfernt aber keine Aktivitäts-Signale mit Perioden < 7 d.

### Harvey-Modell bei SC-Daten
Die Pipeline verwendet für SC-Daten zwei Harvey-Terme. Sehr starke oder
ungewöhnlich breite Oszillationsbuckel können den robusten Hintergrundfit
weiterhin beeinflussen; für Präzisionsanalysen ist ein gemeinsamer
Bayes-Hintergrund- und Hüllkurvenfit vorzuziehen.

### Kepler SC: hohe Ausdünnung der Daten
Bei Kepler SC werden mit dem `hardest`-Qualitätsmask 20–28 % der
Kadenz verworfen. Das verringert den effektiven Duty Cycle auf
~50 % und erhöht das Rauschen im Lomb-Scargle-Spektrum.

### νmax nahe Nyquist
Das Ergebnis wird unsicher, wenn νmax > 45 % der Nyquist-Frequenz ist.
Die Pipeline warnt in diesem Fall und empfiehlt Kurzkadenz-Daten.

---

## Empfohlene Einstiegssterne

| Stern | ID | Teff (K) | νmax (μHz) | Δν (μHz) | Mission | Besonderheit |
|-------|-----|---------|-----------|---------|---------|----------|
| η Serpentis | TIC 234295610 | 4970 | ~173 | ~13.4 | TESS SC | Standardstern, hell, einfach |
| ε Ophiuchi | TIC 272821450 | 4700 | ~59 | ~5.3 | TESS SC | Klarer Oszillationsbuckel |
| KIC 4351319 | KIC 4351319 | 4800 | ~47 | ~4.7 | Kepler LC | Viele Referenzwerte |
| KIC 10963065 | KIC 10963065 | 5900 | ~2000 | ~103 | Kepler SC | Sonnenanaloge, aktiv |

---

## Weiterführende Schritte

### Verbesserte νmax-Bestimmung
- **Gemeinsamer Background-/Oszillationsfit** statt sequenzieller Anpassung
- **Bayes'scher Background-Fit** (z. B. mit DIAMONDS/FAMED)
- **Ensembles**: νmax aus mehreren Sektoren/Quartalen gemittelt

### Verbesserte Δν-Bestimmung
- **Universal Pattern** (Mosser et al. 2011): berücksichtigt
  asymptotische Korrekturen zweiter Ordnung
- **MCMC-Fitting** individueller Moden-Peaks

### Publizierbare Stellar-Parameter
- **BASTA** (BAyesian STellar Algorithm): Gitter-basierte Modellierung
- **AMP** (Asteroseismic Modeling Portal): webbasiert
- Korrigierte Skalenrelationen: Sharma et al. (2016, ApJ 822, 15)

### Weitere Datenquellen

| Archiv | Inhalt |
|--------|--------|
| MAST (mast.stsci.edu) | TESS, Kepler, K2 — via lightkurve |
| KASOC (kasoc.phys.au.dk) | Aufbereitete Kepler-LC für Asteroseismologie |
| Gaia DR3 | Variabilitätskatalog, Teff-Referenzen |

---

## Literatur

- **Chaplin & Miglio (2013**, ARA&A 51, 353) — Solar-like oscillations, Review
- **Stello et al. (2009**, ApJ 700, 1589) — νmax–Δν Skalenrelation
- **Mosser & Appourchaux (2009**, A&A 508, 877) — Autokorrelationsmethode
- **Harvey (1985**, ESA SP-235) — Granulations-Hintergrundmodell
- **Sharma et al. (2016**, ApJ 822, 15) — Korrigierte Skalenrelationen
- **Kjeldsen & Bedding (1995**, A&A 293, 87) — Ursprüngliche Skalenrelationen


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
