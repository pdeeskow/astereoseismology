# Projektvorschlag: Seismischer Atlas der Sternentwicklung

Ein abgegrenztes Vergleichsprojekt, das mit der bestehenden Pipeline die
typischen Modenkonstellationen entlang der Entwicklungssequenz
nachvollziehbar macht.

---

## 1  Leitidee

Bisher wurden Sterne einzeln analysiert. Der Erkenntnissprung entsteht
beim **systematischen Vergleich**: Dieselbe Pipeline, dieselbe
Darstellung, angewandt auf eine Sequenz von Sternen, die zusammen einen
Entwicklungsweg abbilden.

Die Leitfrage:

> Wie verändert sich die Modenstruktur eines Sterns, während er sich von
> der Hauptreihe über das Unterriesenstadium zum Roten Riesen entwickelt —
> und woran erkennt man das Stadium allein am Powerspektrum?

Der Reiz liegt darin, dass jede der vier bisher separat behandelten
Größen (Δν, νmax, δν₀₂, gemischte Moden) sich entlang dieser Sequenz
**systematisch und vorhersagbar** verändert. Der Atlas macht diese
Systematik in einem Bild sichtbar.

---

## 2  Umfang und Abgrenzung

**Was das Projekt ist:**
- 6–8 Sterne, alle mit Kepler-Daten und publizierten Referenzwerten
- Eine Pipeline, identisch angewandt, alle Ergebnisse gecacht
- Vier vergleichende Auswertungen (Abschnitt 5)
- Ergebnis: ein Atlas-PDF plus Datentabelle

**Was das Projekt nicht ist:**
- Keine Sternmodellierung (MESA/GYRE) — das wäre ein Folgeprojekt
- Keine neuen Methoden — die Pipeline bleibt wie sie ist
- Keine Stapelverarbeitung hunderter Sterne

**Geschätzter Aufwand:** Datenbeschaffung und Durchläufe ein bis zwei
Abende (PBjam-Sampling läuft im Hintergrund), Auswertung und Layout ein
weiterer. Der Code existiert bereits — es geht um Anwendung und
Zusammenstellung.

---

## 3  Die Sternsequenz

Ausgewählt nach drei Kriterien: Kepler-Daten vorhanden, publizierte
Referenzwerte zur Validierung, und gemeinsam den Entwicklungsweg
abdeckend.

| # | Stern | Stadium | νmax (μHz) | Δν (μHz) | Kadenz | Status |
|---|---|---|---|---|---|---|
| 0 | **Sonne** | Hauptreihe (Referenz) | 3090 | 135.1 | — | Literaturwerte |
| 1 | **16 Cyg A** (KIC 12069424) | Hauptreihe, sonnenähnlich | ~2200 | ~102 | kurz | zu verifizieren |
| 2 | **16 Cyg B** (KIC 12069449) | Hauptreihe, Begleiter | ~2550 | ~117 | kurz | zu verifizieren |
| 3 | **KIC 10963065** | früher Unterriese | 2145 | 103.11 | kurz | ✓ bereits analysiert |
| 4 | **KIC 10920273** (Scully) | Unterriese | 990 | 57.27 | kurz | Campante 2011 |
| 5 | **KIC 10273246** (Mulder) | entwickelter Unterriese | 838 | 48.89 | kurz | ✓ bereits analysiert |
| 6 | **KIC 4351319** | niedrigleucht. RGB | ~47 | ~4.7 | lang | noch nicht analysiert |
| 7 | *RGB oder Red Clump* | Roter Riese | ~30–50 | ~4 | lang | aus Yu+2018 wählen |

Zwei Sterne sind bereits durchgerechnet. KIC 4351319 ist der erste neue
RGB-Lauf; insgesamt betrifft der Zusatzaufwand fünf bis sechs neue Sterne.

**Hinweis zu 16 Cyg A/B:** Das ist ein Doppelsternsystem mit
gemeinsamer Entstehung, also gleichem Alter und gleicher
Anfangszusammensetzung. Beide sind sehr gut untersucht. Als Paar sind
sie besonders lehrreich, weil Unterschiede in ihren Spektren allein auf
die Massendifferenz zurückgehen. Referenzwerte vor der Analyse aus der
Literatur verifizieren.

**Hinweis zur Sonne:** Sie kommt ohne eigene Analyse in den Atlas —
einfach als Ankerpunkt mit Literaturwerten. Sie ist der bestvermessene
Stern überhaupt und verankert die Sequenz am jungen Ende.

---

## 4  Ablauf je Stern

Der Durchlauf ist für alle identisch — genau das macht den Vergleich
belastbar.

```
1. Lichtkurve laden        (Kadenz nach νmax-Prior: >280 μHz → kurz)
2. PSD via Lomb-Scargle
3. Harvey-Fit → SNR-Spektrum
4. νmax (schätzwertfrei, 4-Stufen-Pipeline)
5. Δν via ACF + Stello-Prior als Suchfenster
6. Δν im Échelle visuell nachjustieren (Ridges senkrecht)
7. PBjam: Modenidentifikation + Peakbagging
8. Qualitätsfilter (Prior/Posterior-Verhältnis > 2)
9. δν₀₂ aus l=0/l=2-Paaren messen
10. Avoided crossings im l=1-Ridge zählen
11. Ergebnisse cachen als JSON
```

**Wichtig:** Schritt 11 nicht überspringen. Bei acht Sternen und
MCMC-Sampling lohnt sich das Caching sofort — sonst rechnet jeder
Plot-Durchlauf alles neu.

Ergebnisstruktur je Stern:

```python
{
  "id": "KIC10273246",
  "label": "Mulder",
  "stage": "entwickelter Unterriese",
  "numax": [856.4, 25.0],
  "dnu": [48.44, 0.15],
  "d02": [4.1, 0.4],
  "logg": 3.90,
  "radius": 2.24,
  "mass": 1.45,
  "n_avoided_crossings": 2,
  "modes": {"freq": [...], "l": [...], "n": [...], "height": [...]},
  "reference": {"numax": 838, "dnu": 48.89, "d02": 4.40}
}
```

---

## 5  Die vier Auswertungen

Hier liegt der eigentliche Erkenntnisgewinn. Jede Auswertung macht eine
andere Systematik sichtbar.

### 5.1  Échelle-Atlas

Ein Raster aller Sterne, nach νmax sortiert — jedes Panel ein Échelle
mit farblich getrennten l-Werten, alle im gleichen Stil.

**Was sichtbar wird:** Die Ridges werden von links nach rechts (junge zu
alte Sterne) zunehmend unregelmäßiger. Bei den Hauptreihensternen stehen
alle drei Ridges gerade und sauber getrennt; beim Unterriesen verbiegt
sich l=1; beim Roten Riesen wird der l=1-Bereich zunehmend von
gemischten Moden bevölkert.

**Umsetzung:** `matplotlib` GridSpec, 2×4 Panels, gemeinsame Legende.
Jedes Panel mit x-Achse 0…Δν (nicht repliziert, sonst zu unruhig) und
y-Achse normiert auf (ν − νmax)/Δν, damit alle Sterne vergleichbar
skaliert erscheinen.

### 5.2  C-D-Diagramm (Christensen-Dalsgaard)

δν₀₂ gegen Δν, doppelt logarithmisch. Das ist das seismische Pendant
zum HR-Diagramm.

**Was sichtbar wird:** Die Sterne ordnen sich entlang einer Sequenz.
Δν gibt im Wesentlichen die mittlere Dichte (also Radius/Masse), δν₀₂
das Alter. Sterne gleicher Masse aber verschiedenen Alters liegen auf
absteigenden Linien — genau darum ist das Diagramm für die
Altersbestimmung so nützlich.

**Bonus:** 16 Cyg A und B müssen im C-D-Diagramm auf *derselben*
Alterslinie liegen, weil sie gleich alt sind. Das ist eine schöne
Konsistenzprüfung der eigenen Messung.

### 5.3  Skalenrelations-Test

Zwei Teilplots:

**(a) Δν gegen νmax** mit der Stello-Relation Δν = 0.263·νmax^0.772 als
Kurve. Alle acht Sterne sollten eng an der Kurve liegen. Abweichungen
zeigen, wo die Relation an ihre Grenzen kommt (bei Hauptreihennähe
leicht, wie du bei KIC 10963065 gesehen hast).

**(b) Gemessene gegen publizierte Werte** für νmax, Δν, M, R. Eine
Identitätslinie mit den Datenpunkten daneben. Das ist die
Validierung der ganzen Pipeline in einem Bild — und wissenschaftlich
das Ehrlichste, was man zeigen kann.

### 5.4  Zensus der gemischten Moden

Zahl der detektierten avoided crossings gegen log g (oder gegen νmax).

**Was sichtbar wird:** Die erwartete Systematik — Hauptreihensterne
zeigen keine, Unterriesen einige, Rote Riesen viele. Der Grund ist die
Kernkontraktion: Je dichter der Kern, desto höher die
g-Moden-Frequenzen, desto mehr Kreuzungen fallen in das beobachtete
Frequenzfenster.

**Zusatzgröße, falls machbar:** Das Periodenspacing ΔΠ₁ der
g-dominierten l=1-Moden. Es ist die g-Moden-Entsprechung zu Δν und
misst direkt die Kernstruktur. Bei Roten Riesen unterscheidet es sogar
wasserstoff- von heliumbrennenden Sternen (RGB vs. Red Clump) — eine
der elegantesten Anwendungen der Asteroseismologie überhaupt.

---

## 6  Deliverables

| Datei | Inhalt |
|---|---|
| `atlas_echelle.pdf` | Raster aller Échelle-Diagramme |
| `atlas_cd_diagram.pdf` | C-D-Diagramm mit allen Sternen |
| `atlas_scaling.pdf` | Skalenrelationen + Validierung gegen Literatur |
| `atlas_mixed_modes.pdf` | Zensus gemischter Moden |
| `results_table.csv` | Alle gemessenen Größen, maschinenlesbar |
| `results_table.md` | Dieselbe Tabelle formatiert, mit Literaturvergleich |
| `README.md` | Methodenbeschreibung, Reproduktionsanleitung |

---

## 7  Projektstruktur

```
seismischer-atlas/
├── config/
│   └── targets.yaml          # Sternliste mit Prioren und Referenzwerten
├── cache/
│   └── <KIC>/                # PBjam-Ergebnisse, JSON + Rohspektren
├── pipeline/                 # bestehender Code, unverändert
├── atlas/
│   ├── run_all.py            # Batch über alle Sterne
│   ├── plot_echelle_grid.py
│   ├── plot_cd_diagram.py
│   ├── plot_scaling.py
│   └── plot_mixed_census.py
├── results/
└── README.md
```

Die Sternliste als YAML hält Konfiguration und Code getrennt:

```yaml
targets:
  - id: KIC10273246
    label: Mulder
    stage: entwickelter Unterriese
    cadence: short
    numax_prior: [838, 50]
    dnu_prior: [48.89, 0.5]
    teff: [6150, 100]
    reference:
      numax: 838
      dnu: 48.89
      d02: 4.40
      source: "Campante et al. 2011, A&A 534, A6"
```

---

## 8  Umsetzungsreihenfolge

**Phase 1 — Konsolidieren (kein neuer Stern)**
Die zwei bereits analysierten Sterne KIC 10273246 und KIC 10963065 in die
neue Struktur überführen, Caching einbauen und die verfügbaren Plot-Skripte
an diesen beiden testen.

**Phase 2 — Sequenz vervollständigen**
Zuerst KIC 4351319 als neuen Long-Cadence-RGB-Fall analysieren und separat
validieren. Danach folgen Scully (analog zu Mulder, gut dokumentiert),
16 Cyg A/B (hohe Frequenzen, aber exzellentes SNR) und zuletzt der zweite
Rote Riese.

**Phase 3 — Auswertung und Layout**
Die vier Vergleichsplots erstellen, Tabelle gegen Literatur prüfen,
README schreiben.

---

## 9  Was das Projekt lehrt

Über die einzelnen Plots hinaus macht der Atlas drei Dinge greifbar,
die bei Einzelanalysen unsichtbar bleiben:

**Die Systematik der Skalenrelationen.** Dass Δν und νmax über zwei
Größenordnungen einer engen Relation folgen, ist beeindruckender, wenn
man es an eigenen Messungen sieht statt in einem Lehrbuchdiagramm.

**Die Grenzen der eigenen Methode.** Bei acht Sternen unterschiedlicher
Helligkeit, Kadenz und Frequenzlage zeigt sich, wo die Pipeline robust
ist und wo nicht. Der gescheiterte TESS-Versuch bei KIC 10963065 war
schon so ein Lernmoment — der Atlas macht solche Grenzen systematisch
sichtbar statt zufällig.

**Die Physik der Kernentwicklung.** Der Übergang von "keine gemischten
Moden" zu "verbogener l=1-Ridge" zu "dichter Wald gemischter Moden" ist
die direkte Beobachtung der Kernkontraktion nach dem
Wasserstoffbrennen. Das in einer eigenen Bilderreihe zu sehen, ist etwas
anderes als es zu lesen.

---

## 10  Mögliche Folgeprojekte

Falls der Atlas Lust auf mehr macht:

- **Modellvergleich mit MESA/GYRE:** Für einen der Sterne ein
  Modellgitter rechnen und das Alter aus der Frequenzanpassung
  bestimmen. Das ist der Schritt von "Parameter messen" zu "Stern
  verstehen".

- **Rotationssplitting:** Bei ausreichendem SNR spalten die m-Komponenten
  auf. Aus der Aufspaltung gemischter Moden lässt sich die
  Kernrotation getrennt von der Hüllenrotation messen — eines der
  spektakulärsten Ergebnisse der Kepler-Ära (Rote Riesen rotieren im
  Kern viel schneller als in der Hülle).

- **TESS-Erweiterung:** Die hellen Roten Riesen aus der Frühphase des
  Projekts (ε Oph, η Ser) mit Halo-Photometrie erschließen und in den
  Atlas aufnehmen. Das verbindet die Kepler-Sequenz mit den hellen
  Nachbarsternen.
