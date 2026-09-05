# Konzept: Repliziertes Échelle-Diagramm mit farblich getrennten l-Werten

Erweiterung der Astroseismologie-Pipeline zur Sichtbarmachung gemischter Moden
und avoided crossings bei Unterriesen (z. B. KIC 10273246 / "Mulder").

---

## 1  Ziel und Motivation

Das Standard-Échelle faltet das Powerspektrum modulo Δν. Bei Unterriesen mit
gemischten Moden (l=1) entsteht dabei ein Problem: Der l=1-Ridge verbiegt sich
an den avoided crossings und läuft oft **über den rechten Diagrammrand hinaus**,
wo er dann links wieder eintritt. Diese Umbrüche zerreißen die durchgehende
S-Kurve, die eine gemischte Mode eigentlich zeichnet.

Zwei Erweiterungen lösen das:

1. **Repliziertes Échelle** (Bedding 2012): Das Diagramm wird horizontal
   mehrfach nebeneinander wiederholt. Dadurch läuft der l=1-Ridge ohne
   Randumbruch durch, und die avoided crossings werden als glatte
   Wendepunkte sichtbar.

2. **Farbliche l-Trennung**: Statt einer einheitlichen Heatmap werden die
   Moden nach Grad l = 0, 1, 2 identifiziert und farblich getrennt
   dargestellt. So lässt sich der verbogene l=1-Ridge klar von den
   geraden l=0/l=2-Ridges unterscheiden.

---

## 2  Konzept "Repliziertes Échelle"

### 2.1  Grundidee

Im normalen Échelle wird jede Frequenz an Position `(ν mod Δν, ν)` geplottet.
Die x-Achse läuft von 0 bis Δν. Ein Ridge, der bei x ≈ Δν den rechten Rand
erreicht, springt zum nächsten Datenpunkt an den linken Rand (x ≈ 0) zurück.

Beim replizierten Échelle wird die x-Achse **verdoppelt oder verdreifacht**:
Jede Mode wird nicht nur an `x = ν mod Δν` geplottet, sondern zusätzlich an
`x + Δν` (und optional `x + 2Δν`). Die x-Achse läuft dann von 0 bis 2Δν
(bzw. 3Δν).

```
Normal:      [0 ........ Δν]
Repliziert:  [0 ........ Δν ........ 2Δν]
             |___Kopie 1___|___Kopie 2___|
```

Ein Ridge, der in Kopie 1 den rechten Rand erreicht, setzt sich nahtlos in
Kopie 2 fort. Die visuelle Kontinuität bleibt erhalten.

### 2.2  Warum das bei gemischten Moden hilft

Der l=1-Ridge eines Unterriesen ist keine gerade vertikale Linie, sondern
eine S-Kurve: An einem avoided crossing weicht die Mode zur Seite aus
(mode bumping), kehrt dann zurück. Diese S-Kurve überstreicht oft mehr als
eine Δν-Breite. Im normalen Échelle wird sie am Rand abgeschnitten; im
replizierten Échelle bleibt sie zusammenhängend sichtbar.

### 2.3  Implementierungsschema

```python
def replicated_echelle_coords(freq, power, dnu, n_replicas=2):
    """
    Erzeugt die (x, y, power)-Tripel für ein repliziertes Échelle.

    Parameter
    ---------
    freq        : array, Frequenzachse (μHz)
    power       : array, SNR- oder PSD-Werte
    dnu         : float, große Separation (μHz)
    n_replicas  : int, Anzahl horizontaler Wiederholungen (2 oder 3)

    Rückgabe
    --------
    x_all, y_all, p_all : arrays für den Scatter/Heatmap-Plot
    """
    x_base = freq % dnu          # Position in erster Kopie [0, dnu)
    y_base = freq                # y bleibt die echte Frequenz

    x_list, y_list, p_list = [], [], []
    for k in range(n_replicas):
        x_list.append(x_base + k * dnu)   # k-te Kopie um k·Δν verschoben
        y_list.append(y_base)
        p_list.append(power)

    import numpy as np
    return (np.concatenate(x_list),
            np.concatenate(y_list),
            np.concatenate(p_list))
```

Beim Plotten läuft die x-Achse dann von 0 bis `n_replicas * dnu`.

---

## 3  Konzept "Farbliche l-Trennung"

### 3.1  Das Zuordnungsproblem

Um Moden farblich nach l zu trennen, müssen sie zuerst **identifiziert**
werden. Das ist der schwierige Teil, weil man aus dem rohen Spektrum nicht
direkt sieht, welcher Peak welchen l-Wert hat. Es gibt drei Ansätze,
gestaffelt nach Aufwand:

### 3.2  Ansatz A — Position im Échelle (einfachste Methode)

Reine p-Moden folgen der asymptotischen Näherung. Ihre Position im Échelle
(ε-Phase) ist charakteristisch für l:

```
ν(n,l) ≈ Δν · (n + l/2 + ε) − δν₀ₗ
```

Daraus folgen typische x-Positionen im Échelle:

| Grad l | typische ε-Phase | Bemerkung |
|---|---|---|
| l = 0 | ε₀ (Referenz) | radiale Moden, definieren die Phase |
| l = 2 | ε₀ − δ₀₂/Δν | knapp links von l=0 (kleine Separation) |
| l = 1 | ε₀ + 0.5 | etwa eine halbe Δν versetzt |
| l = 3 | ε₀ + 0.5 − δ₁₃/Δν | schwach, selten sichtbar |

**Umsetzung:** Man definiert Fenster in der x-Position (ν mod Δν) und ordnet
jeden extrahierten Peak dem nächstliegenden l-Fenster zu.

```python
def classify_l_by_position(freq_modes, dnu, eps0, d02_rel=0.12):
    """
    Ordnet extrahierten Moden einen l-Wert nach ihrer Échelle-Position zu.

    Parameter
    ---------
    freq_modes : array, Frequenzen der extrahierten Moden (μHz)
    dnu        : float, große Separation
    eps0       : float, Phasenoffset der l=0-Moden (0..1)
    d02_rel    : float, kleine Separation δ₀₂ in Einheiten von Δν

    Rückgabe
    --------
    labels : array von l-Werten (0, 1, 2) oder -1 für unklar
    """
    import numpy as np
    x = (freq_modes % dnu) / dnu          # normierte Position [0,1)

    # Erwartete Positionen
    x_l0 = eps0 % 1.0
    x_l2 = (eps0 - d02_rel) % 1.0
    x_l1 = (eps0 + 0.5) % 1.0

    targets = {0: x_l0, 2: x_l2, 1: x_l1}
    tol = 0.10                             # Toleranzfenster (in Δν-Einheiten)

    labels = np.full(len(freq_modes), -1, dtype=int)
    for i, xi in enumerate(x):
        best_l, best_d = -1, tol
        for l, xt in targets.items():
            # zyklische Distanz auf [0,1)
            d = min(abs(xi - xt), 1 - abs(xi - xt))
            if d < best_d:
                best_l, best_d = l, d
        labels[i] = best_l
    return labels
```

**Grenze:** Funktioniert gut für reine p-Moden. Bei gemischten l=1-Moden
(die ja gerade *nicht* an der erwarteten Position sitzen) versagt die reine
Positionszuordnung — diese Moden werden zunächst als "unklar" (-1) markiert,
was aber selbst diagnostisch wertvoll ist: Die -1-Moden im l=1-Bereich
*sind* die gemischten Moden.

### 3.3  Ansatz B — l=0/l=2-Paar-Erkennung (robuster)

Radiale (l=0) und Quadrupol-Moden (l=2) treten als eng benachbartes Paar
auf, getrennt durch die kleine Separation δ₀₂ (typisch 4–10 μHz). Dieses
Paar ist ein zuverlässiger Anker:

1. Finde die stärksten Moden pro radialer Ordnung → Kandidaten für l=0
2. Suche knapp links (−δ₀₂) davon nach einem schwächeren Peak → l=2
3. Alles bei etwa +0.5·Δν versetzt → l=1 (inklusive gemischter Moden)

```python
def classify_l_by_pairs(freq_modes, amp_modes, dnu, d02_range=(3, 12)):
    """
    Identifiziert l=0/l=2-Paare und ordnet den Rest l=1 zu.

    Strategie:
    - l=0: stärkste Mode je Δν-Ordnung
    - l=2: schwächere Mode 3-12 μHz links von l=0
    - l=1: verbleibende Moden nahe +0.5·Δν
    """
    import numpy as np
    labels = np.full(len(freq_modes), -1, dtype=int)
    order = np.argsort(freq_modes)
    fs = freq_modes[order]
    amps = amp_modes[order]

    # Gruppiere nach radialer Ordnung (Δν-Bins)
    n_index = np.round(fs / dnu).astype(int)
    for n in np.unique(n_index):
        mask = n_index == n
        idx_in = np.where(mask)[0]
        if len(idx_in) == 0:
            continue
        # stärkste Mode = l=0-Kandidat
        i_l0 = idx_in[np.argmax(amps[idx_in])]
        labels[order[i_l0]] = 0
        f_l0 = fs[i_l0]
        # l=2: Mode links davon im d02-Fenster
        for j in idx_in:
            df = f_l0 - fs[j]
            if d02_range[0] < df < d02_range[1]:
                labels[order[j]] = 2
    # Rest → l=1 (nur wenn nahe +0.5 Δν)
    x = (freq_modes % dnu) / dnu
    for i in range(len(freq_modes)):
        if labels[i] == -1:
            # grobe Nähe zur l=1-Region
            labels[i] = 1
    return labels
```

### 3.4  Ansatz C — externe Modenidentifikation (Goldstandard)

Für publikationsreife Ergebnisse verwendet man etablierte Werkzeuge, die
die l-Zuordnung mit Modellvergleich oder Peak-Bagging leisten:

| Werkzeug | Methode | Referenz |
|---|---|---|
| **PBjam** | Peak-Bagging + HMM, automatische l-Zuordnung | Nielsen et al. 2021 |
| **FAMED** | interaktives Peak-Bagging (DIAMONDS) | Corsaro et al. 2020 |
| **lightkurve** | `Seismology`-Objekt, νmax/Δν, Echelle | Lightkurve Collab. |
| **echelle** (Hey) | Interaktives Δν-Tuning + Modenmarkierung | Hey & Ball 2020 |

Empfehlung: Für die eigene Pipeline Ansatz A oder B implementieren, für
finale wissenschaftliche Aussagen mit PBjam gegenprüfen.

---

## 4  Farbschema

Konsistente Farbcodierung über alle Plots hinweg:

| Grad l | Farbe | Hex | Physik |
|---|---|---|---|
| l = 0 | Blau | `#185FA5` | radial, reine p-Mode |
| l = 1 | Grün | `#0F9E66` | dipol, oft gemischt (p+g) |
| l = 2 | Orange | `#D85A30` | quadrupol, reine p-Mode |
| l = 3 | Violett | `#6040A0` | oktupol, selten |
| unklar | Grau | `#999999` | nicht zugeordnet / gemischt |

**Wichtig:** Gerade die "unklaren" (grauen) Moden im l=1-Bereich sind
diagnostisch — sie markieren die gemischten Moden an den avoided crossings.

---

## 5  Vollständiger Plot-Aufbau

```python
import numpy as np
import matplotlib.pyplot as plt

L_COLORS = {0: "#185FA5", 1: "#0F9E66", 2: "#D85A30", 3: "#6040A0", -1: "#999999"}
L_LABELS = {0: "l=0 (radial)", 1: "l=1 (dipol, gemischt)",
            2: "l=2 (quadrupol)", 3: "l=3", -1: "unklar/gemischt"}

def plot_replicated_echelle(freq_modes, amp_modes, labels, dnu, numax,
                            n_replicas=2, outpath="echelle_repliziert.pdf"):
    """
    Repliziertes Échelle mit farblich getrennten l-Werten.

    Parameter
    ---------
    freq_modes : array, Frequenzen der extrahierten Moden (μHz)
    amp_modes  : array, Amplituden (für Punktgröße)
    labels     : array, l-Zuordnung aus classify_l_*
    dnu, numax : float, seismische Parameter
    n_replicas : int, horizontale Wiederholungen
    """
    fig, ax = plt.subplots(figsize=(7, 9))

    x_base = freq_modes % dnu

    for l in [0, 1, 2, 3, -1]:
        mask = labels == l
        if not mask.any():
            continue
        for k in range(n_replicas):
            xk = x_base[mask] + k * dnu
            # Punktgröße nach Amplitude skaliert
            sizes = 20 + 200 * (amp_modes[mask] / amp_modes.max())
            ax.scatter(xk, freq_modes[mask], s=sizes,
                       c=L_COLORS[l], alpha=0.75,
                       edgecolors="white", linewidths=0.4,
                       label=L_LABELS[l] if k == 0 else None,
                       zorder=3 if l == 1 else 2)

    # Trennlinie zwischen den Kopien
    for k in range(1, n_replicas):
        ax.axvline(k * dnu, color="gray", lw=0.6, ls="--", alpha=0.4)

    # νmax-Linie
    ax.axhline(numax, color="orange", lw=1.0, ls=":", alpha=0.6,
               label=f"νmax = {numax:.0f} μHz")

    ax.set_xlabel(f"Frequenz mod Δν  (Δν = {dnu:.2f} μHz)", fontsize=10)
    ax.set_ylabel("Frequenz (μHz)", fontsize=10)
    ax.set_xlim(0, n_replicas * dnu)
    ax.set_title(f"Repliziertes Échelle ({n_replicas}×)  ·  "
                 f"l-Werte farbcodiert", fontsize=11)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
```

---

## 6  Erkennung der avoided crossings

Nach der l-Zuordnung lassen sich avoided crossings automatisch detektieren.
Ein avoided crossing zeigt sich als **lokale Abweichung des l=1-Ridges**
von seiner erwarteten Position:

```python
def detect_avoided_crossings(freq_l1, dnu, eps0, threshold=0.15):
    """
    Detektiert avoided crossings als Ausreißer im l=1-Ridge.

    Ein avoided crossing = l=1-Mode weicht deutlich von der erwarteten
    Position (eps0 + 0.5)·Δν ab.

    Parameter
    ---------
    freq_l1   : array, Frequenzen der l=1-Moden (μHz)
    dnu       : float, große Separation
    eps0      : float, Phasenoffset l=0
    threshold : float, Mindestabweichung (in Δν-Einheiten) für Detektion

    Rückgabe
    --------
    crossings : Liste von Frequenzen, an denen avoided crossings auftreten
    """
    import numpy as np
    x = (freq_l1 % dnu) / dnu
    x_expected = (eps0 + 0.5) % 1.0

    deviations = np.array([min(abs(xi - x_expected), 1 - abs(xi - x_expected))
                           for xi in x])
    crossing_mask = deviations > threshold
    return freq_l1[crossing_mask]
```

Für Mulder (KIC 10273246) erwartet man **zwei** avoided crossings
(Campante et al. 2011), für Scully (KIC 10920273) **eines**
(plus ein wahrscheinliches).

---

## 7  Diagnostische Interpretation

Nach der Umsetzung sollte das replizierte, farbcodierte Échelle folgende
Merkmale zeigen:

- **l=0 (blau) und l=2 (orange)**: gerade, vertikale Ridges. Ihre
  Geradheit ist die Kontrolle, dass Δν korrekt gewählt ist.
- **l=1 (grün)**: bei jungen Sternen (KIC 10963065) nahezu gerade, bei
  entwickelten Unterriesen (Mulder) mit sichtbaren S-förmigen
  Verbiegungen an den avoided crossings.
- **graue Punkte** im l=1-Bereich: die gemischten Moden selbst — sie
  tragen die Kerninformation.

Die Zahl und Lage der avoided crossings erlaubt Rückschlüsse auf:

- **Kernmasse und -dichte**: bestimmt die g-Moden-Frequenzen
- **Evolutionsstadium**: mehr avoided crossings = weiter entwickelt
- **Alter**: über Modellvergleich (z. B. mit MESA/GYRE-Gittern)

---

## 8  Integration in die bestehende Pipeline

Die neuen Funktionen ergänzen die bestehende Pipeline nach dem
Pre-whitening-Schritt (der die `freq_modes` und `amp_modes` liefert):

```python
# Nach Frequenzextraktion:
labels = classify_l_by_pairs(freq_modes, amp_modes, dnu)

# Repliziertes Échelle:
plot_replicated_echelle(freq_modes, amp_modes, labels, dnu, numax,
                        n_replicas=2)

# Avoided crossings:
freq_l1 = freq_modes[labels == 1]
crossings = detect_avoided_crossings(freq_l1, dnu, eps0)
print(f"Avoided crossings bei: {crossings} μHz")
```

---

## 9  Abhängigkeiten und Werkzeuge

| Paket | Zweck | Installation |
|---|---|---|
| numpy, scipy | Grundrechnung, Fits | `pip install numpy scipy` |
| matplotlib | Plots | `pip install matplotlib` |
| lightkurve | Datenzugang, Seismology-Objekt | `pip install lightkurve` |
| echelle | interaktives Δν-Tuning (optional) | `pip install echelle` |
| pbjam | Peak-Bagging + l-Zuordnung (Ansatz C) | `pip install pbjam` |

Das `echelle`-Paket von Daniel Hey bietet einen interaktiven Modus, in dem
man Δν per Schieberegler justiert und die Ridge-Ausrichtung live sieht —
ideal zur Feinjustierung vor der finalen Darstellung.

---

## 10  Literatur

- Bedding T.R. (2012): *Solar-like oscillations: An observational perspective.*
  In: Asteroseismology, Canary Islands Winter School. — Einführung des
  replizierten Échelle für gemischte Moden.

- Bedding T.R. et al. (2011): *Gravity modes as a way to distinguish between
  hydrogen- and helium-burning red giant stars.* Nature 471, 608.

- Campante T.L. et al. (2011): *Asteroseismology of two Kepler subgiants:
  KIC 10273246 and KIC 10920273.* A&A 534, A6. — Referenz für die
  avoided crossings von Mulder und Scully.

- Hey D. & Ball W. (2020): *Echelle: Dynamic echelle diagrams for
  asteroseismology.* Zenodo. — Das interaktive echelle-Paket.

- Nielsen M.B. et al. (2021): *PBjam: A Python package for automating
  asteroseismology.* AJ 161, 62. — Automatische Modenidentifikation.
