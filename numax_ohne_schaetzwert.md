# νmax ohne Schätzwert: Selbstständige Bootstrapping-Pipeline

Beschreibung des vierstufigen Verfahrens zur automatischen Bestimmung von νmax
aus dem Powerspektrum — ohne jeden a-priori-Schätzwert.

---

## Hintergrund: Warum kein `argmax` auf der rohen PSD?

Das Powerspektrum setzt sich aus drei überlagerten Komponenten zusammen:

```
PSD(f) = Harvey-Granulation(f) + Weißrauschen + Oszillationsexzess(f)
```

| Komponente | Spektrale Form | Typische Stärke (Roter Riese) |
|---|---|---|
| Harvey-Granulation | `a / (1 + (f/b)^c)` | dominiert bei f < 50 μHz |
| Weißrauschen | flaches Plateau | dominiert bei f > 200 μHz |
| Oszillationsexzess | Gauß-Glocke um νmax | lokales Maximum, überlagert |

Das globale Maximum der rohen PSD liegt fast immer im Granulationsbereich
bei niedrigen Frequenzen — weit entfernt von νmax.

---

## Pipeline-Übersicht

```
PSD(f)
  │
  ├─[1]─ Harvey-Fit im Log-Raum ──────────────────► Hintergrundmodell B(f)
  │
  ├─[2]─ PSD / B(f) ──────────────────────────────► SNR-Spektrum
  │
  ├─[3]─ Breitbandige Glättung + argmax ──────────► erster Kandidat ν₀
  │
  └─[4]─ Gauß-Fit auf [0.55·ν₀, 1.55·ν₀] ────────► finales νmax
```

Jeder Schritt verwendet nur die Ausgabe des vorherigen — kein externer
Schätzwert wird benötigt.

---

## Schritt 1 — Harvey-Hintergrundmodell

### Physikalischer Hintergrund

Konvektive Granulation erzeugt ein charakteristisches `1/f`-ähnliches
Rauschen. Das Harvey-Modell beschreibt es als Lorentzian im Frequenzraum:

```
Harvey(f) = a / (1 + (f / b)^c)
```

- `a` — Amplitude (Stärke der Granulation)
- `b` — charakteristische Frequenz (Abknickpunkt der Granulation)
- `c` — Exponent (typisch 2.0–4.0 für Rote Riesen)

Für Rote Riesen liegt `b` typisch bei 10–50 μHz. Das Gesamthintergrundmodell
umfasst außerdem ein konstantes Weißrauschen `w`:

```
B(f) = Harvey(f) + w = a / (1 + (f/b)^c) + w
```

### Implementierung

```python
import numpy as np
from scipy.optimize import curve_fit

def harvey_background(f, a, b, c, w):
    """
    Vollständiges Harvey-Hintergrundmodell.

    Parameter
    ---------
    f : array, Frequenz in μHz
    a : float, Granulations-Amplitude (ppm²/μHz)
    b : float, charakteristische Frequenz (μHz)
    c : float, Exponent (dimensionslos)
    w : float, Weißrauschen-Plateau (ppm²/μHz)
    """
    return a / (1.0 + (f / b) ** c) + w

def fit_harvey_background(freq, power):
    """
    Fittet das Harvey-Modell im Log-Raum an das gesamte Powerspektrum.

    Der Fit im Log-Raum (log PSD vs. log f) ist robuster, weil er die
    logarithmische Streuung des χ²-verteilten Powerspektrums angleicht
    und den Dynamikbereich von mehreren Dekaden gleichmäßig gewichtet.

    Rückgabe
    --------
    bg_model : array, Hintergrundmodell B(f) für alle Frequenzen
    popt     : tuple, (a, b, c, w) — gefittete Parameter
    """
    log_power = np.log10(np.clip(power, 1e-10, None))
    log_freq  = np.log10(freq)

    def log_harvey(f, log_a, log_b, c, log_w):
        a = 10 ** log_a
        b = 10 ** log_b
        w = 10 ** log_w
        return np.log10(harvey_background(f, a, b, c, w))

    # Startwerte: grobe Schätzung aus Datenbereich
    p0 = [
        np.log10(power[0]),      # log_a: Amplitude ≈ niedrigstes PSD
        np.log10(freq[len(freq)//10]),  # log_b: Abknick bei ~10% der Bandbreite
        2.5,                     # c: typischer Exponent
        np.log10(np.percentile(power, 5)),  # log_w: Weißrauschen ≈ 5. Perzentile
    ]
    bounds = (
        [-2, -1, 0.5, -5],   # untere Schranken
        [ 8,  4, 6.0,  5],   # obere Schranken
    )

    try:
        popt_log, _ = curve_fit(
            log_harvey, freq, log_power,
            p0=p0, bounds=bounds, maxfev=10000
        )
        a  = 10 ** popt_log[0]
        b  = 10 ** popt_log[1]
        c  = popt_log[2]
        w  = 10 ** popt_log[3]
        bg_model = harvey_background(freq, a, b, c, w)
        return bg_model, (a, b, c, w)
    except RuntimeError:
        # Fallback: sehr breite Glättung als Proxy
        from scipy.ndimage import gaussian_filter1d
        sigma = max(1, len(freq) // 10)
        bg_model = gaussian_filter1d(power, sigma=sigma)
        return bg_model, None
```

> **Hinweis:** Für Sterne mit mehreren Granulationsskalen (z. B. Super-Riesen)
> verwendet man zwei Harvey-Terme. Für normale Rote Riesen genügt einer.

---

## Schritt 2 — SNR-Spektrum

Nach Abzug des Hintergrunds wird das Signal-zu-Rausch-Verhältnis gebildet:

```
SNR(f) = PSD(f) / B(f)
```

Im SNR-Spektrum gilt:
- Granulation → SNR ≈ 1 (flach)
- Weißrauschen → SNR ≈ 1 (flach)
- Oszillationsexzess → SNR > 1 (lokale Erhöhung)

Das globale Maximum im SNR-Spektrum liegt jetzt zuverlässig im
Modenbereich — ohne dass man weiß, wo dieser liegt.

```python
def compute_snr_spectrum(power, bg_model):
    """
    SNR-Spektrum: PSD normiert auf den gefitteten Hintergrund.

    Clipping verhindert Division durch sehr kleine Werte am Rand.
    """
    bg_safe = np.clip(bg_model, np.percentile(bg_model, 1), None)
    return power / bg_safe
```

---

## Schritt 3 — Erster νmax-Kandidat

Das SNR-Spektrum wird mit einer **frequenzunabhängigen** Breite geglättet.
Im Gegensatz zur rohen PSD muss die Glättungsbreite hier *nicht* von νmax
abhängen, weil das SNR-Spektrum kein steil abfallendes Kontinuum mehr hat.

Eine Glättungsbreite von 10–20 μHz funktioniert für alle solar-like
oscillators von Hauptreihe bis Riesenpunkt robust.

```python
from scipy.ndimage import gaussian_filter1d

def find_numax_candidate(freq, snr, smooth_width_uHz=15.0, fmin_search=5.0):
    """
    Erster νmax-Kandidat aus dem geglätteten SNR-Spektrum.

    Parameter
    ---------
    smooth_width_uHz : float
        Glättungsbreite in μHz. Frequenzunabhängig — kein Prior nötig.
        Empfehlung: 10–20 μHz für Rote Riesen.
    fmin_search : float
        Untergrenze des Suchbereichs in μHz. Schließt DC-Spike aus.

    Rückgabe
    --------
    numax_cand : float, erster Kandidat in μHz
    snr_smooth : array, geglättetes SNR-Spektrum
    """
    df = freq[1] - freq[0]
    sigma_bins = max(1, round(smooth_width_uHz / df))
    snr_smooth = gaussian_filter1d(snr, sigma=sigma_bins)

    # Suche nur oberhalb fmin_search
    search_mask = freq >= fmin_search
    idx_local   = np.argmax(snr_smooth[search_mask])
    numax_cand  = freq[search_mask][idx_local]

    return numax_cand, snr_smooth
```

---

## Schritt 4 — Gauß-Fit zum finalen νmax

Der erste Kandidat ν₀ definiert automatisch einen Bandpass:

```
Bandpass: [0.55 · ν₀,  1.55 · ν₀]
```

Die Faktoren 0.55 und 1.55 sind so gewählt, dass der Bereich asymmetrisch
ist (die Oszillationshülle ist nicht perfekt gaußförmig und fällt zur
niederfrequenten Seite steiler ab) und gleichzeitig weit genug, um den
gesamten Exzess einzuschließen. Im Bandpass wird ein Gauß + Konstante
gefittet:

```
Model(f) = A · exp(−(f − μ)² / (2σ²)) + c₀
```

μ ergibt das finale νmax, σ die Breite des Exzesses.

```python
from scipy.optimize import curve_fit

def refine_numax(freq, snr_smooth, numax_cand,
                 band_lo_factor=0.55, band_hi_factor=1.55):
    """
    Verfeinert den νmax-Kandidaten durch Gauß-Fit auf das geglättete
    SNR-Spektrum im automatisch bestimmten Bandpass.

    Parameter
    ---------
    band_lo_factor, band_hi_factor : float
        Bandpass-Grenzen relativ zu numax_cand.
        Asymmetrisch (0.55/1.55), weil die Oszillationshülle
        zur niederfrequenten Seite steiler abfällt.

    Rückgabe
    --------
    numax   : float, finales νmax in μHz
    sigma   : float, Breite der Gauß-Hülle in μHz
    """
    band_lo = numax_cand * band_lo_factor
    band_hi = numax_cand * band_hi_factor
    mask    = (freq >= band_lo) & (freq <= band_hi)

    f_fit = freq[mask]
    s_fit = snr_smooth[mask]

    def model(f, A, mu, sigma, c0):
        return A * np.exp(-0.5 * ((f - mu) / sigma) ** 2) + c0

    A0     = s_fit.max() - s_fit.min()
    sigma0 = numax_cand * 0.25
    p0     = [A0, numax_cand, sigma0, s_fit.min()]
    bounds = (
        [0,    band_lo, 1.0,  0   ],
        [1e6,  band_hi, band_hi-band_lo, np.inf],
    )

    try:
        popt, pcov = curve_fit(model, f_fit, s_fit, p0=p0,
                               bounds=bounds, maxfev=5000)
        numax_fit  = float(popt[1])
        sigma_fit  = float(abs(popt[2]))
        sigma_err  = float(np.sqrt(pcov[1, 1])) if pcov is not None else np.nan
    except RuntimeError:
        # Fit nicht konvergiert: Kandidat als Fallback
        numax_fit = numax_cand
        sigma_fit = numax_cand * 0.25
        sigma_err = np.nan

    return numax_fit, sigma_fit
```

---

## Vollständige Pipeline — Aufruf

```python
def estimate_numax_auto(freq, power, verbose=True):
    """
    Vollautomatische νmax-Bestimmung ohne Schätzwert.

    Parameter
    ---------
    freq  : array, Frequenz in μHz (gleichmäßig abgetastet)
    power : array, Powerspektrum in ppm²/μHz

    Rückgabe
    --------
    dict mit:
        numax       : finales νmax (μHz)
        numax_sigma : Breite der Gauß-Hülle (μHz)
        numax_cand  : erster Kandidat vor Gauß-Fit (μHz)
        harvey_bg   : Hintergrundmodell B(f) (array)
        snr         : SNR-Spektrum (array)
        snr_smooth  : geglättetes SNR-Spektrum (array)
    """
    # Schritt 1: Hintergrundmodell
    bg_model, harvey_params = fit_harvey_background(freq, power)
    if verbose and harvey_params:
        a, b, c, w = harvey_params
        print(f"  Harvey-Fit: a={a:.0f} ppm²/μHz,  b={b:.1f} μHz,  "
              f"c={c:.2f},  w={w:.1f} ppm²/μHz")

    # Schritt 2: SNR-Spektrum
    snr = compute_snr_spectrum(power, bg_model)

    # Schritt 3: Kandidat
    numax_cand, snr_smooth = find_numax_candidate(freq, snr)
    if verbose:
        print(f"  νmax-Kandidat (SNR-Maximum): {numax_cand:.1f} μHz")

    # Schritt 4: Gauß-Fit
    numax, numax_sigma = refine_numax(freq, snr_smooth, numax_cand)
    if verbose:
        print(f"  νmax final (Gauß-Fit):       {numax:.1f} μHz  "
              f"(σ = {numax_sigma:.1f} μHz)")

    return {
        "numax"       : numax,
        "numax_sigma" : numax_sigma,
        "numax_cand"  : numax_cand,
        "harvey_bg"   : bg_model,
        "snr"         : snr,
        "snr_smooth"  : snr_smooth,
        "harvey_params": harvey_params,
    }
```

### Verwendung im Hauptskript

```python
# Ersetzt den bisherigen Aufruf mit Schätzwert:
#   numax, sigma, smooth = fit_numax(freq, power, numax_guess=LIT_NUMAX)

result = estimate_numax_auto(freq, power)
numax       = result["numax"]
numax_sigma = result["numax_sigma"]
smooth      = result["snr_smooth"]   # für Échelle-Plots weiterverwendbar
```

---

## Robustheit und Grenzen

### Wann die Pipeline zuverlässig ist

- Solar-like oscillators mit SNR > 3 im Modenbereich
- TESS 2-min-Daten, mindestens ein Sektor (≥ 27 Tage)
- Rote Riesen: νmax zwischen ~5 und ~300 μHz

### Bekannte Schwachstellen

| Situation | Problem | Abhilfe |
|---|---|---|
| Sehr schwaches Signal (SNR < 2) | Schritt 3 findet Rauschmaxima | Mehr Sektoren laden oder `fmin_search` anpassen |
| Zwei Sterne im Pixel | Zwei Exzess-Glocken | Zwei-Gauß-Fit in Schritt 4 |
| δ Scuti-Sterne | Kein glattes Harvey-Kontinuum | Nur Schritt 3/4, andere Glättungsbreite (~0.5 μHz) |
| Sehr kurze Zeitreihe (< 14 d) | Frequenzauflösung zu grob für Δν < 10 μHz | Mindestens 2 Sektoren |

### Vergleich mit etablierten Pipelines

| Pipeline | Ansatz | Öffentlich |
|---|---|---|
| **pySYD** | Harvey-Fit + SNR + iterativer Gauß-Fit | ja (`pip install pysyd`) |
| **A2Z** | Glättung + Autokorrelation | ja (IDL/Python) |
| **CAN** | Bayesianischer Harvey-Fit | teilweise |
| **DIAMONDS** | MCMC auf vollem Modell | ja (`pip install diamonds`) |

Die hier beschriebene Pipeline entspricht im Wesentlichen dem Kern von
`pySYD` — als eigenständige, leicht anpassbare Implementierung.

---

## Erweiterung: Unsicherheit auf νmax

`scipy.optimize.curve_fit` liefert die Kovarianzmatrix `pcov`. Die
formale 1σ-Unsicherheit auf νmax:

```python
popt, pcov = curve_fit(...)
numax_err = np.sqrt(pcov[1, 1])   # Index 1 entspricht μ im Gauß-Modell
```

Für eine realistischere Unsicherheitsabschätzung empfiehlt sich Bootstrap-
Resampling über Frequenz-Subsets oder über mehrere TESS-Sektoren.
