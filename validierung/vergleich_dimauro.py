"""
Vergleich der PBjam-Modenliste mit Di Mauro et al. 2011, Tabelle 3
===================================================================
KIC 4351319 ("Pooh")

Aufruf:
    python vergleich_dimauro.py <pbjam_modes.csv> [--tol 0.5] [--quality-min 0]

Erwartete Spalten in der PBjam-CSV (Namen werden heuristisch erkannt):
    Frequenz : freq | nu | nu_mean | frequency
    Grad     : l | ell | degree
    Ordnung  : n | n_p | order            (optional)
    Qualitaet: quality | prior_posterior_ratio | qual   (optional)
    Fehler   : freq_err | nu_std | sigma  (optional)

Ausgabe:
    - Zuordnungstabelle PBjam <-> Di Mauro
    - Trefferquote je Grad l
    - Liste der von PBjam verpassten Literaturmoden
    - Liste der PBjam-Moden ohne Literaturgegenstueck
    - Wenn Qualitaetsspalte vorhanden: Trefferquote als Funktion der Schwelle
"""

import sys
import csv
import argparse
from pathlib import Path

import numpy as np


REF_FILE = Path(__file__).with_name("dimauro_table3.csv")

COLMAP = {
    "freq":    ["freq", "nu", "nu_mean", "frequency", "freq_uhz", "nu_uhz"],
    "l":       ["l", "ell", "degree", "harmonic_degree"],
    "n":       ["n", "n_p", "order", "radial_order"],
    "quality": ["quality", "prior_posterior_ratio", "qual",
                "width_ratio", "prior_post_ratio"],
    "err":     ["freq_err", "nu_std", "sigma", "freq_sigma", "nu_err"],
}


def sniff_columns(fieldnames):
    """Ordnet die CSV-Spalten den benoetigten Groessen zu."""
    lower = {f.lower().strip(): f for f in fieldnames}
    found = {}
    for key, candidates in COLMAP.items():
        for c in candidates:
            if c in lower:
                found[key] = lower[c]
                break
    return found


def load_pbjam(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh)
        cols = sniff_columns(rdr.fieldnames or [])
        if "freq" not in cols:
            raise SystemExit(
                f"Keine Frequenzspalte erkannt. Vorhanden: {rdr.fieldnames}"
            )
        out = []
        for row in rdr:
            try:
                rec = {"freq": float(row[cols["freq"]])}
            except (TypeError, ValueError):
                continue
            rec["l"] = int(float(row[cols["l"]])) if cols.get("l") and row[cols["l"]] not in ("", None) else None
            rec["n"] = row.get(cols["n"], "") if cols.get("n") else ""
            for k in ("quality", "err"):
                if cols.get(k) and row[cols[k]] not in ("", None):
                    try:
                        rec[k] = float(row[cols[k]])
                    except ValueError:
                        rec[k] = None
                else:
                    rec[k] = None
            out.append(rec)
    return out, cols


def load_reference():
    ref = []
    with open(REF_FILE, newline="") as fh:
        for row in csv.DictReader(fh):
            if not row["l"]:
                continue          # marginale Moden ueberspringen
            ref.append({
                "freq": float(row["freq_uHz"]),
                "err":  float(row["freq_err"]),
                "snr":  float(row["snr"]),
                "l":    int(row["l"]),
                "conf": row["confirmed_all_teams"] == "1",
            })
    return ref


def match(pbjam, ref, tol):
    """Greedy-Zuordnung nach kleinstem Frequenzabstand."""
    pairs, used_ref = [], set()
    for p in sorted(pbjam, key=lambda x: x["freq"]):
        best, best_d = None, tol
        for i, r in enumerate(ref):
            if i in used_ref:
                continue
            d = abs(p["freq"] - r["freq"])
            if d < best_d:
                best, best_d = i, d
        if best is not None:
            used_ref.add(best)
            pairs.append((p, ref[best], best_d))
        else:
            pairs.append((p, None, None))
    missed = [r for i, r in enumerate(ref) if i not in used_ref]
    return pairs, missed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="PBjam-Modenliste")
    ap.add_argument("--tol", type=float, default=0.5,
                    help="Zuordnungstoleranz in uHz (Default 0.5)")
    ap.add_argument("--quality-min", type=float, default=0.0,
                    help="Nur Moden ueber dieser Qualitaetsschwelle")
    args = ap.parse_args()

    pbjam, cols = load_pbjam(args.csv)
    ref = load_reference()

    print("=" * 70)
    print("  KIC 4351319 - PBjam gegen Di Mauro et al. 2011, Tab. 3")
    print("=" * 70)
    print(f"  Erkannte Spalten: {cols}")
    print(f"  PBjam-Moden gesamt : {len(pbjam)}")
    print(f"  Literaturmoden     : {len(ref)}")

    if args.quality_min > 0:
        before = len(pbjam)
        pbjam = [p for p in pbjam
                 if p.get("quality") is not None and p["quality"] >= args.quality_min]
        print(f"  Nach Qualitaetsfilter >= {args.quality_min}: "
              f"{len(pbjam)} von {before}")

    pairs, missed = match(pbjam, ref, args.tol)

    # ---- Zuordnungstabelle -------------------------------------------
    print(f"\n  Zuordnung (Toleranz {args.tol} uHz)")
    print("  " + "-" * 66)
    print(f"  {'PBjam nu':>10} {'l':>2} {'Qual':>6} | "
          f"{'DiMauro':>9} {'l':>2} {'SNR':>5} | {'Delta':>7} {'l ok':>5}")
    print("  " + "-" * 66)
    l_agree = l_total = 0
    for p, r, d in pairs:
        q = f"{p['quality']:.2f}" if p.get("quality") is not None else "  -  "
        if r is None:
            print(f"  {p['freq']:>10.2f} {str(p['l']):>2} {q:>6} | "
                  f"{'---':>9} {'-':>2} {'-':>5} | {'-':>7} {'-':>5}")
        else:
            ok = "-"
            if p["l"] is not None:
                l_total += 1
                ok = "ja" if p["l"] == r["l"] else "NEIN"
                l_agree += (p["l"] == r["l"])
            print(f"  {p['freq']:>10.2f} {str(p['l']):>2} {q:>6} | "
                  f"{r['freq']:>9.2f} {r['l']:>2} {r['snr']:>5.1f} | "
                  f"{d:>+7.3f} {ok:>5}")

    # ---- Statistik ----------------------------------------------------
    matched = [x for x in pairs if x[1] is not None]
    print("\n  " + "-" * 66)
    print(f"  Zugeordnet          : {len(matched)} von {len(pbjam)} PBjam-Moden")
    print(f"  Literaturabdeckung  : {len(matched)} von {len(ref)} "
          f"({len(matched)/max(len(ref),1)*100:.0f} %)")
    if l_total:
        print(f"  l-Zuordnung korrekt : {l_agree} von {l_total} "
              f"({l_agree/l_total*100:.0f} %)")
    if matched:
        dev = np.array([x[2] for x in matched])
        print(f"  Frequenzabweichung  : Median {np.median(np.abs(dev)):.3f} uHz, "
              f"max {np.max(np.abs(dev)):.3f} uHz")

    # ---- pro Grad -----------------------------------------------------
    print("\n  Abdeckung je Grad")
    for l in (0, 1, 2):
        ref_l = [r for r in ref if r["l"] == l]
        hit_l = [x for x in matched if x[1]["l"] == l]
        print(f"    l={l}: {len(hit_l):>2} von {len(ref_l):>2} "
              f"({len(hit_l)/max(len(ref_l),1)*100:>3.0f} %)")

    # ---- verpasste Moden ----------------------------------------------
    if missed:
        print(f"\n  Von PBjam nicht gefunden ({len(missed)}):")
        print(f"    {'nu':>9} {'l':>2} {'SNR':>5}  alle Teams")
        for r in sorted(missed, key=lambda x: (x["l"], x["freq"])):
            print(f"    {r['freq']:>9.2f} {r['l']:>2} {r['snr']:>5.1f}"
                  f"  {'ja' if r['conf'] else 'nein'}")
        hi = [r for r in missed if r["snr"] >= 5.0]
        if hi:
            print(f"    davon {len(hi)} mit SNR >= 5 - diese sollten "
                  f"nachweisbar sein.")

    # ---- Qualitaetsschwelle durchfahren --------------------------------
    if any(p.get("quality") is not None for p in pbjam):
        print("\n  Abdeckung als Funktion der Qualitaetsschwelle")
        print(f"    {'Schwelle':>9} {'Moden':>6} {'Treffer':>8} {'Abdeckung':>10}")
        all_modes, _ = load_pbjam(args.csv)
        for thr in (0.0, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0):
            sub = [p for p in all_modes
                   if p.get("quality") is not None and p["quality"] >= thr]
            pr, _ = match(sub, ref, args.tol)
            hits = sum(1 for x in pr if x[1] is not None)
            print(f"    {thr:>9.2f} {len(sub):>6} {hits:>8} "
                  f"{hits/max(len(ref),1)*100:>9.0f} %")

    print("=" * 70)


if __name__ == "__main__":
    main()
