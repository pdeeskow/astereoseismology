"""Serieller, wiederaufnehmbarer Batchlauf für den seismischen Atlas."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from asteroseismologie import _atlas_output_dir, analyze_target
from atlas_config import load_targets
from atlas_models import ATLAS_SCHEMA_VERSION, TargetConfig


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seismischen Atlas seriell berechnen")
    parser.add_argument("--config", default="config/targets.yaml", help="Target-Konfiguration")
    parser.add_argument("--target", action="append", default=[], help="Nur diese Target-ID ausführen")
    parser.add_argument("--include-pending", action="store_true", help="Noch nicht analysierte Targets einschließen")
    parser.add_argument("--refresh-analysis", action="store_true", help="Vorhandene Atlas-Summary ignorieren")
    parser.add_argument("--refresh-pbjam", action="store_true", help="PBjam trotz Cache neu ausführen")
    parser.add_argument(
        "--reuse-existing-pbjam",
        action="store_true",
        help="Geprüften Alt-Cache trotz abweichender Signatur übernehmen",
    )
    parser.add_argument("--no-legacy-figures", action="store_true", help="Keine Einzelstern-Figuren erzeugen")
    return parser.parse_args()


def _cached_summary(target: TargetConfig) -> dict[str, Any] | None:
    path = _atlas_output_dir(target.id) / "summary.json"
    diagnostics = _atlas_output_dir(target.id) / "diagnostics.npz"
    if not path.exists() or not diagnostics.exists():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    expected = asdict(target)
    cached_target = data.get("target", {})
    for transient_key in ("pbjam_refresh", "pbjam_reuse_existing"):
        cached_target[transient_key] = expected[transient_key]
    if data.get("schema") != ATLAS_SCHEMA_VERSION or cached_target != expected:
        return None
    return data


def _literature_summary(target: TargetConfig) -> dict[str, Any]:
    reference = target.reference
    return {
        "schema": ATLAS_SCHEMA_VERSION,
        "target": asdict(target),
        "numax": reference.get("numax"),
        "numax_sigma": reference.get("numax_sigma"),
        "deltanu_auto": None,
        "deltanu_auto_sigma": None,
        "deltanu_used": reference.get("deltanu"),
        "deltanu_used_sigma": reference.get("deltanu_sigma"),
        "deltanu_source": "literature",
        "d02": reference.get("d02"),
        "d02_sigma": reference.get("d02_sigma"),
        "d02_pairs": 0,
        "stellar_parameters": {
            key: reference[key] for key in ("mass", "radius", "logg") if key in reference
        },
        "crossing_candidates": [],
        "mode_counts": {},
        "qc": [],
        "qc_status": "reference",
        "artifacts": {},
        "provenance": {"source": reference.get("source")},
    }


def _table_row(summary: dict[str, Any]) -> dict[str, Any]:
    target = summary["target"]
    stellar = summary.get("stellar_parameters", {})
    reference = target.get("reference", {})
    return {
        "id": target["id"],
        "label": target["label"],
        "stage": target["stage"],
        "numax": summary.get("numax"),
        "numax_sigma": summary.get("numax_sigma"),
        "deltanu_auto": summary.get("deltanu_auto"),
        "deltanu_used": summary.get("deltanu_used"),
        "deltanu_used_sigma": summary.get("deltanu_used_sigma"),
        "deltanu_source": summary.get("deltanu_source"),
        "d02": summary.get("d02"),
        "d02_sigma": summary.get("d02_sigma"),
        "mass": stellar.get("mass"),
        "radius": stellar.get("radius"),
        "logg": stellar.get("logg"),
        "crossing_candidates": len(summary.get("crossing_candidates", [])),
        "qc_status": summary.get("qc_status"),
        "reference_numax": reference.get("numax"),
        "reference_deltanu": reference.get("deltanu"),
        "reference_d02": reference.get("d02"),
        "reference_source": reference.get("source"),
    }


def _stored_summaries(root: str | Path = "results/atlas") -> list[dict[str, Any]]:
    return [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(Path(root).glob("*/summary.json"))
    ]


def write_results_tables(summaries: list[dict[str, Any]]) -> None:
    if not summaries:
        return
    rows = [_table_row(summary) for summary in summaries]
    output_dir = Path("results/tables")
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "results_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    headers = list(rows[0])
    markdown = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    markdown.extend("| " + " | ".join(str(row[key]) if row[key] is not None else "" for key in headers) + " |" for row in rows)
    (output_dir / "results_table.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")


def run_all(args: argparse.Namespace) -> int:
    requested = {" ".join(value.upper().split()) for value in args.target}
    summaries: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []
    for status, target in load_targets(args.config):
        normalized_id = " ".join(target.id.upper().split())
        if requested and normalized_id not in requested:
            continue
        if target.literature_only:
            print(f"[Literatur] {target.label}: keine Pipelineanalyse")
            summary = _literature_summary(target)
            summary_path = _atlas_output_dir(target.id) / "summary.json"
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
            summaries.append(summary)
            continue
        if status == "pending" and not args.include_pending and normalized_id not in requested:
            print(f"[Ausstehend] {target.label}: mit --include-pending einschließen")
            continue
        cached = None if args.refresh_analysis else _cached_summary(target)
        if cached is not None and not args.refresh_pbjam:
            print(f"[Cache] {target.label}")
            summaries.append(cached)
            continue
        try:
            run_target = replace(
                target,
                pbjam_refresh=args.refresh_pbjam,
                pbjam_reuse_existing=args.reuse_existing_pbjam,
            )
            result = analyze_target(
                run_target,
                write_atlas=True,
                legacy_outputs=not args.no_legacy_figures,
            )
            summaries.append(result.to_dict())
        except Exception as exc:
            failures.append((target.id, str(exc)))
            print(f"[Fehler] {target.label}: {exc}")

    write_results_tables(_stored_summaries())
    print(f"Atlaslauf: {len(summaries)} erfolgreich, {len(failures)} fehlgeschlagen")
    return 1 if failures else 0


def main() -> None:
    raise SystemExit(run_all(_parse_args()))


if __name__ == "__main__":
    main()