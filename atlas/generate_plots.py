"""Erzeugt alle möglichen Atlasplots ausschließlich aus Cache-Artefakten."""

from __future__ import annotations

from pathlib import Path

from atlas.plot_cd_diagram import plot_cd_diagram
from atlas.plot_common import load_summaries
from atlas.plot_echelle_grid import plot_echelle_grid
from atlas.plot_mixed_census import plot_mixed_census
from atlas.plot_scaling import plot_scaling


def main() -> None:
    summaries = load_summaries()
    if not summaries:
        raise SystemExit("Keine Atlas-Summaries unter results/atlas gefunden.")
    plots = [
        (plot_echelle_grid, "atlas_echelle.pdf"),
        (plot_cd_diagram, "atlas_cd_diagram.pdf"),
        (plot_scaling, "atlas_scaling.pdf"),
        (plot_mixed_census, "atlas_mixed_modes.pdf"),
    ]
    for plot, filename in plots:
        try:
            output = plot(summaries, Path("results/figures") / filename)
            print(f"Gespeichert: {output}")
        except ValueError as exc:
            print(f"Übersprungen: {filename} ({exc})")


if __name__ == "__main__":
    main()