"""Recompute the density/velocity/pressure/Mach POD audit from all 15 snapshots.

The pipeline is loaded from ``src/`` and the 15 clean Tecplot files named
P=<pressure>.dat are read from ``data/``.
The physical and registered matrices use the same per-case jump normalization.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


PRESSURES = ["15", "16", "18", "19", "20", "22", "23", "24", "25",
             "26", "27", "28", "29", "30", "33"]
FIELDS = ["Rho", "U", "Pressure", "Mach"]
COORDINATES = ["x_physical", "xi_jump"]


def load_pipeline(script_path: Path):
    spec = importlib.util.spec_from_file_location("shock_centered_pipeline", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def first_mode_count(cumulative, threshold: float) -> int:
    for index, value in enumerate(cumulative, start=1):
        if value >= threshold:
            return index
    return len(cumulative)


def main() -> None:
    base = Path(__file__).resolve().parent
    repo_root = base.parent
    data_dir = repo_root / "data"
    pipeline_path = base / "5_Shock_Centered_ROM_MultiMetrics_PRF.py"
    if not pipeline_path.exists():
        raise FileNotFoundError(f"Missing recovered pipeline: {pipeline_path.name}")

    files = [data_dir / f"P={pressure}.dat" for pressure in PRESSURES]
    missing = [path.name for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing clean snapshots: " + ", ".join(missing))

    pipeline = load_pipeline(pipeline_path)
    cases = [pipeline.collect_case_profiles(path) for path in files]

    rows = []
    for field in FIELDS:
        for coordinate in COORDINATES:
            pod = pipeline.pod_metrics(cases, coordinate, field=field)
            energy = pod["energy_frac"]
            cumulative = pod["cumulative"]
            row = {
                "field": field,
                "coordinate": coordinate,
                "Ns": len(pod["labels"]),
                "E1_percent": 100.0 * float(energy[0]),
                "E1_to_E2_percent": 100.0 * float(cumulative[1]),
                "E1_to_E3_percent": 100.0 * float(cumulative[2]),
                "N99": first_mode_count(cumulative, 0.99),
            }
            rows.append(row)

    outdir = repo_root / "generated" / "Figures_6_9_Ns15"
    outdir.mkdir(exist_ok=True)
    csv_path = outdir / "Table_IV_Ns15_consistent_normalization.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    text_path = outdir / "Table_IV_Ns15_consistent_normalization.txt"
    lines = [
        "All-15-snapshot POD audit",
        "Same per-case jump normalization in physical and registered coordinates",
        "",
        "Field     Coordinate      E1 (%)   E1:2 (%)  E1:3 (%)  N99",
    ]
    for row in rows:
        lines.append(
            f"{row['field']:<9} {row['coordinate']:<15} "
            f"{row['E1_percent']:8.3f} {row['E1_to_E2_percent']:10.3f} "
            f"{row['E1_to_E3_percent']:10.3f} {row['N99']:4d}"
        )
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print(f"\nWrote {csv_path}")
    print(f"Wrote {text_path}")


if __name__ == "__main__":
    main()
