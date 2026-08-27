# Shock-Centered Micro-Nozzle POD Reproducibility

Reproducibility package for manuscript Figures 6 and 9 and the associated POD audit in **Shock-Centered Low-Rank Structure and Shock-Aligned Surrogate Modeling of Rarefied Micro-Nozzle Flows** by Ehsan Roohi and Amirmehran Mahdavi.

The package resolves a provenance inconsistency in the earlier figures: the legacy plots used five pressure snapshots while later manuscript text referred to `N_s=15`. The primary workflow uses all 15 clean DSMC snapshots and applies the same per-case jump normalization in physical and shock-centered coordinates.

## Contents

- `src/regenerate_figures_6_and_9_ns15.py` — primary all-15-snapshot workflow for Figures 6 and 9.
- `src/audit_table_iv_ns15.py` — density, velocity, pressure, and Mach POD audit.
- `src/4_*`, `src/5_*`, and `src/7_*` — recovered legacy and intermediate analysis workflows retained for provenance.
- `data/` — 15 Tecplot ASCII pressure snapshots from 15 to 33 kPa.
- `reference_outputs/` — verified figures and machine-readable tables from the reference run.
- `docs/` — detailed Persian-language audit and proof-correction notes.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

Activate the virtual environment using the command appropriate for your operating system.

## Reproduce Figures 6 and 9

From the repository root, run:

```bash
python src/regenerate_figures_6_and_9_ns15.py
```

Generated files are written to `generated/Figures_6_9_Ns15/`. The workflow stops if any required snapshot is missing.

## Reproduce the Table IV audit

```bash
python src/audit_table_iv_ns15.py
```

The expected density results are:

| Coordinate | E1 (%) | E1:2 (%) | E1:3 (%) | N99 |
|---|---:|---:|---:|---:|
| Physical x | 78.922 | 88.467 | 93.762 | 8 |
| Shock-centered | 97.986 | 99.182 | 99.865 | 2 |

The first three individual shock-centered density-mode energies used in Figure 9 are 97.986%, 1.197%, and 0.683%.

## Reproducibility notes

- The primary analysis requires exactly the 15 snapshots listed by the manuscript.
- Figure 6 and Figure 9 are produced from the same POD calculation.
- The corrected shock-centered coordinate is `xi_j = (x - x_s) / delta_j`.
- The reference numerical outputs are included so a new run can be compared without relying on visual inspection alone.

## Authorship and licensing

Authorship and contribution information is recorded in [AUTHORS.md](AUTHORS.md). This public repository is a reproducibility record; no general open-source or open-data license is granted. See [LICENSE](LICENSE) and contact the corresponding author before redistribution or derivative use.

## Citation

Citation metadata is provided in [CITATION.cff](CITATION.cff). Please cite both the associated manuscript and the archived GitHub release.
