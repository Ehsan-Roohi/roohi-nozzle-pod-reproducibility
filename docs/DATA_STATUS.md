# Nozzle data status — original-export recovery, 2026-09-06 UTC

The fifteen historical paper snapshots remain available unchanged under `data/`.
The accompanying [original ASCII exports](../data/original_exports/) were recovered
from the corresponding author's existing research archive and added for provenance.
Their complete numerical correspondence to the historical files is checked by
`src/verify_original_exports.py`, not inferred from filenames or a plotted example.
This addition is not a new simulation campaign, a corrected dataset or a new release.

## Known limitation

The stated symmetry row contains nonzero exported transverse velocity in both
the recovered exports and the historical paper snapshots. Streamwise mass-flow
diagnostics also require examination, with explicit cell/node, boundary and
integration conventions. Their variation alone does not establish a solver
conservation failure. Numerical agreement between two exports verifies lineage;
it does not validate the underlying physical solution.

The producing solver/exporter revision, input deck, raw accumulated moments and
sampling record are still unresolved. An exporter defect is a hypothesis, not
a proven diagnosis. Recovery of a similar nozzle source is not sufficient to
attribute a specific defect to these runs.

## Before any corrected-data claim

1. Link the exact producing code, inputs and raw records to these fifteen cases.
2. Reproduce the original export and isolate the source of the discrepancy.
3. Correct at source and re-export from valid moments, or rerun if necessary.
4. Check symmetry, wall conditions, integrated flux and sampling/grid sensitivity;
   recompute affected paper figures and metrics with their original protocols.
5. Have the authors assess any publication correction; archive corrected data as
   a distinct version and preserve the historical version.

No corrected-data DOI or author/journal correction decision is implied here.
