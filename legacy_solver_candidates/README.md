# Recovered legacy nozzle solver candidates

This directory preserves two small source packages supplied by the corresponding
author on 2026-09-06. They establish the solver family and several case inputs,
but neither package is claimed to be the exact executable/source revision that
produced all fifteen article snapshots.

## `nozzle_2009`

The source identifies itself as `Nozzle.FOR`, two-dimensional version 1.0
(August 2008), and builds program `DSMC2`. Its retained input is a nitrogen
micro-nozzle case with `PIN=100 kPa`, `POUT=25 kPa`, `FTMP=500 K`, a `100 x 30`
main cell grid and a `30 x 40` buffer grid. The resulting nodal Tecplot grid is
therefore `101 x 31` plus `31 x 41`, matching the two-zone dimensions of the
article exports. This is strong evidence for the solver family and the meaning
of `P=25`: **25 kPa back pressure**, not pressure ratio 25.

It is still a candidate rather than the exact producer. Its `OUT2` routine labels
the output moments as translational temperature, rotational temperature and
spanwise velocity, whereas the retained article exports contain `QX`, `QY` and
`Txy`. The exact source revision used for the article must contain those
additional moment calculations and labels.

## `ghs_2026`

This is a later generalized-hard-sphere (GHS) extension prepared in 2026. Its
input uses `POUT=7 kPa`, `FTMP=300 K` and `IMOLMODEL=2`. It documents later model
development and is not substituted for the article's earlier producer.

## Integrity and scope

[`manifest.json`](manifest.json) records the hashes of the supplied archives and
every retained file. The large `Pb=7-Tw=Tin=300.rar` heat-flux study is not copied
into Git: it is a separate 7-kPa/300-K campaign containing Tecplot layouts,
figures, executables, outputs and restart records. Its archive hash and the
relevant executable/source fingerprints are recorded in the manifest and in the
[recovery report](../docs/LEGACY_SOLVER_RECOVERY.md).

These legacy source artifacts are retained for scientific provenance. They are
not covered by the repository's root MIT license unless their original and
modified-code licensing is separately confirmed. The repository's data license
continues to apply only to the published data and reference outputs identified
there.
