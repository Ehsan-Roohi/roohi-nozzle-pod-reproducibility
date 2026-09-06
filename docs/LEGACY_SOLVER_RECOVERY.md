# Legacy nozzle solver recovery — 2026-09-06

Three author-supplied archives were inspected against the fifteen retained
back-pressure exports. The result materially narrows the missing lineage, but it
does not yet identify the exact article executable.

## What is now established

The compact `nozzle.rar` package contains the modified Bird-family Fortran
program `DSMC2`, source header `Nozzle.FOR`, 2-D version 1.0 (August 2008), and
the files it reads directly: `common.txt`, `Property.txt` and `Inputdata.txt`.
The source writes `DSMC2.RES`, `HISTORY.PLT`, `Convergence.plt`,
`Contour_2_zone.PLT` and wall-property files.

The retained input defines:

| Quantity | Value |
|---|---:|
| Inlet temperature | 500 K |
| Inlet pressure | 100 kPa |
| Back pressure | 25 kPa |
| Main cells | 100 x 30 |
| Buffer cells | 30 x 40 |
| Nozzle length | 205 micrometres |
| Outlet half-height | 69 micrometres |
| Throat half-height | 15 micrometres |
| Inlet half-height | 34 micrometres |

The cell counts generate the exact `101 x 31` main and `31 x 41` buffer zone
dimensions in the recovered exports. The 25-kPa input also establishes that the
numbers in `P=<value>` denote **back pressure in kPa**, not pressure ratio.

The separate `Pb=7-Tw=Tin=300.rar` archive contains ten wall-heat-flux
directories, a repeated 327,680-byte Windows executable, matching Fortran source
copies, Tecplot layouts/figures, convergence histories and 88,915,116-byte
`DSMC2.RES` restart records. The executable has PE timestamp
`2018-11-09 20:22:44 UTC` and carries a Compaq Fortran runtime marker. This is
valuable run-chain evidence for the later heat-flux campaign, but its 7-kPa back
pressure and 300-K inlet temperature do not match the fifteen-case article sweep.

The 2026 GHS package is later model-development work and is not the article
producer.

## Exact relationship of the full-domain files

The fifteen original two-zone `.PLT` exports in this repository were compared
with the fifteen `P=<kPa>full.dat` files recovered on Unity and published in
FlowMLLab commit `3f6a40a`. For all 15 cases and all 4,402 points in the original
main and buffer zones:

- columns other than `Y` are exactly equal after the `SINGLE`/float32 Tecplot
  conversion;
- `Y_full = Y_original - 92 micrometres` within float32 rounding;
- the full files then add two Tecplot mirror zones;
- those mirror zones share variables 1 and 3--12 with the source zones and write
  only reflected `Y` values.

This proves the deterministic original-export-to-full-domain mapping even though
the exact Tecplot layout/macro that performed the operation is not yet present.
It also identifies a physical problem in that transformation: `QY`, `V` and
`Txy` are odd under reflection about the centerline, but the Tecplot
`VARSHARELIST` shares them unchanged. Their mirrored signs are therefore wrong.

The original centerline row is also not a boundary value forced to `V=0`. In the
recovered Fortran exporter, the `j=NCY+1` nodal row is populated from the nearest
cell-centered sampled moments. Across the fifteen retained exports, the maximum
absolute value on that row is 4.257--4.740 m/s. This supports a cell-to-node
export explanation; it is not evidence that the particle solver violated its
symmetry boundary condition.

## Why the exact producer is still open

The 2009 candidate's output routine labels its ten properties as density,
translational temperature, rotational temperature, temperature, three velocity
components, Mach number, pressure and Knudsen number. The article exports instead
contain density, `QX`, `QY`, temperature, `U`, `V`, `Txy`, Mach number, pressure
and Knudsen number. The exact producing source must therefore be a different
revision with heat-flux and shear-stress moment sampling/export.

To close the final gap, search for a Fortran file containing all of the literal
strings `QX`, `QY`, `Txy` and `Contour_2_zone.PLT`, together with the fifteen
matching input directories or restart records. Until that revision is found, the
recovered 2009 code is labeled as a solver-family/P=25 candidate rather than the
article producer.
