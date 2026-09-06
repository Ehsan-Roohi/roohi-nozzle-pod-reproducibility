# Recovered original ASCII exports — fifteen back-pressure cases

These files were recovered from Ehsan Roohi's pre-existing research archive
`NozzleMahdavi/BackPressure`. They are research simulation exports, not new
simulations generated for FlowMLLab and not AI-generated physical fields.
The filenames and file bytes, including original CRLF line endings, are preserved.

## Associated research

Ehsan Roohi and Amirmehran Mahdavi, micro-nozzle study:
[Physics of Fluids DOI 10.1063/5.0343101](https://doi.org/10.1063/5.0343101).
The earlier preprint is
[Shock-Centered Low-Rank Structure and Neural-Operator Representation of Rarefied Micro-Nozzle Flows](https://arxiv.org/abs/2605.12723).
The existing repository README uses the later manuscript title. Preprint and
later manuscript numerical summaries should not be assumed interchangeable.

The data retain the repository's [CC BY 4.0 license](../../DATA_LICENSE.md).
Credit both authors and cite the paper and the repository revision used.

## Relationship to the already-published files

- `BackPressure/Pout=<kPa>/<original-name>.PLT` contains the recovered export.
- `../P=<kPa>.dat` is the historical Tecplot re-save, already present at commit
  `e1b234ba499408d3b6224633972f939f3b2301d6`.
- [manifest.json](manifest.json) maps every file pair and records SHA-256 values.
- Each file contains a main `101 x 31` zone and a buffer `31 x 41` zone, with
  the same twelve columns: `X,Y,Density,QX,QY,T,U,V,Txy,Mach,Pressure,Kundsen`.
  Original labels are retained; they are not a complete unit specification.
- All **792,360 values** across all fifteen cases and both zones compare exactly
  after conversion to `float32`, the historical `.dat` files' `DT=SINGLE` precision.
  The ASCII text formats and their byte hashes are different. This is not a claim
  of identical text or exact equality between the original decimal strings.

Verify from the repository root (NumPy is the only dependency of this check):

```bash
python src/verify_original_exports.py
python -m unittest discover -s tests -v
```

The verifier reads retained evidence without modifying it. Its optional
`--output generated/export-lineage.json` writes a new report and refuses to
overwrite an existing file. The check also accepts CRLF working copies of the
historical `.dat` files by hashing their LF-normalized bytes; the recovered
`.PLT` files must match their original bytes exactly.

## Scientific limits — not corrected data

The nonzero exported transverse velocity on the stated symmetry row is present
in these earlier exports too. Thus the discrepancy did not originate in the
later `.dat` re-save or FlowMLLab import. This does **not** by itself identify
whether the cause is boundary treatment, sampling, cell-to-node export or indexing.
See the [data-status note](../../docs/DATA_STATUS.md).

These are macroscopic field exports, **not raw particle records or accumulated
DSMC moments**. The exact producing solver/exporter revision, matching input
deck, restart state and sampling metadata have not yet been linked to these
fifteen cases. The numerical parts of the filenames are preserved as archive
identifiers, not asserted to be verified timestep counts. Similar legacy nozzle
Fortran code is not substituted for the missing producing revision.

No field has been corrected, zeroed, interpolated or resimulated in this recovery.
The original `data/P=*.dat` files and retained reference outputs remain unchanged.
