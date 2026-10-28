# Canonical Numerical Results

These JSON files are the analysis outputs cited by the unified manuscript. They include per-run values behind reported means, split differences, effective-rank/CKA summaries, and per-structure RMSD measurements. They are not regenerated from rounded manuscript tables.

## Directory Map

- `e1_geometry/`: canonical E1 eight-arm 500-mutant table, full-set RBD geometry, and E7 recycle-gating geometry.
- `e1_supervised/`: E1 dose/content and E7 supervised results.
- `e49_rbd_geometry/`, `e49_rbd_supervised/`: four-arm RBD no-MSA/MSA/template/joint comparisons.
- `rbd_canonical/`: RBD canonical four-condition and fusion matrices, including gradient diagnostics.
- `casp3/`, `ci/`: canonical MSA analyses and channel-standardized fusion matrices.
- `ab_validation/`: dual-validation, strict train-only standardization, input-noise, and learning-rate interventions.
- `d1/`: shallow/deep distillation ladders, single-modality and one-hot controls, CKA-objective variants, and functional recovery.
- `e427/`: contiguous-block cross-validation and lightweight-probe results, including all four RBD template arms.
- `e429/`: two independent processes with eight same-split/same-initialization training replicates each.
- `e2/`: range-matched 13-entry distogram results and the separate GA88/GB88 residue-universe scoring flow.
- `casp3_interface/`: paired and chain-B-shuffled unpaired per-mutant interface statistics.
- `e428/results.json`: 480 sample-level RMSD records used for all 30-cell coordinate summaries.

## Version Policy

The canonical tree intentionally excludes incomplete early E1 passes, position-0 artifacts, pre-range-matching scores, scalar-standardization backups, superseded pilots, debug files, and `*_smoke.json` stubs. The original claim-level provenance notes are in `provenance/`.

Where a result directory contains both a full-set and a 500-mutant analysis, the mutant universe is part of the protocol definition and must not be changed when comparing values. The two GA88/GB88 scoring flows are both retained intentionally; they use different residue universes and are not alternate versions of one number.

RBD metadata order is significant. Feature/label alignment must follow the metadata file order and must not be reconstructed by sorting mutant identifiers.
