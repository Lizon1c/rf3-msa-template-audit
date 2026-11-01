# Zenodo Analysis-Ready Package

The external dataset package for the RF3 MSA/template audit is published at [10.5281/zenodo.23859152](https://doi.org/10.5281/zenodo.23859152).

## Record

- **Title:** RF3 MSA/template input-channel audit: structures and analysis-ready features
- **Resource type:** Dataset
- **Creator:** Cheng Ding, School of Life Science and Technology, ShanghaiTech University
- **Related code:** https://github.com/Lizon1c/rf3-msa-template-audit

The original package is generated without modifying any numerical result:

```bash
python scripts/prepare_zenodo.py \
  --feature-root /path/to/output_heads \
  --output /path/to/rf3_audit_zenodo_upload
```

The feature root must contain `e2/*/*/disto.pt` and `e427_blocksplit_lightprobe/cache/*.npy`.

## Published Split Package

The published record uses 21 independent ZIP archives, each no larger than approximately 90 MB, rather than the original four large archives:

- `e428_sample_structures.zip`: all 480 unique sample-level E428 CIFs and 240 repetition metadata files, organized by target, arm, and repetition.
- Fifteen `e427_*.zip` archives: the exact row-sliced NumPy feature caches used by E427 contiguous-block cross-validation and lightweight probes. Large arrays are losslessly sharded along the sample axis; complete arrays and shard locations are recorded in `E427_MANIFEST.json`.
- Four `e2_distograms_part*.zip` archives: all 53 original per-arm distogram tensors underlying the range-matched two-state analysis, with paths and hashes recorded in `E2_MANIFEST.json`.
- `canonical_results.zip`: the complete canonical JSON result tree used by the manuscript.
- `MANIFEST.json` and `SHA256SUMS.txt`: archive sizes, file counts, and integrity hashes for the split package.

To restore the original ten E427 arrays, keep the E427 archives, `E427_MANIFEST.json`, and `reassemble_e427.py` in one directory and run:

```bash
python reassemble_e427.py --output ./e427_feature_cache
```

The reassembled arrays were verified against every original array with exact equality. Extract all four E2 archives into the same destination to recover the original `e2/.../disto.pt` tree. The 480 E428 structures support independent recomputation of every E428 RMSD summary.

## Scope and Exclusions

This is an **analysis-ready subset**, not a copy of the approximately 250 GB unprocessed per-mutant feature library. It excludes RF3/ESM checkpoints, ranked duplicate E428 structures, smoke runs, superseded position-0 or pre-correction results, and regenerable intermediate files not needed by E427 or the distogram analysis.

Full regeneration of other feature stores (E1/E7 variants, D1 inputs, FAESM fusion inputs, and full CASP3 pair matrices) requires the upstream RF3/ESM environments and the protocol scripts in `analysis/`. Per-run numerical outputs for those protocols are nevertheless included both in the GitHub repository and in `canonical_results.zip`. Do not describe this package as the complete raw feature archive.

Original analysis artifacts use CC BY 4.0 where permitted; public DMS/PDB inputs and third-party software outputs retain their original terms and citations. The GitHub repository is publicly available alongside the Zenodo record.
