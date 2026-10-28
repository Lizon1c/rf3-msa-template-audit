# Zenodo Analysis-Ready Package

This directory documents the external dataset package for the RF3 MSA/template audit. The upload itself is performed by the author; no Zenodo DOI exists yet.

## Intended Record

- **Title:** RF3 MSA/template input-channel audit: structures and analysis-ready features
- **Resource type:** Dataset
- **Creator:** Cheng Ding, School of Life Science and Technology, ShanghaiTech University
- **Related code:** https://github.com/Lizon1c/rf3-msa-template-audit

The package is generated without modifying any numerical result:

```bash
python scripts/prepare_zenodo.py \
  --feature-root /path/to/output_heads \
  --output /path/to/rf3_audit_zenodo_upload
```

The feature root must contain `e2/*/*/disto.pt` and `e427_blocksplit_lightprobe/cache/*.npy`.

## Package Contents

- `e428_sample_structures.zip`: all 480 unique sample-level E428 CIFs and 240 repetition metadata files, organized by target, arm, and repetition.
- `e427_feature_cache.zip`: the exact row-sliced NumPy feature caches used by E427 contiguous-block cross-validation and lightweight probes, including the RBD four template/MSA arms and the RBD/CASP3/CI no-MSA and MSA arms.
- `e2_distogram_outputs.zip`: all 53 per-arm distogram tensors underlying the range-matched two-state analysis.
- `canonical_results.zip`: the complete canonical JSON result tree used by the manuscript.
- `MANIFEST.json` and `SHA256SUMS.txt`: archive sizes, file counts, and integrity hashes.

The three larger archives make the key E427 and distogram analyses directly reproducible without downloading the full raw feature library. The 480 structures support independent recomputation of every E428 RMSD summary.

## Scope and Exclusions

This is an **analysis-ready subset**, not a copy of the approximately 250 GB unprocessed per-mutant feature library. It excludes RF3/ESM checkpoints, ranked duplicate E428 structures, smoke runs, superseded position-0 or pre-correction results, and regenerable intermediate files not needed by E427 or the distogram analysis.

Full regeneration of other feature stores (E1/E7 variants, D1 inputs, FAESM fusion inputs, and full CASP3 pair matrices) requires the upstream RF3/ESM environments and the protocol scripts in `analysis/`. Per-run numerical outputs for those protocols are nevertheless included both in the GitHub repository and in `canonical_results.zip`. Do not describe this package as the complete raw feature archive.

## Author Upload Checklist

1. Create a new Zenodo record using the title and creator information above.
2. Upload the four ZIP archives plus `MANIFEST.json` and `SHA256SUMS.txt`.
3. Add the GitHub repository as a related identifier and link the manuscript when a public DOI or URL is available.
4. Use CC BY 4.0 for original analysis artifacts where permitted; retain the original terms and citations for public DMS/PDB inputs and third-party software outputs.
5. Verify the uploaded file sizes and hashes against `SHA256SUMS.txt`.
6. Publish the record, then replace the Zenodo placeholder in `manuscript/rf3_msa_template_audit_manuscript.tex` with the real DOI and rebuild `M12.pdf`.

The repository must be public at submission. A private preparation repository and an unpublished Zenodo draft do not satisfy the manuscript's final Data Availability requirement.
