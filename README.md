# RF3 MSA / Template Input-Channel Audit

Analysis code, evaluation protocols, per-run numerical outputs, structures, and figure generation for **An input-channel audit of RoseTTAFold3: MSA writing, mutation-effect prediction, and conformational-state selection**.

## Overview

AlphaFold3-style biomolecular structure predictors discard the explicit multiple sequence alignment (MSA) representation after a compact four-block module within each recycle, leaving the residue-pair representation ($\mathbf{Z}_{\mathrm{II}}$) to carry evolutionary information through the trunk. This repository contains the analysis code, evaluation protocols, per-run numerical outputs, and figure generation scripts for our multi-level audit of **RoseTTAFold3 (RF3)** across:

1. **Pair-representation geometry & one-shot MSA writing:** Effective rank (`eff90`), linear CKA, and mutant directional consistency (`dCos`) across SARS-CoV-2 RBD, human CASP3, and $\lambda$ repressor CI, plus $K$-recycle gating (Protocol E7) and alignment dose $\times$ content decomposition (Protocol E1).
2. **Deep mutational scanning (DMS) de-biasing & contiguous-block cross-validation:** Canonical transformer probes (`SinglePredictorV2`), dual-validation and strict train-only standardization de-biasing, 5-fold contiguous sequence-block cross-validation and lightweight MLP/linear probes (Protocol E427), three-stream fusion with ESM-2 650M (`FAESM`), paired training-noise floor audits (Protocol E429), and pair-representation distillation (Protocol D1).
3. **Conformational-state selection & explicit 3D coordinate sampling:** Range-matched distogram preference shifts across a 13-entry two-state panel (Protocol E2) and explicit 3D coordinate sampling across 30 condition cells ($N = 480$ structures; Protocol E428) spanning natural fold-switchers (KaiB, RfaH, Mad2, XCL1) and the designed GA88/GB88 pair.
4. **Template vs. MSA conflict resolution & readout decoupling:** Direct comparison of internal distogram negative log-likelihood (NLL) preferences against sampled 3D coordinates ($\Delta\mathrm{RMSD}$), alongside four-arm RBD template vs. MSA DMS benchmarks (Protocols E49-RBD and E427).

---

## Key Summary Findings

| Audit Layer | Experimental Protocol | Principal Finding |
| :--- | :--- | :--- |
| **MSA Write & Recycle Gating** | Three-protein geometry + Protocol E7 ($K \in \{1,2,3,10\}$) | MSA input compresses mutant ensemble effective rank (`eff90`) by 16–42% without erasing mutant directions; a single MSA recycle ($K=1$) installs ~74% of geometric compression and 85% of DMS functional gain, with nine subsequent MSA-free recycles erasing only ~25% of the geometric shift. |
| **Alignment Content vs. Rank** | Protocol E1 (9-arm RBD 500-mutant subset) + Protocol E427 | Uninformative alignments (`random`, `rowshuffle`) compress rank without functional gain. Relative to the single-row MSA pathway baseline (`depth01`), site-wise conservation contributes ~2/3 (`+0.095`) and covariation ~1/3 (`+0.042`) of the biological DMS gain, saturating near 8 rows. Contiguous 5-fold block-CV confirms robust net gains on RBD (`+0.256`) and CASP3 (`+0.264`). |
| **Conformational Selection** | 13-entry distogram panel + Protocol E428 ($N=480$ 3D structures) | Large MSA shifts concentrate in fold-switching proteins (2 away from crystal: KaiB, RfaH; 2 toward crystal: Mad2, XCL1; 1 own-fold sharpening: GA88/GB88; 8 within noise), consistent with family-consensus state selection and verified by explicit 3D coordinate sampling. |
| **Template Priors & Readout Decoupling** | Protocols E428, E49-RBD, and E427 block-CV | Wild-type template `6M0J` impairs RBD DMS prediction (`-0.049` subset; `-0.129` block-CV) and overrides the MSA gain when combined. In 3D coordinate sampling, templates dominate opposing MSAs in KaiB/RfaH/Mad2, whereas native MSAs preserve the native fold in GA88 (16/16) and yield a 50/50 split in GB88. Several template arms (`Mad2 tplA`, `XCL1 tplA`, `RfaH tplB`, `GA88/GB88 MSA+tplOther`) exhibit sign inversions or decoupling between distogram NLL and sampled 3D coordinates. |

## Repository Contents

- `analysis/`: protocol implementations and shared model/readout code.
- `scripts/`: figure generation and repository validation utilities.
- `data/inputs/`: MSA arms, two-state inputs, template mmCIFs, DMS labels, and metadata.
- `data/results/`: canonical per-run and per-mutant numerical outputs used by the manuscript.
- `data/structures/e428/`: all 480 unique sample-level E428 structures and 240 repetition metadata files.
- `manuscript/`: the unified manuscript, standalone S1 Text, Cover Letter, and S1--S4 Table workbooks.
- `figures/`: Fig 1, Fig 2, S1 Fig, and S2 Fig as vector PDF, 600-dpi PNG, and 600-dpi RGB LZW TIFF.
- `provenance/`: claim-to-source traceability ledgers from the two precursor studies.
- `zenodo/`: metadata and documentation for the external analysis-ready data package.
- `submission/`: the complete PLOS Computational Biology initial-submission directory, including upload files, source copies, metadata, and checksums.

Only sample-level `*sample-*_model.cif` files are included in E428. Ranked top-level duplicates are intentionally excluded to prevent double counting.

## Quick Start: Validate and Reproduce the Figures and Tables

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/validate_repository.py
python scripts/make_figures.py
python scripts/make_supporting_figures.py
python scripts/make_supporting_tables.py
python scripts/assemble_submission.py
python scripts/validate_submission.py
```

To recompute the E428 summary from the 480 deposited structures:

```bash
python analysis/protocols/e428/e428_structure_sampling.py analyze
```

The E428 analysis uses Needleman–Wunsch alignment and Kabsch superposition over the same shared residue subset for both references. It does not require a GPU.

## Full Analysis and External Dependencies

The lightweight reproduction path above does not rerun RF3 or ESM-2 inference. Full feature regeneration additionally requires:

- the upstream RF3/AtomWorks environment and a compatible checkpoint;
- ESM-2 650M (`esm2_t33_650M_UR50D`) weights;
- the per-mutant feature stores or the corresponding Zenodo analysis-ready package;
- protocol-specific external feature roots, as documented in `analysis/README.md`.

Model checkpoints, ESM weights, regenerable caches, and approximately 250 GB of unprocessed per-mutant feature stores are not committed to Git. See `zenodo/README.md` for the exact scope of the external package. Numerical results underlying every reported summary statistic are retained in `data/results/` even when full feature regeneration requires external inputs.

## Manuscript

- `manuscript/rf3_msa_template_audit_manuscript.tex` and `manuscript/M12.pdf`
- `manuscript/S1_text.tex` and `manuscript/S1_text.pdf`
- `manuscript/cover_letter.tex` and `manuscript/cover_letter.pdf`
- `manuscript/supporting_information/S1_Table.xlsx` through `S4_Table.xlsx`
- `submission/Manuscript.pdf`, `submission/Cover_Letter.pdf`, and the named upload files described in `submission/README.md`

The manuscript currently targets **PLOS Computational Biology**. The initial-submission package contains 11 upload files, LaTeX source copies, submission metadata, SHA-256 checksums, and a 13-check validation report. The analysis-ready external package is published at [10.5281/zenodo.23859152](https://doi.org/10.5281/zenodo.23859152). This GitHub repository and the Zenodo record are publicly available.

## Result Conventions

- Absolute Spearman $\rho$ values are compared only within the same feature universe, split scheme, normalization, readout head, and checkpoint-selection protocol.
- Distogram NLL and coordinate $\Delta\mathrm{RMSD}$ are different readouts. Their signs may be compared after orienting both toward the same state, but their units and uncertainties are not interchangeable.
- The $\pm0.02$ floor is reserved for unpaired independent-cell comparisons; paired interventions use split consistency and paired standard errors.
- CI did not undergo the complete DMS de-biasing program applied to RBD and CASP3.
- No universal template–MSA–sequence hierarchy is claimed. D1's unrecovered ~40% refers to MSA-specific functional advantage at the tested capacities, not the entire MSA write.

## License and Provenance

Original code and documentation are released under the MIT License; see `LICENSE`. Public input datasets and structures retain their original terms and citations. RF3, ESM-2, and third-party software remain subject to their respective licenses.

The original M1/M2 claim-level source ledgers are preserved in `provenance/`. Superseded position-0, pre-range-matching, pre-standardization-correction, and smoke-run results are excluded from the canonical result tree.

## Citation

Cheng Ding (2026). *An input-channel audit of RoseTTAFold3: MSA writing, mutation-effect prediction, and conformational-state selection*. Manuscript and accompanying analysis repository.
