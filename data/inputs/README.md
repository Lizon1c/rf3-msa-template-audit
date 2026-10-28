# Analysis Inputs

- `e1_msa_arms/`: query-only, real-depth, consensus, column-shuffled, row-shuffled, and random alignment arms for Protocol E1, with the arm manifest.
- `e2_targets/`: target definitions, MSAs, state distance maps, and preparation notes for the 13-entry two-state panel.
- `e49_templates/`: template mmCIF inputs used for the natural fold-switchers, GA88/GB88, RBD, and E428.
- `GOLD_MSA/`: the curated RBD 15-row, CASP3 69-row Swiss-Prot, and CI 512-row UniRef30 alignments used by the DMS analyses.
- `metadata/`: RBD and CASP3 mutation metadata and labels. Preserve file order when pairing features and labels.
- `dms/`: public source tables used for CASP3, ProteinGym lambda-repressor CI, and the orthogonal one-hot D1 control.

Public DMS datasets and PDB-derived templates retain their original licenses, terms, and citations. The MIT License in the repository root applies to original code and documentation, not to third-party datasets or model weights.

RF3 and ESM-2 checkpoints are external dependencies and are not included. The inputs here are sufficient to reproduce figures and summary statistics from `data/results/`; full feature regeneration additionally requires the external feature stores or the corresponding Zenodo package.
