# analysis/ — RF3 MSA/template audit scripts

This directory holds the analysis and protocol scripts behind the RF3
MSA/template audit manuscript. Two tiers exist:

1. **Lightweight reproduction (portable, no GPU)** — figures, E428 structure
   scoring, and repository validation run entirely from the shipped
   `data/` tree.
2. **Heavy inference (not portable)** — scripts that run RF3 forward passes
   or model training still expect the original foundry environment; their
   external requirements are documented below.

## Python dependencies

The root `requirements.txt` lists the lightweight and common analysis
dependencies. The core packages below are sufficient for figure generation,
E428 scoring, and repository validation in a Python >= 3.10 environment:

| Package      | Used by                                        |
|--------------|------------------------------------------------|
| `numpy`      | all analysis scripts                           |
| `matplotlib` | `scripts/make_figures.py`                      |
| `Pillow`     | `scripts/make_figures.py` (TIFF post-processing)|
| `biotite`    | E428 analyze mode (mmCIF parsing)              |

`torch` is needed for model training and RF3 *run* modes, but not for figure
generation, repository validation, or E428 `analyze`. RF3 run modes
additionally need `hydra-core`, `omegaconf`, and the RF3 code stack itself
(see below).

## Lightweight reproduction

```bash
# Figures (writes figures/ by default; use --output-dir for a scratch dir)
python scripts/make_figures.py \
    --data-dir data --output-dir /tmp/rf3_fig_check

# E428 structure scoring: Needleman-Wunsch + Kabsch over the 480 archived
# CIFs; writes data/results/e428/results.json by default. GPU-free.
python analysis/protocols/e428/e428_structure_sampling.py analyze \
    --output /tmp/rf3_e428_check/results.json

# Repository integrity check (counts, no recompute)
python scripts/validate_repository.py
```

Canonical result files read by the lightweight path
(all under `data/results/`):

| Experiment | File                                  | Shape                |
|------------|---------------------------------------|----------------------|
| E2         | `e2/e2_rangematched_results.json`     | 13 entries           |
| E2 (E49b)  | `e2/e49b_results.json`                | GA88/GB88 template   |
| E428       | `e428/results.json`                   | 480 rows = 30x16     |
| E49-RBD    | `e49_rbd_supervised/results_e49.json` | 60 runs (4 arms x15) |
| E427       | `e427/results_A2.json`                | paired block-CV      |
| E429       | `e429/results_{A,B}.json`             | 2 protocols x 8 runs |

## Shared path configuration

`analysis/repo_paths.py` resolves all repo-relative paths from its own
location. Environment variables override the defaults:

| Variable           | Default                              | Purpose                        |
|--------------------|--------------------------------------|--------------------------------|
| `E428_MSA_DIR`     | `data/inputs/e2_targets/msas`        | E428 query MSAs                |
| `E428_TPL_DIR`     | `data/inputs/e49_templates`          | E428 state A/B reference CIFs  |
| `E428_STRUCT_DIR`  | `data/structures/e428`               | sampled-structure CIF root     |
| `E428_RESULTS_OUT` | `data/results/e428/results.json`     | analyze results destination    |
| `E428_RUN_OUT`     | = `E428_STRUCT_DIR`                  | run-mode CIF destination       |

## RF3 inference configuration (run modes only)

Run modes load the RF3 diffusion stack; nothing of it ships with this repo.
Configure it with:

| Variable         | Meaning                                                     |
|------------------|-------------------------------------------------------------|
| `RF3_ROOT`       | RF3/foundry model code root; expects `src/` and `models/rf3/` under it |
| `RF3_CHECKPOINT` | path to the remapped RF3 checkpoint (`.ckpt`)               |
| `RF3_CONFIG_DIR` | hydra configs (default: `$RF3_ROOT/models/rf3/configs`)      |
| `CUDA_VISIBLE_DEVICES` | GPU choice is the caller's; the script never pins an index |
| `E428_EXPECT_GPU` | optional card-model guard for run mode (e.g. `6000`); empty = skip |

## External feature-root requirements of the other protocols

The remaining scripts were copied from the larger project with their
original hardcoded paths left in place. They are **analysis provenance**,
not part of the lightweight path; to re-run them, provide the external
roots they reference (originally `/mnt/j/...` on the foundry box and
`/mnt/k/output_heads` for features):

- `analysis/core/e2_four_condition.py`, `analysis/core/fusion_v3.py` —
  need the project feature root (`output_heads/rbd` with per-condition
  `zii.pt`/`disto.pt`) and the project code root; `fusion_v3.py` also
  references the training checkpoint/log tree (`logs/ckpts`, `logs`).
- `analysis/protocols/e1/` — attention/embedding extraction for RBD arms.
  `attn_extract_unified.py` additionally references multiple model roots
  (RF3, FAESM, Evo2) and a worker script outside this repo;
  `launch_e1_msa_arms.py` / `launch_e7_recycle_gate.py` run RF3 inference
  (need `RF3_ROOT`/`RF3_CHECKPOINT` equivalents). `e1_analyze.py`,
  `e1_supervised.py`, `fix_e1_wt_replicates.py` read precomputed features
  from the RBD feature root.
- `analysis/protocols/e2/` — `launch_e2_targets.py`, `launch_e49_template.py`,
  `e2_06_ga_gb.py`, `e49b_ga_gb_template.py` run RF3 inference;
  `e49_filter_templates.py` writes into the project inputs tree;
  `e2_rangematched_analyze.py` / `e49_analyze_e2.py` read the E2 feature
  root (`output_heads/e2`).
- `analysis/protocols/e3/`, `analysis/protocols/canonical/` — read the
  casp3 / ci feature roots (`output_heads/casp3`, `output_heads/ci`) plus
  the project code root.
- `analysis/protocols/d1/` — RBD distillation suite; needs the RBD feature
  root and the project code root (GPU training).
- `analysis/protocols/e427/` — block-split light probe; imports project code
  and reads/writes the E427 feature root.
- `analysis/protocols/e429/` — paired noise floor; writes to its own
  output root (set via editing `OUT` at the top of the script).

Canonical JSON schemas are unchanged from the foundry originals. The
portable figure and E428 scripts read from the repository's `data/` tree;
the other provenance scripts retain their original external output roots.
The canonical result copies under `data/results/` are the inputs consumed by
`scripts/make_figures.py` and `scripts/validate_repository.py`.
