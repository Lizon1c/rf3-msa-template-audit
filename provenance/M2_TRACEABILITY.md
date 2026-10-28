# TRACEABILITY.md — M2 paper claim → source map

Paper: `paper.tex` (sidepapers/m2_state_selection/). Every quantitative claim below maps to
an archived data file or a ledger line. `AGENTS.md` = `/mnt/j/conda_envs/foundry/DMS_Project/AGENTS.md`
(line numbers as of 2026-09-07). `findings.md` =
`msa_role_properties/docs/findings.md`. All conformational numbers are **range-matched**
(`e2_rangematched_results.json`); pre-matching values are superseded (see Contradictions).

## E2 (13 targets, Fig 1, Table 1)

| Claim in paper | Value | Source |
|---|---|---|
| KaiB shift, sign flip | pref_no +0.72 → msa −1.57, shift −2.29 | `e2_rangematched_results.json` E2-01 (`pref_nomsa` 0.7201, `pref_msa` −1.5691, `shift` −2.2892); AGENTS.md l.288 |
| RfaH shift | +0.10 → −2.89, shift −2.99 | same file E2-02 (0.1020, −2.8861, −2.9881); AGENTS.md l.288 |
| Mad2 toward (crystal = state B) | −0.05 → −0.96, shift −0.91 | E2-03 (−0.0486, −0.9586, −0.9101); AGENTS.md l.288; findings.md §6 |
| XCL1 toward | −0.78 → +2.09, shift +2.87 | E2-04 (−0.7782, +2.0899, +2.8681); AGENTS.md l.288 |
| selecase null | +0.13 → +0.13, shift −0.01 | E2-05 (0.1340, 0.1283, −0.0057); findings.md §6 |
| GA88/GB88 in-flow sharpening | +0.52 → +3.07, shift +2.56 | E2-06 (0.5146, 3.0742, 2.5596); AGENTS.md l.288 |
| GA88/GB88 per-sequence (nomsa→msa) | GA88 +0.38→+2.74; GB88 −0.23→−2.66 | `E2-06/e49b_results.json`; AGENTS.md l.302 (`e2_06_ga_gb.py`) |
| β2AR/Abl1/EGFR/LacY/ubiq/GB1/TrxA nulls | +0.15, +0.01, −0.02, −0.08, −0.03, +0.09, +0.13 | E2-07..13 `shift` fields; findings.md §6 tally |
| LacY/EGFR pre-matching collapse | +0.94→−0.08; +0.31→−0.02 | findings.md §6 "Range-matched rerun"; AGENTS.md l.288 |
| Noise band ±0.14–0.27 | — | findings.md §6 caveats (archived repeat runs) |
| Shared-universe sizes (39–388 res) | Table 1 "res." column | `shared_res` fields, same JSON |
| MSA depths (rows incl. query) | 257/160/2/24/89 | `inputs/e2_targets/msas/E2-*.a3m` header counts |
| RfaH compression bug (101–114, 14-res shift) | — | findings.md §6; `e2_rangematched_analyze.py` docstring |

## E49 / E49b (template vs MSA, Figs 2–3, Table 2)

| Claim | Value | Source |
|---|---|---|
| KaiB tpl arms | tplA +3.59, tplB −3.00, msa_tplA +3.36, steer −6.59, conflict −0.23 | E2-01 rows of `e2_rangematched_results.json` |
| RfaH tpl arms | +1.49 / +0.30 / +1.45, steer −1.19, conflict −0.04 | E2-02 rows |
| Mad2 tpl arms (incl. tplA anomaly) | tplA −0.17 vs nomsa −0.05; tplB −1.87; conflict +0.001 | E2-03 rows; anomaly: AGENTS.md l.304(d), l.310 |
| XCL1 tpl arms | tplA −1.49, tplB −0.86, msa_tplA −1.44, steer +0.63 | E2-04 rows; findings.md §6 (anomaly resolved) |
| LacY tpl arms | +0.54 / +0.08 / +0.51, steer −0.47, conflict −0.03 | E2-10 rows |
| Conflict ≈ 0 on all 5 (\|Δ\| ≤ 0.23) | — | `conflict` fields; findings.md §6 |
| GA88/GB88 four arms + combined | GA88: +0.38/+2.74/+4.88/−4.30/−4.24; GB88: −0.23/−2.66/−4.34/+4.75/+4.58 | `e2/E2-06/e49b_results.json` (verbatim); AGENTS.md l.306 |
| Hierarchy template ±4–5 > MSA ±2.7 > sequence ±0.2–0.4 | — | AGENTS.md l.306 (same JSON) |
| AF2 hierarchy (MSA > template) | quote | Kovalevskiy et al. 2024 PNAS 121(34):e2315002121, PMC11348012 ("If the MSA signal is strong, AlphaFold tends to ignore structural information from the template") |
| E49 v1 template-dropped run quarantined | — | AGENTS.md l.304, l.533 (`/mnt/k/output_heads/e49_quarantine_templatedropped/`) |

## E49-RBD (MSA ≠ template; Fig 5, Table 3)

| Claim | Value | Source |
|---|---|---|
| eff90 181→139 (MSA) / 181→160 (tpl) / 159 (msatpl) | — | `rbd/e1_msa_arms_analysis/results_e49.json` arms `noMSA_ref`/`origMSA_ref`/`zii_e49_tpl`/`zii_e49_msatpl` |
| CKA vs no-MSA 0.568 / 0.443 / 0.448; CKA(tpl,MSA)=0.482 | — | same file (`cka_vs_noMSA`, `cka_vs_origMSA` of tpl arm) |
| Δ-direction median cos 0.48 (MSA) vs 0.26 (tpl) | — | `delta_cos_vs_noMSA_median` (0.4826 / 0.2576 / 0.2611) |
| Δ-norm ratio 1.21 (MSA) vs 1.51 (tpl) | — | `delta_norm_ratio_median` (1.2061 / 1.5143 / 1.5164) |
| msatpl ≡ tpl (Δ ≤ 0.01) | — | same file; AGENTS.md l.308(c) |
| Supervised: no 0.415, full15 0.638 (+0.223), tpl 0.366 (−0.049), msatpl 0.370 | — | `rbd/e1_supervised/results_e49.json` `summary` (0.41514/0.63794/0.36609/0.37029); AGENTS.md l.310 |
| Honest full-set MSA effect +0.14–0.15 | — | findings.md §1 (val2@same +0.141 / best_val2 +0.143 / strict +0.146) |
| 3-level falsification statement | — | AGENTS.md l.310 |

## E3 (homodimer interface amplification, Fig 4)

| Claim | Value | Source |
|---|---|---|
| Paired (n=392): AA 0.712/1.054, BB 0.712/1.053, AB 0.749/1.280, AB−AA +0.037, p=5.5e-66 | — | `casp3/e3_blocks/results.json` `summary` (medians: cos_AA 0.7124, normratio_AA 1.0536, cos_BB 0.7116, 1.0528, cos_AB 0.7491, 1.2798; Wilcoxon p 5.54e-66); AGENTS.md l.286 |
| Unpaired (n=100): AA 0.713/1.055, BB 0.712/1.054, AB 0.749/1.280, p=3.9e-18 | — | `casp3/e3_blocks/zii_msa_unpaired/results.json` `summary` (0.71258/1.05500, 0.71172/1.05439, 0.74889/1.27992; p 3.897e-18) |
| Per-mutant paired-vs-unpaired cos 0.969 | — | AGENTS.md l.286 |
| Unpaired protocol = chain-B homolog rows shuffled | — | AGENTS.md l.286 |
| Consistency with ColabFold/ESMFold unpaired complexes | — | AGENTS.md l.286; lit_review.md §3.1 |

## Methods / model

| Claim | Source |
|---|---|
| RF3 inference: 10 recycles, 50 diffusion steps, batch 2, checkpoint rf3_foundry_01_24_latest_remapped | `msa_role_properties/scripts/launchers/launch_e2_targets.py` `init_engine()` |
| Distogram bins: 64 over 2–22 Å + no-contact | `e2_rangematched_analyze.py` (`EDGES`), `rf3/metrics/distogram.py` |
| Experimental maps: Cα distance maps | `e2_analyze.py` docstring; `inputs/e2_targets/structs/*.json` |
| Template filtering (single chain, original chain id kept; verified nonzero templated atoms) | `scripts/e49_filter_templates.py` docstring; AGENTS.md l.533 |
| Template tokens appended after query; slicing by component length | AGENTS.md l.533 |
| RBD DMS dataset (3998 mut) | AGENTS.md §2; Starr et al. 2020 |
| CASP3 DMS (homodimer, 488 tokens) | AGENTS.md §2; Roychowdhury & Romero 2022 |
| UniRef30 MSAs, ≤256 homologs + query | `experiment_index.md` MSA-build row |

## Literature provenance

All 33 references: [1]–[25]-equivalents verified in `msa_role_properties/docs/lit_review.md`
(§6 citation map, local PDFs/PMC). Newly web-verified 2026-09-07 for this paper:
Kovalevskiy 2024 (PMC11348012, quote confirmed), Alexander 2007 (PNAS 104:11963),
Porter & Looger 2018 (PNAS 115:5968), Chang 2015 (Science 349:324), Burmann 2012
(Cell 150:291), Luo 2004 (NSMB 11:338), Tuinstra 2008 (PNAS 105:5057),
Starr 2020 (Cell 182:1295, PMC7310626), Roychowdhury & Romero 2022 (Cell Death
Discov., doi:10.1038/s41420-021-00799-0), Lan 2020 (Nature 581:215, 6M0J).

## Known source-data contradictions / caveats (for the record)

1. **E3 unpaired N**: AGENTS.md l.286 says "400 mutants, zii_msa_unpaired" but the archived
   `zii_msa_unpaired/results.json` has `n=100`. Paper reports n=100 (what the file contains).
2. **CASP3 DMS year**: AGENTS.md §2 says "Roychowdhury 2020"; the publication is 2022
   (Cell Death Discovery, doi:10.1038/s41420-021-00799-0). Paper cites 2022.
3. **E2 pre-matching values superseded**: XCL1 +4.81, RfaH −1.42, KaiB −2.62, Mad2 −2.26,
   LacY +0.94, EGFR +0.31 (old `e2_results.json`) are invalid (span/compression bugs);
   paper uses only range-matched values. fig7_e2_conformation.png in msa_role_properties
   still shows the old numbers — it was NOT used in this paper.
4. **GA88/GB88 two scoring conventions**: e2_06 flow gives GA88 +0.38/+2.74 (also used by
   E49b); the range-matched main flow gives +0.52/+3.07 (shift +2.56). Both are reported
   and labeled in the paper; they differ because the residue universes differ.
5. **E49 conflict bound**: AGENTS.md l.304 states "|conflict| ≤ 0.20" (pre-matching);
   range-matched KaiB conflict is −0.23. Paper states ≤ 0.23.
6. **E49-RBD supervised protocol**: 500-subset + best-epoch val (optimistically biased);
   paper pairs it with the de-biased full-set estimate (+0.14–0.15) for the MSA arm.
7. **RBD MSA gain (+0.223)** is within the 500-subset protocol only; not directly comparable
   to the canonical full-set number.

## E427 contiguous-block robustness (added 2026-10-27, end of §2.3 + Methods §4.4)

Source: `e427/results_A2.json` in this directory (canonical copy:
`/mnt/k/output_heads/e427_blocksplit_lightprobe/`); script
`e427/e427_blocksplit_lightprobe.py`. RBD 500-subset, all four arms, contiguous
5-fold position-block CV, SinglePredictorV2 + per-channel std (canonical head,
not the E1 head; absolute levels compared within-protocol only).

- MSA gain **+0.259** (5/5 folds: +0.121/+0.373/+0.321/+0.223/+0.258); arm means
  no 0.3828 / msa 0.6419 (≈ E1-protocol 0.638).
- Template deficit **−0.129** (4/5 folds negative: −0.250/−0.152/−0.101/−0.138;
  fifth fold −0.004 = inside noise floor).
- MSA+template vs template alone **+0.003** (null; per-fold +0.025/−0.042/+0.001/
  +0.014/+0.016); MSA+template vs no-MSA −0.126.
- Editorial (no new data): CKA expanded at first use (§2.3); "honest caveat" →
  "explicit caveat" (Mad2 paragraph).
