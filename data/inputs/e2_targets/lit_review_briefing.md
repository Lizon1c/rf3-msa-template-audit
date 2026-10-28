# Briefing: MSA-bias literature review (for expansion into msa_bias_lit_review.md)

## Context (our project, 2 paragraphs max)
We study RF3's pair representation (Z_II) for DMS mutation-effect prediction, with and without MSA input. Clean corrected findings (2026-07-25, after discovering the June extraction was artifacted — all files sliced at position 0; the original claims of "rank collapse 231→30 and supervised −40%" were invalid and have been retracted internally):
- With MSA (15-row shallow GOLD_MSA): effective rank 181→139 (−23%), CKA(Z_msa, Z_no)=0.568, per-mutant delta direction preserved ~0.48; supervised cross-position Spearman IMPROVES 0.414→0.638 (+0.22; +0.15 over the honest depth-1 baseline).
- Dose × content arms (depth01/04/08/full15, colshuffle, consensus, rowshuffle, random): collapse has a generic shared-context component (even random sequences collapse rank to ~150) but the functional gain comes only from real biological content, mostly per-column conservation (colshuffle≈2/3 of the gain) with covariation secondary; saturates at ~8 sequences; one recycle of MSA installs ~75% of the geometric effect and ~85% of the functional gain (recycle-gated experiment E7); the trunk cannot self-correct afterwards.
- Framing we are developing: in AF3-style trunks the MSA is injected as a family-consensus mean (OuterProductMean over rows = the only MSA→Z channel in RF3), i.e. a one-shot "retrieval-key / family-consensus write", whose effect is to re-encode mutation signals toward family evolutionary statistics — which is exactly what cross-position DMS generalization needs.

## Papers to verify against LOCAL PDFs (extract exact quotes + section/page)
- /mnt/l/download/s41586-021-03819-2.pdf — AlphaFold2 (Jumper 2021). Key claims to verify: (a) "the MSA information is needed to coarsely find the correct structure within the early stages of the network, but refinement ... does not depend crucially on the MSA information" (MSA-depth section); (b) "AlphaFold is trained to produce the protein structure most likely to appear as part of a PDB structure"; (c) accuracy drop below ~30 MSA sequences; (d) MSA/template as separate input tracks; masked-MSA loss rationale.
- /mnt/l/download/s41586-024-07487-w.pdf — AlphaFold3 (Abramson 2024). Verify: "MSA processing is substantially de-emphasized" + 4-block MSA module with pair-weighted averaging; "the MSA representation is not retained and all information passes through the pair representation"; "they typically predict static structures as seen in the PDB, not the dynamical behaviour of biomolecular systems in solution"; E3-ligase "exclusively predicts the closed state"; interactions "without dependence on MSAs".
- /mnt/l/download/2022.07.20.500902.full.pdf — ESMFold/ESM-2 (Lin et al.). Verify: removal of MSA+templates; "language models are learning information similar to the contents of MSAs"; AF2/RoseTTAFold without MSA degrade below ESMFold; speed 60×; orphan/shallow-MSA proteins.
- /mnt/l/download/2024.07.01.600583v2.full.pdf — ESM3 (Hayes et al.). Verify: no MSA; EVMutation (pure MSA coevolution) beats ESM3 on viral DMS; AF2-MSA advantage concentrated on hard targets (CAMEO vs CASP14/15 table).
- /mnt/l/download/2025.08.14.670328v1.full.pdf — RF3/AtomWorks (Corley et al. 2025). Verify: implements AF3 algorithms; MSA only as featurized input; template is a separate track; NO MSA ablation in the paper.
- /mnt/l/download/41586_2023_Article_6415.pdf — RFdiffusion (Watson 2023). Verify: MSA never used; pretrained RF weights carry the prior.
- /mnt/l/download/s41586-018-0509-0.pdf — NOTE: this is Dou et al. 2018 (de novo β-barrel design), NOT AlphaFold1; flag as mislabeled; the real AF1 paper is Senior et al. 2020 (s41586-019-1923-7).

## Papers to cover via WEB (fetch abstracts/key statements, cite with URL)
- Del Alamo et al. 2022 eLife 75751 (MSA subsampling → alternative conformations).
- Stein & Mchaourab 2022 PLoS Comput Biol 18:e1010483 (SPEACH_AF).
- Wayment-Steele et al. 2024 Nature 625:832 (AF-Cluster).
- Kalakoti & Wallner 2025 (AFsample2, PMC11882827) and Wallner 2023 (AFsample).
- Chakravarty et al. 2024 (PMC11344769, fold-switch memorization) and Chakravarty & Porter 2022 (fold switching fails).
- Saldaño et al. 2022 Bioinformatics (conformational diversity vs AlphaFold).
- Sala et al. 2023 Frontiers (GPCR active/inactive bias quote).
- Heo & Feig 2022 Proteins (shallow MSA + templates control state).
- Mirdita et al. 2022 Nat Methods (ColabFold; paired MSA "surprisingly often not needed").
- Notin et al. 2023 NeurIPS (ProteinGym; Table A6 performance by MSA depth).
- MSA Pairformer 2025 bioRxiv 2025.08.02.668173 (contact-vs-variant-effect trade-off).
- The RSC Med Chem review d4md00869 quote: "AlphaFold2 has learnt a MSA-structure relationship, not a sequence-structure relationship as claimed" — try to identify their ref 71.

## Structure for the output doc (write in English, ~6-10 pages, with a citation map table)
1. Why MSA was introduced (coevolution lineage: DCA/EVcouplings → AF1 → AF2 Evoformer).
2. The de-emphasis/removal arc (ESMFold → ESM3 → AF3 → RF3) — with the exact quotes.
3. Evidence that deep MSA pins predictors to single/crystallized conformations (MSA-perturbation family, GPCR/kinase bias, memorization, "static structures as seen in the PDB").
4. The counter-evidence (MSA carries real evolutionary statistics: DCA, ProteinGym depth-dependence, EVMutation > ESM3 viral).
5. The gap: nobody has (a) quantified MSA's effect at the representation level (rank/CKA/direction), (b) decomposed it into shared-context vs conservation vs covariation channels, (c) measured its causal write-once dynamics (recycle gating), or (d) tested its net functional effect on mutation prediction with correct controls — our work fills this.
6. Citation map: table of claim → paper → exact quote → location.
