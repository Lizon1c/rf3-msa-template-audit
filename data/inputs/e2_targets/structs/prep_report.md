# E2 structural reference prep report

Generated: 2026-07-26  
Source: `inputs/e2_targets/e2_targets.json` (13 targets x 2 states = 26 structures)  
Pipeline: RCSB download (legacy PDB via proxy) -> model-1 CA atoms of first specified chain, restricted to parsed `construct.boundaries` -> NxN CA distance map (float32, angstrom) saved as `{ID}_{stateA|stateB}_{PDB}.npy`; per-target metadata in `{ID}.json`.

**Totals: 26 structures OK, 0 failed.**

## Summary table

| Target | State | PDB | Chain (spec/used) | Boundaries parsed | Range used | N CA | Full chain | Status |
|---|---|---|---|---|---|---|---|---|
| E2-01 | stateA | 2QKE | A-F/A | 1-108 | 1-98 | 98 | 98 | OK |
| E2-01 | stateB | 5JYT | A/A | 1-108 | 1-106 | 106 | 106 | OK |
| E2-02 | stateA | 2OUG | A-D/A | 1-162 | 2-156 | 141 | 141 | OK |
| E2-02 | stateB | 2LCL | A/A | 1-162 | 97-162 | 66 | 66 | OK |
| E2-03 | stateA | 1DUJ | A/A | 1-205 | 11-197 | 187 | 187 | OK |
| E2-03 | stateB | 1S2H | A/A | 1-205 | 1-205 | 205 | 205 | OK |
| E2-04 | stateA | 1J8I | A/A | 22-114 | 22-93 | 72 | 93 | OK |
| E2-04 | stateB | 2JP1 | A,B/A | 22-114 | 22-60 | 39 | 60 | OK |
| E2-05 | stateA | 4QHF | A/A | 1-110 | 1-109 | 109 | 109 | OK |
| E2-05 | stateB | 4QHJ | A,B/A | 1-110 | 1-108 | 108 | 108 | OK |
| E2-06 | stateA | 2JWS | A/A | none (whole chain) | 1-56 | 56 | 56 | OK |
| E2-06 | stateB | 2JWU | A/A | none (whole chain) | 1-56 | 56 | 56 | OK |
| E2-07 | stateA | 2RH1 | A/A | 1-413 | 29-342 | 282 | 442 | OK |
| E2-07 | stateB | 3SN6 | R/R | 1-413 | 30-341 | 284 | 443 | OK |
| E2-08 | stateA | 2GQG | A,B/A | 228-500 | 228-500 | 272 | 277 | OK |
| E2-08 | stateB | 2HYY | A-D/A | 228-500 | 235-498 | 263 | 263 | OK |
| E2-09 | stateA | 1M17 | A/A | 695-1022 | 695-995 | 289 | 312 | OK |
| E2-09 | stateB | 1XKK | A/A | 695-1022 | 702-1018 | 289 | 289 | OK |
| E2-10 | stateA | 1PV6 | A,B/A | 1-417 | 1-417 | 417 | 417 | OK |
| E2-10 | stateB | 4OAA | A,B/A | 1-417 | 4-410 | 391 | 391 | OK |
| E2-11 | stateA | 1UBQ | A/A | 1-76 | 1-76 | 76 | 76 | OK |
| E2-11 | stateB | 1D3Z | A/A | 1-76 | 1-76 | 76 | 76 | OK |
| E2-12 | stateA | 1PGB | A/A | 228-282 | 1-56 | 56 | 56 | OK |
| E2-12 | stateB | 2GB1 | A/A | 228-282 | 1-56 | 56 | 56 | OK |
| E2-13 | stateA | 2TRX | A,B/A | 1-108 | 1-108 | 108 | 108 | OK |
| E2-13 | stateB | 1XOA | A/A | 1-108 | 1-108 | 108 | 108 | OK |

## Per-target details and flags

### E2-01 — KaiB (UniProt Q79V61)

Construct boundaries (raw): `1-108`

- **stateA** `2QKE` chain A (spec `A-F`): range used 1-98, 98 CA (full chain 98), saved `E2-01_stateA_2QKE.npy` (98x98).
  - flag: multi-chain spec 'A-F': used first chain 'A'
  - flag: partial coverage: range 1-108 expects 108 residues, extracted 98 (auth 1-98); chain auth span 1-98, full length 98
- **stateB** `5JYT` chain A (spec `A`): range used 1-106, 106 CA (full chain 106), saved `E2-01_stateB_5JYT.npy` (106x106).
  - flag: partial coverage: range 1-108 expects 108 residues, extracted 106 (auth 1-106); chain auth span 1-106, full length 106

### E2-02 — RfaH (UniProt P0AFW0)

Construct boundaries (raw): `1-162`

- **stateA** `2OUG` chain A (spec `A-D`): range used 2-156, 141 CA (full chain 141), saved `E2-02_stateA_2OUG.npy` (141x141).
  - flag: multi-chain spec 'A-D': used first chain 'A'
  - flag: partial coverage: range 1-162 expects 162 residues, extracted 141 (auth 2-156); chain auth span 2-156, full length 141
- **stateB** `2LCL` chain A (spec `A`): range used 97-162, 66 CA (full chain 66), saved `E2-02_stateB_2LCL.npy` (66x66).
  - flag: partial coverage: range 1-162 expects 162 residues, extracted 66 (auth 97-162); chain auth span 97-162, full length 66

### E2-03 — Mad2 (MAD2A) (UniProt Q13257)

Construct boundaries (raw): `1-205`

- **stateA** `1DUJ` chain A (spec `A`): range used 11-197, 187 CA (full chain 187), saved `E2-03_stateA_1DUJ.npy` (187x187).
  - flag: partial coverage: range 1-205 expects 205 residues, extracted 187 (auth 11-197); chain auth span 11-197, full length 187
- **stateB** `1S2H` chain A (spec `A`): range used 1-205, 205 CA (full chain 205), saved `E2-03_stateB_1S2H.npy` (205x205).

### E2-04 — XCL1 / lymphotactin (UniProt P47992)

Construct boundaries (raw): `22-114`

- **stateA** `1J8I` chain A (spec `A`): range used 22-93, 72 CA (full chain 93), saved `E2-04_stateA_1J8I.npy` (72x72).
  - flag: partial coverage: range 22-114 expects 93 residues, extracted 72 (auth 22-93); chain auth span 1-93, full length 93
- **stateB** `2JP1` chain A (spec `A,B`): range used 22-60, 39 CA (full chain 60), saved `E2-04_stateB_2JP1.npy` (39x39).
  - flag: multi-chain spec 'A,B': used first chain 'A'
  - flag: partial coverage: range 22-114 expects 93 residues, extracted 39 (auth 22-60); chain auth span 1-60, full length 60

### E2-05 — Selecase (MJ1213) (UniProt Q58610)

Construct boundaries (raw): `1-110`

- **stateA** `4QHF` chain A (spec `A`): range used 1-109, 109 CA (full chain 109), saved `E2-05_stateA_4QHF.npy` (109x109).
  - flag: partial coverage: range 1-110 expects 110 residues, extracted 109 (auth 1-109); chain auth span 1-109, full length 109
- **stateB** `4QHJ` chain A (spec `A,B`): range used 1-108, 108 CA (full chain 108), saved `E2-05_stateB_4QHJ.npy` (108x108).
  - flag: multi-chain spec 'A,B': used first chain 'A'
  - flag: partial coverage: range 1-110 expects 110 residues, extracted 108 (auth 1-108); chain auth span 1-108, full length 108

### E2-06 — GA88 / GB88 (designed) (UniProt None)

Construct boundaries (raw): `56 aa (two sequences, 88% identical)`

- **stateA** `2JWS` chain A (spec `A`): range used 1-56, 56 CA (full chain 56), saved `E2-06_stateA_2JWS.npy` (56x56).
  - flag: no numeric range in boundaries ('56 aa (two sequences, 88% identical)'): used whole chain
- **stateB** `2JWU` chain A (spec `A`): range used 1-56, 56 CA (full chain 56), saved `E2-06_stateB_2JWU.npy` (56x56).
  - flag: no numeric range in boundaries ('56 aa (two sequences, 88% identical)'): used whole chain

### E2-07 — beta2-adrenergic receptor (beta2AR) (UniProt P07550)

Construct boundaries (raw): `1-413`

- **stateA** `2RH1` chain A (spec `A`): range used 29-342, 282 CA (full chain 442), saved `E2-07_stateA_2RH1.npy` (282x282).
  - flag: partial coverage: range 1-413 expects 413 residues, extracted 282 (auth 29-342); chain auth span 29-342, full length 442
- **stateB** `3SN6` chain R (spec `R`): range used 30-341, 284 CA (full chain 443), saved `E2-07_stateB_3SN6.npy` (284x284).
  - flag: partial coverage: range 1-413 expects 413 residues, extracted 284 (auth 30-341); chain auth span 1002-341, full length 443

### E2-08 — Abl1 kinase (UniProt P00519)

Construct boundaries (raw): `228-500`

- **stateA** `2GQG` chain A (spec `A,B`): range used 228-500, 272 CA (full chain 277), saved `E2-08_stateA_2GQG.npy` (272x272).
  - flag: multi-chain spec 'A,B': used first chain 'A'
  - flag: partial coverage: range 228-500 expects 273 residues, extracted 272 (auth 228-500); chain auth span 223-500, full length 277
- **stateB** `2HYY` chain A (spec `A-D`): range used 235-498, 263 CA (full chain 263), saved `E2-08_stateB_2HYY.npy` (263x263).
  - flag: multi-chain spec 'A-D': used first chain 'A'
  - flag: partial coverage: range 228-500 expects 273 residues, extracted 263 (auth 235-498); chain auth span 235-498, full length 263

### E2-09 — EGFR kinase (UniProt P00533)

Construct boundaries (raw): `695-1022`

- **stateA** `1M17` chain A (spec `A`): range used 695-995, 289 CA (full chain 312), saved `E2-09_stateA_1M17.npy` (289x289).
  - flag: partial coverage: range 695-1022 expects 328 residues, extracted 289 (auth 695-995); chain auth span 672-995, full length 312
- **stateB** `1XKK` chain A (spec `A`): range used 702-1018, 289 CA (full chain 289), saved `E2-09_stateB_1XKK.npy` (289x289).
  - flag: partial coverage: range 695-1022 expects 328 residues, extracted 289 (auth 702-1018); chain auth span 702-1018, full length 289

### E2-10 — LacY (lactose permease) (UniProt P02920)

Construct boundaries (raw): `1-417`

- **stateA** `1PV6` chain A (spec `A,B`): range used 1-417, 417 CA (full chain 417), saved `E2-10_stateA_1PV6.npy` (417x417).
  - flag: multi-chain spec 'A,B': used first chain 'A'
- **stateB** `4OAA` chain A (spec `A,B`): range used 4-410, 391 CA (full chain 391), saved `E2-10_stateB_4OAA.npy` (391x391).
  - flag: multi-chain spec 'A,B': used first chain 'A'
  - flag: partial coverage: range 1-417 expects 417 residues, extracted 391 (auth 4-410); chain auth span 4-410, full length 391

### E2-11 — Ubiquitin (UniProt P0CG48)

Construct boundaries (raw): `1-76`

- **stateA** `1UBQ` chain A (spec `A`): range used 1-76, 76 CA (full chain 76), saved `E2-11_stateA_1UBQ.npy` (76x76).
- **stateB** `1D3Z` chain A (spec `A`): range used 1-76, 76 CA (full chain 76), saved `E2-11_stateB_1D3Z.npy` (76x76).

### E2-12 — GB1 (protein G B1 domain) (UniProt P06654)

Construct boundaries (raw): `56-aa B1 domain (UniProt 228-282; MTYKLI...TVTE)`

- **stateA** `1PGB` chain A (spec `A`): range used 1-56, 56 CA (full chain 56), saved `E2-12_stateA_1PGB.npy` (56x56).
  - flag: EMPTY intersection: range 228-282 vs chain auth 1-56 -> fell back to whole chain
- **stateB** `2GB1` chain A (spec `A`): range used 1-56, 56 CA (full chain 56), saved `E2-12_stateB_2GB1.npy` (56x56).
  - flag: EMPTY intersection: range 228-282 vs chain auth 1-56 -> fell back to whole chain

### E2-13 — Thioredoxin (TrxA) (UniProt P0AA25)

Construct boundaries (raw): `1-108 (auth) / 2-109 (UniProt)`

- **stateA** `2TRX` chain A (spec `A,B`): range used 1-108, 108 CA (full chain 108), saved `E2-13_stateA_2TRX.npy` (108x108).
  - flag: multi-chain spec 'A,B': used first chain 'A'
- **stateB** `1XOA` chain A (spec `A`): range used 1-108, 108 CA (full chain 108), saved `E2-13_stateB_1XOA.npy` (108x108).

## E2-06 GA88/GB88 sequences

Sequences extracted from the deposited PDB SEQRES records of 2JWS (GA88) and 2JWU (GB88); literature cross-check noted below.

- GA88 (2JWS): `TTYKLILNLKQAKEEAIKELVDAGIAEKYIKLIANAKTVEGVWTLKDEILTFTVTE` (56 aa)
- GB88 (2JWU): `TTYKLILNLKQAKEEAIKELVDAATAEKYFKLYANAKTVEGVWTYKDETKTFTVTE` (56 aa)

FASTA written to `/tmp/e2/seqs/E2-06_GA88.fa` and `/tmp/e2/seqs/E2-06_GB88.fa` (headers `>E2-06|GA88` / `>E2-06|GB88`).

**Literature verification (high confidence).** The sequences appear in figure form in the
primary literature — [Alexander et al. 2007, PNAS 104:11963-11968](https://pmc.ncbi.nlm.nih.gov/articles/PMC1906725/)
(design paper, Fig. 2) and [He et al. 2008, PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC2567172/)
(NMR structures of GA88/GB88, Fig. 1; PDB 2JWS/2JWU; BMRB 15535/15537). The SEQRES records
of 2JWS/2JWU are the authors' own machine-readable deposition of those same sequences, and
they satisfy every residue-level constraint stated in the paper texts: GA88 core/labelled
residues (K13, A16, L20, I33, V42, K46, I49, F52; disordered tails 1-8/54-56; PSD-1
mutations A6I, N7L, S8N, A52F, L53T) and GB88 (Y3, L5, L7, L9, A20, A26, F30, W43, Y45,
F52, V54; GB1 mutations G9L, L12A, G14E, G38T, T18K, A20L, V29Y, Q32L, N37K, D40E).
The two sequences differ at exactly 7 positions — 24 (G/A), 25 (I/T), 30 (I/F), 33 (I/Y),
45 (L/Y), 49 (I/T), 50 (L/K) — i.e. 49/56 = 88% identity, matching the published "88%"
and the seven fold-determining substitutions (positions match the GA95/GB95 design
discussion in He et al. 2008). The GA88 sequence also matches the one already recorded in
the `construct.notes` field of `e2_targets.json`.


## Construct-boundary vs PDB consistency review (task 5)

General rules applied: multi-chain specs (`A-F`, `A-D`, `A,B`) -> first chain used;
NMR ensembles -> model 1 only; altLoc A/blank only; filtering on auth residue number.
All 26 downloads succeeded; all 26 distance maps are square, symmetric, float32, with
zero diagonal. Per-state metadata (range used, chain length, flags, sequence) is in
`{ID}.json`.

### Major inconsistencies (affect cross-state comparability)

- **E2-04 XCL1 — numbering mismatch + missing density.** `boundaries` "22-114" is
  UniProt numbering; DBREF in both PDBs maps auth 1-93 = UniProt 22-114. Applied
  literally to auth numbers, the filter extracts auth 22-93 from 1J8I (72 CA) and
  22-60 from 2JP1 (39 CA). Additionally:
  - 2JP1 (beta-sandwich): residues 61-93 are declared missing in ALL 20 NMR models
    (REMARK 465) — only auth 1-60 has coordinates. The target notes already cite
    2NYZ chD/E (X-ray, 2.60 A) as full-coverage support for the beta state.
  - 1J8I (chemokine state): auth 1-93 is present, but the C-terminal ~64-93 is poorly
    defined in the ensemble (model 1: max CA-CA distance 81 A; C-tail residues ~60 A
    from the centroid). The saved map uses deposited model-1 coordinates as-is.
  - Recommendation: for a like-for-like map pair, restrict both states to auth 22-60,
    or substitute 2NYZ for stateB and use auth 1-93.
- **E2-07 beta2AR — fusion constructs, receptor-only extraction.** Both PDBs are T4L
  fusions. T4L residues are numbered >=1002 (2RH1: 1002-1161; 3SN6: 1002-1160), so the
  1-413 filter automatically excluded them: the npy files contain receptor CA only —
  2RH1 auth 29-342 (282 CA), 3SN6 chain R auth 30-341 (284 CA). ICL3 (replaced by T4L),
  termini, and flexible loops are absent; native 1-413 is not fully represented. The
  two extractions overlap on auth ~30-341.
- **E2-02 RfaH — stateB is CTD-only.** 2LCL covers auth 97-162 (66 CA, isolated CTD
  beta-barrel, as designed); 2OUG covers auth 2-156 (141 CA, disordered linker gap).
  The two maps span different residue ranges and are not index-comparable without
  alignment.
- **E2-03 Mad2 — stateA truncated.** 1DUJ covers auth 11-197 (187 CA, N-terminal
  truncation, as noted in the target JSON); 1S2H is complete 1-205.
- **E2-12 GB1 — non-numeric boundaries.** The first parseable range in "56-aa B1
  domain (UniProt 228-282; ...)" is 228-282 (UniProt), which has EMPTY intersection
  with auth 1-56 -> fell back to whole chain. This is exactly the intended 56-residue
  B1 domain (both PDBs cover auth 1-56, 56 CA each). Note: target notes say 1PGB
  "auth 2-56", but the downloaded legacy file is numbered 1-56.

### Minor inconsistencies (terminal/loop gaps)

- **E2-01 KaiB**: 2QKE chain A auth 1-98 (99-108 disordered); 5JYT auth 1-106
  (fold-switch-stabilized variant).
- **E2-05 Selecase**: 4QHF auth 1-109, 4QHJ auth 1-108 — the missing C-terminal
  residue(s) sit in the helix<->strand swap region that IS the switch element;
  interpret with care.
- **E2-08 Abl1**: 2GQG 228-500 nearly complete (272/273 CA, one internal gap);
  2HYY auth 235-498 (fraying termini missing).
- **E2-09 EGFR**: 1M17 auth 695-995 (996-1022 absent; chain also carries 672-694
  outside the construct range, excluded by the filter); 1XKK auth 702-1018.
- **E2-10 LacY**: 1PV6 complete 1-417; 4OAA auth 4-410 (391 CA; termini/loops missing).

### Clean targets

- **E2-06 GA88/GB88**: whole chain 1-56 used for both (no numeric range in
  boundaries; by design — two distinct 56-aa sequences).
- **E2-11 Ubiquitin**: 1-76 complete in both 1UBQ and 1D3Z.
- **E2-13 Thioredoxin**: 1-108 complete in both 2TRX and 1XOA (1XOA uses auth 1-108,
  so the auth interpretation of the boundaries applies cleanly).
