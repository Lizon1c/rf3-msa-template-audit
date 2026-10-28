# E2 — Conformational-State-Selection Target List (RF3)

**Experiment E2:** test whether the RF3 structure predictor shows a systematic bias toward
one experimentally observed conformational state when a single sequence can adopt two
well-characterized states. Each target below has TWO experimentally resolved states.
For every target we record which state dominates the PDB (the "crystallized" / over-
represented state) so we can ask whether RF3 recapitulates that state preferentially.

- **Generated:** 2026-07-25
- **N targets:** 13 (within the requested 10–14)
- **All PDB IDs verified** against RCSB (`data.rcsb.org/rest/v1/core/entry/{id}`) and
  PDBe SIFTS mappings on 2026-07-25 — see *Verification* at the bottom.
- MSA-depth labels are rough a-priori estimates; confirm with `jackhmmer`/`hhblits` at runtime.
- "No DMS found" claims are from a literature/web search (July 2026), not exhaustive.

## Coverage of requested categories

| Category | Targets |
|---|---|
| (a) Fold-switchers | E2-01 KaiB, E2-02 RfaH, E2-03 Mad2, E2-04 XCL1/lymphotactin, E2-05 selecase, E2-06 GA88/GB88 |
| (b) GPCR active/inactive | E2-07 β2AR |
| (c) Kinase DFG/αC states | E2-08 Abl1 (DFG-in/out), E2-09 EGFR (αC-in/out) |
| (d) Transporter inward/outward | E2-10 LacY |
| (e) Same-construct X-ray vs NMR | E2-11 ubiquitin, E2-12 GB1, E2-13 thioredoxin |

## Key references

- **AF-Cluster:** Wayment-Steele et al., *Nature* **625**:832–839 (2024) — MSA-clustering
  recovers multiple conformational states; their KaiB-TV validation structure is PDB 8UBH.
- **Fold-switching census:** Porter & Looger, *PNAS* **115**:5968 (2018);
  Dishman & Volkman, *Science* (2021, PMC8017559).
- **Selecase / AF2 failure:** Schafer et al., *bioRxiv* 617857 (2024) — AF2 cannot predict
  the selecase metamorphic switch.
- **Conformational diversity databases:** CoDNaS-Q (Saldaño et al., *Bioinformatics*
  **38**:4959, 2022); CoDNaS 2.0 (Monzon et al., *Database* baw038, 2016).

## Main table

| ID | Protein | Organism | UniProt | Category | State A (PDB / chain / res / method) | State B (PDB / chain / res / method) | Construct (residues) | Crystallized (dominant) state | MSA depth | Known DMS |
|----|---------|----------|---------|----------|--------------------------------------|--------------------------------------|----------------------|-------------------------------|-----------|-----------|
| E2-01 | KaiB | *Thermosynechococcus vestitus* BP-1 | Q79V61 | fold-switcher | **2QKE** / A–F / 2.70 Å / X-ray (ground state, GS) | **5JYT** / A / NMR / (fold-switched, fs); alt 5JWO chB 1.80 Å X-ray | 1–108 (WT) | GS (ground state) | shallow | none found |
| E2-02 | RfaH | *Escherichia coli* K-12 | P0AFW0 | fold-switcher | **2OUG** / A–D / 2.10 Å / X-ray (α-helical autoinhibited) | **2LCL** / A / NMR / (β-barrel CTD) | 1–162 | α-helical (autoinhibited) | moderate | none found |
| E2-03 | Mad2 (MAD2A) | *Homo sapiens* | Q13257 | fold-switcher | **1DUJ** / A / NMR / (O-Mad2; UniProt 11–195) | **1S2H** / A / NMR / (C-Mad2; 1–205) | 1–205 | C-Mad2 (1GO4 2.05 Å, 2VFX 1.95 Å) | moderate–deep | none found |
| E2-04 | XCL1 / lymphotactin | *Homo sapiens* | P47992 | fold-switcher | **1J8I** / A / NMR / (chemokine fold; mature 22–114) | **2JP1** / A,B / NMR / (β-sandwich); X-ray support 2NYZ chD/E 2.60 Å | 22–114 (mature) | chemokine fold | shallow | none found |
| E2-05 | Selecase (MJ1213) | *Methanocaldococcus jannaschii* | Q58610 | fold-switcher (metamorphic) | **4QHF** / A / 2.10 Å / X-ray (Slc-1 monomer) | **4QHJ** / A,B / 1.75 Å / X-ray (oligomeric autoinhibited; I100F/H107F mutant) | 1–110 | monomer | shallow | none found |
| E2-06 | GA88 / GB88 (designed) | designed (no organism) | N/A (designed) | fold-switcher (designed) | **2JWS** / A / NMR / (GA88, 3-α-helix, 56 aa) | **2JWU** / A / NMR / (GB88, α/β, 56 aa) | 56 aa (two sequences, 88% id) | — (two distinct sequences) | shallow | GB1 parent (Olson 2014) |
| E2-07 | β2-adrenergic receptor (β2AR) | *Homo sapiens* | P07550 | GPCR | **2RH1** / A / 2.40 Å / X-ray (inactive; T4L in ICL3) | **3SN6** / R / 3.20 Å / X-ray (active, Gs-bound); alt 3P0G 3.50 Å | 1–413 (native) | inactive | deep | **Jones et al., eLife 2020 (PMC7707821)** — ~7,800 singles, Gs/CRE signaling |
| E2-08 | Abl1 kinase | *Homo sapiens* | P00519 | kinase (DFG) | **2GQG** / A,B / 2.40 Å / X-ray (DFG-in active; +dasatinib) | **2HYY** / A–D / 2.40 Å / X-ray (DFG-out inactive; +imatinib); alt 1IEP (mouse) 2.10 Å | 228–500 (kinase domain) | DFG-in (active) | deep | none comprehensive; clinical resistance catalogs + Lyczek 2021 (PMC8609647) |
| E2-09 | EGFR kinase | *Homo sapiens* | P00533 | kinase (αC) | **1M17** / A / 2.60 Å / X-ray (αC-in active; +erlotinib) | **1XKK** / A / 2.40 Å / X-ray (αC-out inactive; +lapatinib) | 695–1022 (kinase domain) | active (αC-in) | deep | **Hayes 2024 Nat Commun (s41467-024-45594-4)** ~22,500 variants; + 2025 bioRxiv 646429 KD-DMS |
| E2-10 | LacY (lactose permease) | *Escherichia coli* K-12 | P02920 | transporter (MFS) | **1PV6** / A,B / 3.50 Å / X-ray (inward-facing; C154G) | **4OAA** / A,B / 3.50 Å / X-ray (outward-facing; G46W/G262W + sugar); alt 2V8N WT 3.60 Å, 5GXB 3.30 Å | 1–417 | inward-facing | deep (MFS) | none found (cf. Kaback classical scanning) |
| E2-11 | Ubiquitin | *Homo sapiens* | P0CG48 | X-ray vs NMR | **1UBQ** / A / 1.80 Å / X-ray | **1D3Z** / A / NMR | 1–76 | crystal (1UBQ) | very deep | **Roscoe 2013 JMB 425:1583 (PMC3615125)** all point mutants; + Roscoe & Bolon 2014 JMB 426:2854 |
| E2-12 | GB1 (protein G B1 domain) | *Streptococcus* group G | P06654 | X-ray vs NMR | **1PGB** / A / 1.92 Å / X-ray (auth 2–56) | **2GB1** / A / NMR; alt 3GB1 | 56-aa B1 domain (UniProt 228–282) | crystal (1PGB) | shallow | **Olson et al. 2014** singles+pairs, IgG-Fc binding; + Wu 2016 4-site combinatorial |
| E2-13 | Thioredoxin (TrxA) | *Escherichia coli* K-12 | P0AA25 | X-ray vs NMR | **2TRX** / A,B / 1.68 Å / X-ray (oxidized; UniProt 2–109) | **1XOA** / A / NMR (oxidized); 1XOB = reduced | 1–108 (auth) / 2–109 (UniProt) | crystal (2TRX) | deep | none found |

## Per-target notes & caveats

- **E2-01 KaiB.** 5JYT is a *fold-switch-stabilized variant* (≈6 substitutions vs WT, incl.
  Y8A and N30A — approximate, from manual alignment) and is missing residues 100–108. For RF3
  use the **WT 1–108** sequence; the fs state is only sparsely populated in WT. The AF-Cluster
  validation structure for KaiB-TV is 8UBH. Dominant PDB state = ground state (GS).
- **E2-02 RfaH.** The two states are the full-length α-helical autoinhibited form (2OUG) and
  the refolded β-barrel CTD (2LCL, residues 101–162). The β state also appears in cryo-EM
  elongation complexes (6C6S/6C6T). Construct = full length 1–162.
- **E2-03 Mad2.** Both reference states are NMR (O-Mad2 1DUJ, C-Mad2 1S2H). 1DUJ is an
  N-terminal truncation (UniProt 11–195); use full 1–205 for RF3. The C state dominates the
  PDB (crystal forms 1GO4, 2VFX and many checkpoint complexes).
- **E2-04 XCL1.** Both states solved by NMR; the β-sandwich (2JP1) also has X-ray support in
  the M3 chemokine-binding complex (2NYZ, chains D/E, 2.60 Å). Use the mature chemokine
  22–114 (auth 1–93).
- **E2-05 Selecase.** Concentration-dependent monomer↔oligomer metamorphic switch; the
  C-terminal helix refolds to a strand and swaps. The oligomeric structure 4QHJ is the
  **I100F/H107F mutant** (sequence-verified) that stabilizes the autoinhibited oligomer; the
  exact oligomeric order (dimer vs tetramer) was not pinned down here — treat as "oligomeric
  autoinhibited". AF2 is known to fail on this switch (Schafer 2024). Construct 1–110.
- **E2-06 GA88/GB88.** Designed proteins (He/Alexander/Bryan/Orban, *PNAS* 2008); 88%
  identical, 7 fold-determining substitutions, but they are **two different sequences**, so
  this is a designed fold-switch pair rather than a single-sequence switch. No UniProt
  accession. Parent GB1 has DMS data (Olson 2014). Use with the two-sequence caveat in mind.
- **E2-07 β2AR.** Both states are T4-lysozyme / fusion constructs. 2RH1 has T4L (P00720)
  replacing ICL3 (receptor ≈1–230 + 264–365); 3SN6 is the Gs-complex active structure
  (chain R, receptor 29–365 + T4L). For RF3 use the **native 1–413** sequence. Inactive state
  dominates the PDB. Strong DMS dataset (Jones eLife 2020).
- **E2-08 Abl1.** DFG-in active = 2GQG (+dasatinib); DFG-out inactive = 2HYY (+imatinib).
  Earlier candidate 1OPJ was rejected because it is *also* imatinib-bound (DFG-out), not a
  DFG-in partner. 1IEP (mouse Abl, P00520, 2.10 Å, 229–515) is an alternative DFG-out.
  Construct = human kinase domain 228–500. No comprehensive DMS; clinical resistance mutation
  catalogs and Lyczek 2021 (94 mutations) are the closest data.
- **E2-09 EGFR.** αC-in active = 1M17 (+erlotinib); αC-out inactive = 1XKK (+lapatinib), both
  on the same 695–1022 construct. (2GS6 is active + ATP-analog; 5C8M/5CAV are allosteric-site
  mutant structures, not used.) Active state dominates. Excellent full-length DMS (Hayes 2024).
- **E2-10 LacY.** Inward = 1PV6 (C154G mutant, sequence-verified; WT alternative 2V8N 3.60 Å);
  outward = 4OAA (G46W/G262W + sugar; no PDBe SIFTS — RCSB maps it to engineered UniProt
  B1XBJ1). Alternatives: 5GXB 3.30 Å (nanobody complex). Construct 1–417. Inward dominates.
  No DMS found (Kaback classical cysteine-scanning is the legacy data).
- **E2-11 Ubiquitin.** Cleanest crystal-vs-solution bias check at high resolution (1UBQ 1.80 Å
  vs 1D3Z NMR). Bias note: the C-terminal tail (72–76) is ordered in the crystal but flexible
  in solution. Very deep MSA; comprehensive DMS (Roscoe 2013/2014).
- **E2-12 GB1.** 56-residue B1 domain (UniProt P06654 residues 228–282; construct sequence
  MTYKLI…TVTE). 1PGB (1.92 Å) vs 2GB1 NMR. Cleanest single-domain bias control. Shallow MSA.
  DMS: Olson 2014 (singles at 55 positions + pairwise, IgG-Fc binding — exact journal likely
  *Curr. Biol.* 24:2641, flagged approximate) + Wu 2016 4-site combinatorial.
- **E2-13 Thioredoxin.** Both reference structures are the **oxidized** form (2TRX 1.68 Å vs
  1XOA NMR), so redox state is matched and the comparison isolates crystal-vs-solution bias
  (1XOB is the reduced NMR form — not used). Construct 1–108 (auth) = UniProt 2–109.

## Runners-up (considered, rejected — with reasons)

| Candidate | Reason rejected |
|---|---|
| **CcdB** (DNA gyrase poison) | 1VUB / 3VUB / 4VUB are three crystal forms of the *same* fold (Loris 1999), all dimers — no second conformational state in the PDB. |
| **MsbA** (ABC transporter) | All available structures are ≥4.5 Å (low resolution); listed as a transporter runner-up behind LacY. |
| **AdiC** (arginine/agmatine antiporter) | Only outward-facing structures exist (3L1L, 3LRB, 3LRC, 3OB6) — no inward-facing state to pair. |
| **Calbindin D9k** | 3ICB (2.3 Å, holo) vs 1CLB (NMR, apo) — the difference is confounded by Ca²⁺ ligation, not a pure conformational switch. |

## Verification

- **PDB existence / method / resolution:** RCSB REST
  `https://data.rcsb.org/rest/v1/core/entry/{id}` (re-checked 2026-07-25; all 32 IDs used
  above return valid entries with the methods/resolutions shown).
- **UniProt mappings & residue ranges:** PDBe SIFTS
  `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}` (fields `unp_start`/`unp_end`,
  segment lists per chain) and PDBe assembly summaries.
- **Caveats carried forward (honest uncertainties):** 5JYT mutation list is approximate;
  4QHJ oligomeric order undetermined; Olson 2014 GB1 journal citation is approximate;
  "no DMS found" for LacY/TrxA/Abl/fold-switchers is search-based (July 2026), not exhaustive;
  MSA depths are a-priori estimates to be confirmed at runtime.
