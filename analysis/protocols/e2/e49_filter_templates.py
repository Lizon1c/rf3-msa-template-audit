#!/usr/bin/env python3
"""Filter multi-chain mmCIF templates to a single chain and keep its original chain id.

Keeps all non-atom records; ATOM/HETATM/ANISOU lines are kept only for the
requested auth chain, with both label_asym_id and auth_asym_id rewritten to 'Z'
(re-chaining breaks atomworks assembly building — oper lists reference original asym ids and the chain gets silently dropped).
"""
import re, sys
from pathlib import Path

SRC = Path("/tmp/e49/templates")
DST = Path("/mnt/j/conda_envs/foundry/DMS_Project/inputs/e49_templates")
DST.mkdir(parents=True, exist_ok=True)

JOBS = {  # pdb -> (field, chain to keep, out name); field: 'auth' or 'label' asym id
    "2QKE": ("auth", "A", "E2-01_tplA_2QKE"),
    "5JYT": ("auth", "A", "E2-01_tplB_5JYT"),
    "2OUG": ("auth", "A", "E2-02_tplA_2OUG"),
    "2LCL": ("auth", "A", "E2-02_tplB_2LCL"),
    "1DUJ": ("auth", "A", "E2-03_tplA_1DUJ"),
    "1S2H": ("auth", "A", "E2-03_tplB_1S2H"),
    "1J8I": ("auth", "A", "E2-04_tplA_1J8I"),
    "2JP1": ("auth", "A", "E2-04_tplB_2JP1"),
    "1PV6": ("auth", "A", "E2-10_tplA_1PV6"),
    "4OAA": ("auth", "A", "E2-10_tplB_4OAA"),
    "6M0J": ("label", "B", "RBD_tpl_6M0J"),  # auth E maps to TWO label copies (B, I); keep one
}

def col_indices(txt):
    m = re.search(r'loop_\s*\n((?:\s*_atom_site\.\S+\s*\n)+)', txt)
    cols = [c.strip() for c in m.group(1).splitlines() if c.strip()]
    idx = {c: i for i, c in enumerate(cols)}
    return idx

for pdb, (field, keep, out) in JOBS.items():
    txt = open(SRC / f"{pdb}.cif").read()
    idx = col_indices(txt)
    ia = idx.get("_atom_site.auth_asym_id") if field == "auth" else idx.get("_atom_site.label_asym_id")
    lines = txt.splitlines()
    out_lines = []
    n_in = n_out = 0
    for line in lines:
        if line.startswith(("ATOM ", "HETATM ", "ANISOU ")):
            n_in += 1
            p = line.split()
            if ia is not None and len(p) > ia and p[ia] == keep:
                out_lines.append(line)   # keep ORIGINAL chain id — re-chaining breaks
                n_out += 1               # atomworks assembly building (oper lists
            # drop other chains' atoms   # reference original asym ids; renamed chains
        else:                            # get silently dropped => 0 templated atoms)
            out_lines.append(line)
    dst = DST / f"{out}.cif"
    dst.write_text("\n".join(out_lines) + "\n")
    print(f"{pdb} chain {keep}: {n_in} -> {n_out} atom records -> {dst}")
