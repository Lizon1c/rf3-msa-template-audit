#!/usr/bin/env python3
"""E428: structure-level validation of M2's distogram state-preference claims.

M2's hierarchy results (template > MSA > sequence) rest on distogram NLL
preference. Here we SAMPLE actual structures with RF3's diffusion module
(50 steps, diffusion_batch_size=2, 8 reps = 16 structures per arm) and score
each sampled structure against BOTH experimental state references (C-alpha
RMSD after Needleman-Wunsch sequence alignment + Kabsch superposition).
Sign convention matches M2: positive = closer to state A (GA88-fold for the
designed pair), so structure DeltaRMSD = RMSD_B - RMSD_A should track pref.

Cells (30 total = 240 forwards):
  E2-01 KaiB / E2-02 RfaH / E2-03 Mad2 / E2-04 XCL1:
      arms nomsa, msa, tplA, tplB, msa_tplA   (state A/B = E49 template CIFs)
  GA88 / GB88 (E2-06):
      arms nomsa, msa, tplOwn, tplOther, msa_tplOther
      (references 2JWS = GA-fold = A, 2JWU = GB-fold = B for both sequences)

PORTABLE LAYOUT (this repository):
  analyze mode is GPU-free and reads/writes repo-relative paths by default:
    MSAs        data/inputs/e2_targets/msas
    templates   data/inputs/e49_templates
    structures  data/structures/e428/{tid}/{arm}/rep{r}/*.cif
    results     data/results/e428/results.json  (override: --output)
  Optional overrides: --msa-dir, --tpl-dir, --structs-dir, --output, or the
  E428_MSA_DIR / E428_TPL_DIR / E428_STRUCT_DIR environment variables.

RUN MODE (RF3 inference, requires GPU) is NOT portable by itself; it needs:
    RF3_ROOT        RF3/foundry model code root (src/ and models/rf3 under it)
    RF3_CHECKPOINT  RF3 remapped checkpoint (.ckpt)
    RF3_CONFIG_DIR  hydra configs (default: $RF3_ROOT/models/rf3/configs)
  GPU choice is the caller's: export CUDA_VISIBLE_DEVICES yourself. To keep
  the original card-model safety guard, set E428_EXPECT_GPU (e.g. "6000").
  Resume: cells with 16 CIFs are skipped.

Run:    python -B e428_structure_sampling.py run [--only E2-01/nomsa] [--reps 8]
Score:  python -B e428_structure_sampling.py analyze [--output PATH]
"""
import argparse, json, os, shutil, sys, time
from pathlib import Path

_ANALYSIS_DIR = Path(__file__).resolve().parents[2]
if str(_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(_ANALYSIS_DIR))
import repo_paths

import numpy as np

try:  # torch is only needed for run mode (RF3 inference); analyze is GPU- and torch-free
    import torch
except ImportError:  # pragma: no cover
    torch = None

MSA_DIR = repo_paths.env_path("E428_MSA_DIR", repo_paths.E2_MSA_DIR)
TPL_DIR = repo_paths.env_path("E428_TPL_DIR", repo_paths.E49_TPL_DIR)
STRUCTS = repo_paths.env_path("E428_STRUCT_DIR", repo_paths.E428_STRUCTURES_DIR)
RESULTS_OUT = repo_paths.env_path("E428_RESULTS_OUT", repo_paths.E428_RESULTS)
OUT = repo_paths.env_path("E428_RUN_OUT", STRUCTS)  # run-mode CIF destination
t0 = time.time()


def log(m):
    print(f'[{time.time()-t0:7.1f}s] {m}', flush=True)


E2_TPL = {
    "E2-01": ("E2-01_tplA_2QKE.cif", "E2-01_tplB_5JYT.cif"),
    "E2-02": ("E2-02_tplA_2OUG.cif", "E2-02_tplB_2LCL.cif"),
    "E2-03": ("E2-03_tplA_1DUJ.cif", "E2-03_tplB_1S2H.cif"),
    "E2-04": ("E2-04_tplA_1J8I.cif", "E2-04_tplB_2JP1.cif"),
}
# M2 distogram preferences (range-matched, Table 2 of m2_state_selection) for sign comparison
EXPECTED = {
    "E2-01": dict(nomsa=+0.72, msa=-1.57, tplA=+3.59, tplB=-3.00, msa_tplA=+3.36),
    "E2-02": dict(nomsa=+0.10, msa=-2.89, tplA=+1.49, tplB=+0.30, msa_tplA=+1.45),
    "E2-03": dict(nomsa=-0.05, msa=-0.96, tplA=-0.17, tplB=-1.87, msa_tplA=-0.17),
    "E2-04": dict(nomsa=-0.78, msa=+2.09, tplA=-1.49, tplB=-0.86, msa_tplA=-1.44),
    "GA88": dict(nomsa=+0.38, msa=+2.74, tplOwn=+4.88, tplOther=-4.30, msa_tplOther=-4.24),
    "GB88": dict(nomsa=-0.23, msa=-2.66, tplOwn=-4.34, tplOther=+4.75, msa_tplOther=+4.58),
}


def init_engine():
    from omegaconf import OmegaConf
    from hydra import compose, initialize_config_dir
    from hydra.core.global_hydra import GlobalHydra
    from hydra.utils import instantiate
    rf3_root, ckpt = repo_paths.require_rf3()
    for p in [str(rf3_root / "models" / "rf3" / "src"), str(rf3_root / "src")]:
        if p not in sys.path:
            sys.path.insert(0, p)
    config_dir = repo_paths.rf3_config_dir()
    if not config_dir or not Path(config_dir).is_dir():
        raise RuntimeError(f"RF3 config dir not found: {config_dir!r}; set RF3_CONFIG_DIR")
    GlobalHydra.instance().clear()
    overrides = ["num_steps=50", "diffusion_batch_size=2", "n_recycles=10",
                 f"ckpt_path={ckpt}"]
    with initialize_config_dir(config_dir=config_dir, version_base="1.3"):
        cfg = compose(config_name="inference", overrides=overrides)
    run_keys = {"inputs", "out_dir", "dump_predictions", "dump_trajectories",
                "one_model_per_file", "annotate_b_factor_with_plddt",
                "sharding_pattern", "skip_existing", "template_selection",
                "ground_truth_conformer_selection", "cyclic_chains", "add_missing_atoms"}
    cfg_dict = OmegaConf.to_container(cfg, resolve=True)
    init_cfg = OmegaConf.create({k: v for k, v in cfg_dict.items() if k not in run_keys})
    engine = instantiate(init_cfg, _convert_="partial", _recursive_=False)
    engine.initialize()
    log(f"RF3 engine ready on cuda:{torch.cuda.current_device()} = "
        f"{torch.cuda.get_device_name(torch.cuda.current_device())}")
    return engine


def cell_list():
    cells = []
    for tid, (cifa, cifb) in E2_TPL.items():
        for arm in ["nomsa", "msa", "tplA", "tplB", "msa_tplA"]:
            cells.append((tid, arm))
    for tid in ["GA88", "GB88"]:
        for arm in ["nomsa", "msa", "tplOwn", "tplOther", "msa_tplOther"]:
            cells.append((tid, arm))
    return cells


def get_seq(tid):
    """Query sequence: a3m row 1 (GA88 = E2-06); GB88 from its reference CIF."""
    if tid == "GB88":
        return cif_sequence(TPL_DIR / "E2-06_tplGB_2JWU.cif")
    a3m = "E2-06.a3m" if tid == "GA88" else f"{tid}.a3m"
    return open(MSA_DIR / a3m).readlines()[1].strip()


def get_msa(tid):
    lines = open(MSA_DIR / ("E2-06.a3m" if tid in ("GA88", "GB88") else f"{tid}.a3m")).readlines()
    if tid == "GB88":
        lines = lines.copy()
        lines[0] = ">GB88\n"
        lines[1] = get_seq("GB88") + "\n"
    return lines


def get_tpl(tid, arm):
    if arm not in ("tplA", "tplB", "msa_tplA", "tplOwn", "tplOther", "msa_tplOther"):
        return None
    if tid in E2_TPL:
        cif = E2_TPL[tid][0] if arm in ("tplA", "msa_tplA") else E2_TPL[tid][1]
        return TPL_DIR / cif
    own = TPL_DIR / ("E2-06_tplGA_2JWS.cif" if tid == "GA88" else "E2-06_tplGB_2JWU.cif")
    other = TPL_DIR / ("E2-06_tplGB_2JWU.cif" if tid == "GA88" else "E2-06_tplGA_2JWS.cif")
    if arm in ("tplA", "msa_tplA", "tplOwn"):
        return own
    return other


def run_cell(engine, tid, arm, rep):
    dest = OUT / tid / arm / f"rep{rep}"
    if dest.exists() and len(list(dest.glob("*.cif"))) >= 2:
        return
    dest.mkdir(parents=True, exist_ok=True)
    seq = get_seq(tid)
    work = Path(f"/tmp/e428w_{tid}_{arm}_{rep}")
    work.mkdir(parents=True, exist_ok=True)
    try:
        comp = {"chain_id": "Q", "seq": seq}
        tpl = get_tpl(tid, arm)
        use_msa = arm in ("msa", "msa_tplA", "msa_tplOther")
        if use_msa:
            d = work / "Q"
            d.mkdir(exist_ok=True)
            mf = d / "t000_.msa0.a3m"
            with open(mf, "w") as f:
                f.writelines(get_msa(tid))
            comp["msa_path"] = str(mf)
        components = [comp] + ([{"path": str(tpl)}] if tpl else [])
        spec = {"name": f"e428_{tid}_{arm}_r{rep}", "components": components}
        if tpl:
            spec["template_selection"] = ["A"]
        jp = work / "input.json"
        with open(jp, "w") as f:
            json.dump([spec], f, indent=2)
        engine.run(inputs=str(jp), out_dir=str(work / "out"),
                   dump_predictions=True, dump_trajectories=False,
                   one_model_per_file=False, skip_existing=False,
                   template_selection=["A"] if tpl else [])
        cifs = list((work / "out").rglob("*.cif"))
        if not cifs:
            raise RuntimeError("no CIF predictions dumped")
        for c in cifs:
            shutil.copy(c, dest / c.name)
        with open(dest / "meta.json", "w") as f:
            json.dump({"tid": tid, "arm": arm, "rep": rep, "L": len(seq),
                       "n_cif": len(cifs), "template": str(tpl) if tpl else None}, f)
        log(f"{tid}/{arm}/rep{rep}: {len(cifs)} structures (L={len(seq)})")
    except Exception as e:
        log(f"{tid}/{arm}/rep{rep} FAILED: {e}")
        raise
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if torch is not None:
            torch.cuda.empty_cache()


# ── scoring ──
AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def read_ca(path, chain=None):
    from biotite.structure.io.pdbx import CIFFile, get_structure
    st = get_structure(CIFFile.read(str(path)), model=1)
    if chain is not None:
        st = st[st.chain_id == chain]
    st = st[(st.atom_name == "CA")]
    seq = ''.join(AA3.get(r, 'X') for r in st.res_name)
    return seq, st.coord.astype(np.float64), list(st.res_id)


def cif_sequence(path):
    seq, _, _ = read_ca(path)
    return seq


def nw_align(a, b):
    """Needleman-Wunsch, match +2 / mismatch -1 / gap -2. Returns index pairs."""
    n, m = len(a), len(b)
    F = np.full((n + 1, m + 1), -10**9, dtype=np.int32)
    P = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag 1 up 2 left
    F[0, :] = -2 * np.arange(m + 1)
    F[:, 0] = -2 * np.arange(n + 1)
    P[1:, 0] = 1
    P[0, 1:] = 2
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = F[i - 1, j - 1] + (2 if ai == b[j - 1] else -1)
            up = F[i - 1, j] - 2
            left = F[i, j - 1] - 2
            best = max(diag, up, left)
            F[i, j] = best
            P[i, j] = 0 if best == diag else (1 if best == up else 2)
    pairs = []
    i, j = n, m
    while i > 0 or j > 0:
        p = P[i, j]
        if i > 0 and j > 0 and p == 0:
            pairs.append((i - 1, j - 1))
            i -= 1; j -= 1
        elif i > 0 and (j == 0 or p == 1):
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def kabsch_rmsd(P, Q):
    P = P - P.mean(0)
    Q = Q - Q.mean(0)
    C = P.T @ Q
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    U = V @ D @ Wt
    diff = P @ U - Q
    return float(np.sqrt((diff ** 2).sum() / len(P)))


def rmsd_to_ref(pred_seq, pred_coord, ref_seq, ref_coord):
    pairs = nw_align(pred_seq, ref_seq)
    P = np.array([pred_coord[i] for i, j in pairs])
    Q = np.array([ref_coord[j] for i, j in pairs])
    ident = sum(1 for i, j in pairs if pred_seq[i] == ref_seq[j]) / max(len(pairs), 1)
    return kabsch_rmsd(P, Q), len(pairs), ident


def rmsd_shared(pred_seq, pred_coord, seqA, coA, seqB, coB):
    """RMSD to both references on the SAME residue subset (positions aligned
    to both references), mirroring M2's shared-residue-universe convention."""
    mapA = dict(nw_align(pred_seq, seqA))
    mapB = dict(nw_align(pred_seq, seqB))
    shared = sorted(set(mapA) & set(mapB))
    if len(shared) < 10:
        return None
    P = pred_coord[shared]
    QA = coA[[mapA[i] for i in shared]]
    QB = coB[[mapB[i] for i in shared]]
    ra = kabsch_rmsd(P, QA)
    rb = kabsch_rmsd(P, QB)
    ident = float(np.mean([pred_seq[i] == seqA[mapA[i]] for i in shared]))
    return ra, rb, len(shared), ident


def analyze(output=None):
    rows = []
    for tid, arm in cell_list():
        refA = TPL_DIR / (E2_TPL[tid][0] if tid in E2_TPL else "E2-06_tplGA_2JWS.cif")
        refB = TPL_DIR / (E2_TPL[tid][1] if tid in E2_TPL else "E2-06_tplGB_2JWU.cif")
        seqA, coA, _ = read_ca(refA)
        seqB, coB, _ = read_ca(refB)
        for rep_dir in sorted((STRUCTS / tid / arm).glob("rep*")):
            for cif in sorted(rep_dir.glob("*sample-*_model.cif")):
                ps, pc, _ = read_ca(cif, chain="Q")
                out = rmsd_shared(ps, pc, seqA, coA, seqB, coB)
                if out is None:
                    continue
                ra, rb, ns, ids = out
                rows.append(dict(tid=tid, arm=arm, rep=rep_dir.name, cif=cif.name,
                                 rmsdA=ra, rmsdB=rb, n_shared=ns, ident=round(ids, 3),
                                 closer_to='A' if ra < rb else 'B'))
    out_path = Path(output) if output else RESULTS_OUT
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(rows, f, indent=1)
    print(f"\n{'target':<7}{'arm':<14}{'n':>3} {'RMSD_A':>8} {'RMSD_B':>8} "
          f"{'dRMSD':>8} {'%A':>5} {'pref(M2)':>9}  sign_match")
    for tid, arm in cell_list():
        rs = [r for r in rows if r['tid'] == tid and r['arm'] == arm]
        if not rs:
            continue
        ra = np.mean([r['rmsdA'] for r in rs])
        rb = np.mean([r['rmsdB'] for r in rs])
        d = rb - ra
        fa = np.mean([r['closer_to'] == 'A' for r in rs])
        ns = np.mean([r['n_shared'] for r in rs])
        pref = EXPECTED[tid][arm]
        match = (d > 0) == (pref > 0)
        print(f"{tid:<7}{arm:<14}{len(rs):>3} {ra:>8.2f} {rb:>8.2f} {d:>+8.2f} "
              f"{fa:>5.2f} {pref:>+9.2f}  {'OK' if match else 'MISMATCH'} (n_shared={ns:.0f})")
    log(f"saved {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "analyze"])
    ap.add_argument("--only", default="", help="e.g. E2-01/nomsa (single cell)")
    ap.add_argument("--targets", default="", help="comma list, e.g. GA88,GB88")
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--output", default=None,
                    help="analyze: results.json destination (default: repo data/results/e428/results.json)")
    args = ap.parse_args()
    if args.mode == "analyze":
        analyze(output=args.output)
        return
    repo_paths.require_rf3()
    if torch is None:
        raise RuntimeError("run mode needs torch (with CUDA); analyze mode does not")
    expect_gpu = os.environ.get("E428_EXPECT_GPU", "")
    if expect_gpu:
        name = torch.cuda.get_device_name(torch.cuda.current_device())
        if expect_gpu not in name:
            raise RuntimeError(f"wrong GPU: {name}; want {expect_gpu} "
                               f"(set E428_EXPECT_GPU to change or clear to skip the check)")
    if not torch.cuda.is_available():
        raise RuntimeError("run mode needs a CUDA GPU (analyze mode does not)")
    engine = init_engine()
    for tid, arm in cell_list():
        if args.only and f"{tid}/{arm}" != args.only:
            continue
        if args.targets and tid not in args.targets.split(','):
            continue
        for rep in range(args.reps):
            run_cell(engine, tid, arm, rep)
    log("E428 run DONE")


if __name__ == "__main__":
    main()
