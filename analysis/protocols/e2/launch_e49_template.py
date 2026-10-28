#!/usr/bin/env python3
"""E49: MSA vs TEMPLATE — does an actual crystallographic template do to Z_II what
the MSA does? (2026-08-01, night queue after E47/E48)

Motivation: the corrected MSA story is "one-shot family-consensus write" — mild
eff-rank compression + functional gain + context-dependent state selection. The
retracted v1 story was "MSA = crystallographic template". E49 tests the residual
grain of truth in the old story directly: feed RF3 an ACTUAL template structure
(token-level template track, 64-bin pairwise distances) and ask whether the
template write resembles the MSA write, at two levels:

  RBD mode (geometry): does a WT-crystal template compress effective rank /
    rewrite mutation deltas like the MSA does?
    arms: tpl (6M0J chain E only) / msatpl (full15 a3m + 6M0J)
    subset: identical to E1 (stride-8, ~500 + WT), same 4-chain batching, same
    row+col slicing -> outputs are directly comparable to zii_e1_<arm> and zii/.
    out: /mnt/k/output_heads/rbd/zii_e49_{tpl,msatpl}/

  E2 mode (conformational state selection): on the 5 targets with both state
    structures (KaiB, RfaH, Mad2, XCL1, LacY), conditions tplA / tplB / msa_tplA
    (nomsa and msa already archived in /mnt/k/output_heads/e2/<tid>/).
    Readout downstream: distogram NLL state preference via e2_analyze.py.
    Template = single-chain re-chained ('Z') CIF in inputs/e49_templates/.
    out: /mnt/k/output_heads/e2/<tid>/<cond>/

Run: CUDA_VISIBLE_DEVICES=0 python -B launch_e49_template.py [rbd|e2] [arms/conds...]
PILOT: E49_PILOT=1 -> 1 RBD batch + E2-01 tplA only, prints token shapes.
"""
import os, sys, json, time, shutil, logging
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")  # 5090D (inference); must precede torch import
for p in ["/mnt/j/foundry-production/models/rf3/src", "/mnt/j/foundry-production/src"]:
    sys.path.insert(0, p)

import torch

PROJ = Path("/mnt/j/conda_envs/foundry/DMS_Project")
ARMS_DIR = PROJ / "inputs" / "e1_msa_arms"
TPL_DIR = PROJ / "inputs" / "e49_templates"
MANIFEST = PROJ / "data" / "sarbecovirus" / "manifest_SARS_CoV_2_WH1.csv"
TGT_JSON = PROJ / "inputs" / "e2_targets" / "e2_targets.json"
MSA_DIR = PROJ / "inputs" / "e2_targets" / "msas"
OUT_RBD = Path("/mnt/k/output_heads/rbd")
OUT_E2 = Path("/mnt/k/output_heads/e2")
STRIDE = int(os.environ.get("E1_STRIDE", "8"))
BATCH = 4
PILOT = os.environ.get("E49_PILOT", "0") == "1"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler("/tmp/e49_extract.log", mode="a")])
log = logging.getLogger("e49")

# E2 targets with both state structures downloaded (template chain re-chained to 'Z')
E2_TPL = {
    "E2-01": ("E2-01_tplA_2QKE.cif", "E2-01_tplB_5JYT.cif"),
    "E2-02": ("E2-02_tplA_2OUG.cif", "E2-02_tplB_2LCL.cif"),
    "E2-03": ("E2-03_tplA_1DUJ.cif", "E2-03_tplB_1S2H.cif"),
    "E2-04": ("E2-04_tplA_1J8I.cif", "E2-04_tplB_2JP1.cif"),
    "E2-10": ("E2-10_tplA_1PV6.cif", "E2-10_tplB_4OAA.cif"),
}


def init_engine():
    from omegaconf import OmegaConf
    from hydra import compose, initialize_config_dir
    from hydra.core.global_hydra import GlobalHydra
    from hydra.utils import instantiate
    config_dir = "/mnt/j/foundry-production/models/rf3/configs"
    GlobalHydra.instance().clear()
    overrides = ["num_steps=50", "diffusion_batch_size=2", "n_recycles=10",
                 "ckpt_path=/mnt/j/foundry_checkpoints/rf3_foundry_01_24_latest_remapped.ckpt"]
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
    log.info("RF3 engine ready.")
    return engine


def mut_pos(seq, wt):
    return next((j for j in range(min(len(seq), len(wt))) if seq[j] != wt[j]), -1)


# ---------------- RBD mode ----------------

def load_rbd_subset():
    import csv
    rows = list(csv.DictReader(open(MANIFEST)))
    sub = rows[::STRIDE]
    with open(ARMS_DIR / "full15.a3m") as f:
        f.readline()
        wt_seq = f.readline().strip()
    sub = [dict(mutant_id="WT", protein=wt_seq)] + sub
    return sub, wt_seq


def run_rbd_arm(engine, arm, subset, wt_seq):
    """arm: 'tpl' (template only) or 'msatpl' (full15 MSA + template)."""
    out_dir = OUT_RBD / f"zii_e49_{arm}"
    out_dir.mkdir(parents=True, exist_ok=True)
    gold_lines = open(ARMS_DIR / "full15.a3m").readlines() if arm == "msatpl" else None
    tpl_cif = TPL_DIR / "RBD_tpl_6M0J.cif"
    todo = [m for m in subset if not (out_dir / f"{m['mutant_id']}_zii.pt").exists()]
    if PILOT:
        todo = todo[:BATCH]
    log.info(f"[rbd:{arm}] {len(subset) - len(todo)} done, {len(todo)} to go -> {out_dir}")
    n_done, t0 = 0, time.time()
    for b in range(0, len(todo), BATCH):
        batch = todo[b:b + BATCH]
        work = Path(f"/tmp/e49_{arm}_{b:04d}")
        work.mkdir(parents=True, exist_ok=True)
        try:
            chain_ids = ['P', 'Q', 'R', 'S'][:len(batch)]   # avoid template chain 'B' (and any assembly-generated ids)
            components = []
            for cid, m in zip(chain_ids, batch):
                comp = {"chain_id": cid, "seq": m["protein"]}
                if gold_lines is not None:
                    lines = gold_lines.copy()
                    lines[0] = f">{m['mutant_id']}\n"
                    lines[1] = f"{m['protein']}\n"
                    d = work / cid
                    d.mkdir(parents=True, exist_ok=True)
                    msa_file = d / "t000_.msa0.a3m"
                    with open(msa_file, "w") as f:
                        f.writelines(lines)
                    comp["msa_path"] = str(msa_file)
                components.append(comp)
            components.append({"path": str(tpl_cif)})
            spec = {"name": f"e49_{arm}_{b:04d}", "components": components,
                    "template_selection": ["B"]}   # 6M0J filtered to label chain B = single RBD copy
            jp = work / "input.json"
            with open(jp, "w") as f:
                json.dump([spec], f, indent=2)
            os.environ["RF3_ZII_COUNTER"] = "0"
            os.environ["RF3_ZII_PATH"] = str(work / "zii.pt")
            os.environ["RF3_ZII_FINAL_PATH"] = str(work / "zii_final.pt")
            engine.run(inputs=str(jp), out_dir=str(work / "out"),
                       dump_predictions=False, dump_trajectories=False,
                       one_model_per_file=False, skip_existing=False,
                       template_selection=["B"])
            zf = work / "zii_final.pt"
            if not zf.exists():
                raise RuntimeError("no zii_final.pt")
            z = torch.load(zf, map_location="cpu", weights_only=True)
            lens = [len(m["protein"]) for m in batch]
            total = z.shape[-2]
            exp = sum(lens)
            if b == 0:
                log.info(f"[rbd:{arm}] zii tokens {total} = query {exp} + template {total - exp}")
            assert total >= exp, f"token mismatch {total} < {exp}"
            # query chains occupy the first `exp` tokens (component order); template
            # tokens follow and are excluded from the per-mutant slices
            bounds = [0]
            for sl in lens:
                bounds.append(bounds[-1] + sl)
            for i, m in enumerate(batch):
                blk = z[bounds[i]:bounds[i + 1], bounds[i]:bounds[i + 1], :]
                p = mut_pos(m["protein"], wt_seq)
                if m['mutant_id'] == 'WT':
                    torch.save(blk.half(), out_dir / f"{m['mutant_id']}_zii.pt")
                elif p < 0:  # WT-identical replicate: position 0 (old-pipeline convention)
                    torch.save(torch.stack([blk[0, :, :], blk[:, 0, :]]).half(), out_dir / f"{m['mutant_id']}_zii.pt")
                else:
                    torch.save(torch.stack([blk[p, :, :], blk[:, p, :]]).half(), out_dir / f"{m['mutant_id']}_zii.pt")
            n_done += len(batch)
            if n_done % 40 < BATCH:
                rate = n_done / max(time.time() - t0, 1e-6)
                log.info(f"[rbd:{arm}] {n_done}/{len(todo)} {rate:.2f} mut/s "
                         f"ETA {(len(todo) - n_done) / max(rate, 1e-6) / 60:.0f}min")
        except Exception as e:
            log.error(f"[rbd:{arm}] batch {b} FAILED: {e}", exc_info=True)
            if PILOT:
                raise
        finally:
            shutil.rmtree(work, ignore_errors=True)
            torch.cuda.empty_cache()
    log.info(f"[rbd:{arm}] DONE {n_done}/{len(todo)} in {(time.time() - t0) / 60:.0f}min")


# ---------------- E2 mode ----------------

def parse_range(s):
    import re
    m = re.match(r'\s*(\d+)\s*-\s*(\d+)', s or '')
    return (int(m.group(1)) - 1, int(m.group(2))) if m else None


def run_e2_one(engine, tid, seq, a3m_lines, tpl_cif, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir / 'disto.pt').exists():
        return
    work = Path(f"/tmp/e49w_{tid}_{out_dir.name}")
    work.mkdir(parents=True, exist_ok=True)
    try:
        comp = {"chain_id": "Q", "seq": seq}   # 'Q' avoids chain-id clash with template chain 'A' (hard ValueError)
        if a3m_lines is not None:
            d = work / "Q"
            d.mkdir(exist_ok=True)
            mf = d / "t000_.msa0.a3m"
            with open(mf, "w") as f:
                f.writelines(a3m_lines)
            comp["msa_path"] = str(mf)
        components = [comp, {"path": str(tpl_cif)}]
        jp = work / "input.json"
        with open(jp, "w") as f:
            json.dump([{"name": tid, "components": components,
                        "template_selection": ["A"]}], f, indent=2)   # all E2 templates filtered to auth chain A
        os.environ["RF3_ZII_COUNTER"] = "0"
        os.environ["RF3_ZII_PATH"] = str(work / "zii.pt")
        os.environ["RF3_ZII_FINAL_PATH"] = str(work / "zii_final.pt")
        os.environ["RF3_DISTOGRAM_PATH"] = str(work / "disto.pt")
        engine.run(inputs=str(jp), out_dir=str(work / "out"),
                   dump_predictions=False, dump_trajectories=False,
                   one_model_per_file=False, skip_existing=False,
                   template_selection=["A"])
        zf, dp = work / "zii_final.pt", work / "disto.pt"
        if not zf.exists() or not dp.exists():
            raise RuntimeError(f"missing outputs: zii={zf.exists()} disto={dp.exists()}")
        z = torch.load(zf, map_location="cpu", weights_only=True)
        if PILOT:
            log.info(f"[PILOT e2:{tid}/{out_dir.name}] zii shape {tuple(z.shape)}, query L={len(seq)}")
        shutil.copy(zf, out_dir / 'zii.pt')
        shutil.copy(dp, out_dir / 'disto.pt')
        with open(out_dir / 'meta.json', 'w') as f:
            json.dump({'target': tid, 'cond': out_dir.name, 'L': len(seq),
                       'template': str(tpl_cif)}, f)
        log.info(f'{tid}/{out_dir.name} done (L={len(seq)})')
    except Exception as e:
        log.error(f'{tid}/{out_dir.name} FAILED: {e}', exc_info=True)
        if PILOT:
            raise
    finally:
        os.environ.pop("RF3_DISTOGRAM_PATH", None)
        shutil.rmtree(work, ignore_errors=True)
        torch.cuda.empty_cache()


def run_e2(engine, conds):
    targets = {t['id']: t for t in json.load(open(TGT_JSON))['targets']}
    for tid, (cifa, cifb) in E2_TPL.items():
        t = targets[tid]
        a3m_file = MSA_DIR / f'{tid}.a3m'
        a3m_lines = open(a3m_file).readlines()
        full_seq = a3m_lines[1].strip()
        dom = parse_range(t.get('construct', {}).get('boundaries', ''))
        seq = full_seq
        if len(full_seq) > 600 and dom:
            lo, hi = dom
            seq = full_seq[lo:hi]
            a3m_lines = [a3m_lines[0]] + [l[lo:hi] + '\n' if not l.startswith('>')
                                          else l for l in a3m_lines[1:]]
            log.info(f'{tid}: domain-sliced to {lo+1}-{hi} (L={len(seq)})')
        out_t = OUT_E2 / tid
        todo_conds = conds or ["tplA", "tplB", "msa_tplA"]
        if PILOT:
            todo_conds = ["tplA"]
            if tid != "E2-01":
                continue
        for cond in todo_conds:
            if cond == "tplA":
                run_e2_one(engine, tid, seq, None, TPL_DIR / cifa, out_t / 'tplA')
            elif cond == "tplB":
                run_e2_one(engine, tid, seq, None, TPL_DIR / cifb, out_t / 'tplB')
            elif cond == "msa_tplA":
                run_e2_one(engine, tid, seq, a3m_lines, TPL_DIR / cifa, out_t / 'msa_tplA')


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "rbd"
    rest = sys.argv[2:]
    engine = init_engine()
    if mode == "rbd":
        subset, wt_seq = load_rbd_subset()
        log.info(f"subset: {len(subset)} entries (incl WT), stride={STRIDE}")
        for arm in (rest or ["tpl", "msatpl"]):
            run_rbd_arm(engine, arm, subset, wt_seq)
    elif mode == "e2":
        run_e2(engine, rest)
    log.info("E49 DONE")


if __name__ == "__main__":
    main()
