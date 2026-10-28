#!/usr/bin/env python3
"""E2 runner: RF3 forward per E2 target, no-MSA and with-MSA, save Z_II + distogram.

Per target (inputs/e2_targets/e2_targets.json + msas/{id}.a3m from UniRef30 build):
  - fold sequence: UniProt WT; if len > 600 (kinases Abl/EGFR), slice to the
    construct.domain_range field (fallback: construct.boundaries) and slice the
    a3m columns correspondingly.
  - conditions: nomsa (seq-only), msa (a3m with row0 = query)
  - saves /mnt/k/output_heads/e2/{id}/{cond}/zii.pt + disto.pt + meta.json
Run: CUDA_VISIBLE_DEVICES=1 python -B launch_e2_targets.py [E2-01 E2-03 ...]
"""
import os, sys, json, time, shutil, logging
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "1")
for p in ["/mnt/j/foundry-production/models/rf3/src", "/mnt/j/foundry-production/src"]:
    sys.path.insert(0, p)

import torch

PROJ = Path("/mnt/j/conda_envs/foundry/DMS_Project")
TGT_JSON = PROJ / "inputs" / "e2_targets" / "e2_targets.json"
MSA_DIR = PROJ / "inputs" / "e2_targets" / "msas"
SEQ_DIR = Path("/tmp/e2/seqs")
OUT_BASE = Path("/mnt/k/output_heads/e2")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler("/tmp/e2_run.log", mode="a")])
log = logging.getLogger("e2")


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


def parse_range(s):
    """'1-108' -> (0, 108) 0-indexed half-open; returns None if unparseable."""
    import re
    m = re.match(r'\s*(\d+)\s*-\s*(\d+)', s or '')
    return (int(m.group(1)) - 1, int(m.group(2))) if m else None


def run_one(engine, tid, seq, a3m_lines, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir / 'disto.pt').exists():
        return
    work = Path(f"/tmp/e2w_{tid}_{out_dir.name}")
    work.mkdir(parents=True, exist_ok=True)
    try:
        comp = {"chain_id": "A", "seq": seq}
        if a3m_lines is not None:
            d = work / "A"
            d.mkdir(exist_ok=True)
            mf = d / "t000_.msa0.a3m"
            with open(mf, "w") as f:
                f.writelines(a3m_lines)
            comp["msa_path"] = str(mf)
        jp = work / "input.json"
        with open(jp, "w") as f:
            json.dump([{"name": tid, "components": [comp]}], f, indent=2)
        os.environ["RF3_ZII_COUNTER"] = "0"
        os.environ["RF3_ZII_PATH"] = str(work / "zii.pt")
        os.environ["RF3_ZII_FINAL_PATH"] = str(work / "zii_final.pt")
        os.environ["RF3_DISTOGRAM_PATH"] = str(work / "disto.pt")
        engine.run(inputs=str(jp), out_dir=str(work / "out"),
                   dump_predictions=False, dump_trajectories=False,
                   one_model_per_file=False, skip_existing=False,
                   template_selection=[])
        zf, dp = work / "zii_final.pt", work / "disto.pt"
        if not zf.exists() or not dp.exists():
            raise RuntimeError(f"missing outputs: zii={zf.exists()} disto={dp.exists()}")
        shutil.copy(zf, out_dir / 'zii.pt')
        shutil.copy(dp, out_dir / 'disto.pt')
        with open(out_dir / 'meta.json', 'w') as f:
            json.dump({'target': tid, 'cond': out_dir.name, 'L': len(seq)}, f)
        log.info(f'{tid}/{out_dir.name} done (L={len(seq)})')
    except Exception as e:
        log.error(f'{tid}/{out_dir.name} FAILED: {e}', exc_info=True)
    finally:
        os.environ.pop("RF3_DISTOGRAM_PATH", None)
        shutil.rmtree(work, ignore_errors=True)
        torch.cuda.empty_cache()


def main():
    targets = {t['id']: t for t in json.load(open(TGT_JSON))['targets']}
    ids = sys.argv[1:] or list(targets.keys())
    engine = init_engine()
    for tid in ids:
        t = targets.get(tid)
        if t is None:
            continue
        fa = list(SEQ_DIR.glob(f'{tid}_*.fa'))
        if not fa:
            log.warning(f'{tid}: no sequence, skip')
            continue
        lines = open(fa[0]).read().split('\n')
        full_seq = lines[1].strip()
        dom = parse_range(t.get('construct', {}).get('boundaries', ''))
        seq = full_seq
        a3m_file = MSA_DIR / f'{tid}.a3m'
        a3m_lines = open(a3m_file).readlines() if a3m_file.exists() else None
        if len(full_seq) > 600 and dom:
            lo, hi = dom
            seq = full_seq[lo:hi]
            if a3m_lines:
                a3m_lines = [a3m_lines[0]] + [l[lo:hi] + '\n' if not l.startswith('>')
                                              else l for l in a3m_lines[1:]]
            log.info(f'{tid}: domain-sliced to {lo+1}-{hi} (L={len(seq)})')
        out_t = OUT_BASE / tid
        run_one(engine, tid, seq, None, out_t / 'nomsa')
        if a3m_lines:
            run_one(engine, tid, seq, a3m_lines, out_t / 'msa')


if __name__ == "__main__":
    main()
