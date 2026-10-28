#!/usr/bin/env python3
"""E49b: GA88/GB88 template arms — sequence vs template conflict on the designed
fold-switch pair (2026-08-01, extends E49-E2 and e2_06_ga_gb).

Archived: GA88/GB88 each fold into their OWN fold (sequence dominates); MSA
SHARPENS the intrinsic preference (GA88 +0.38->+2.74, GB88 -0.23->-2.66).
E49-E2 (5 targets): ground-truth template DOMINATES the MSA in conflicts.
Here the conflict is sequence vs template:
  GA88_tplOwn   : GA88 + tpl 2JWS (own fold)  — sharpening like MSA?
  GA88_tplOther : GA88 + tpl 2JWU (OTHER fold) — can a template flip the sequence?
  GA88_msa_tplOther : GA88 + GA-family MSA + tpl 2JWU — seq+MSA vs template
  (symmetric three for GB88)
Readout downstream: distogram NLL preference vs 2JWS/2JWU distance npys.
Run: CUDA_VISIBLE_DEVICES=0 python -B e49b_ga_gb_template.py
"""
import os, sys, json, shutil, logging
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
for p in ["/mnt/j/foundry-production/models/rf3/src", "/mnt/j/foundry-production/src"]:
    sys.path.insert(0, p)
import torch

PROJ = Path("/mnt/j/conda_envs/foundry/DMS_Project")
A3M = PROJ / "inputs" / "e2_targets" / "msas" / "E2-06.a3m"
TPL_DIR = PROJ / "inputs" / "e49_templates"
SEQ_DIR = Path("/tmp/e2/seqs")
OUT_BASE = Path("/mnt/k/output_heads/e2/E2-06")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler("/tmp/e49b_run.log", mode="a")])
log = logging.getLogger("e49b")


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


def run(engine, tag, seq, a3m_lines, tpl_cif, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir / 'disto.pt').exists():
        log.info(f'{tag} exists, skip')
        return
    work = Path(f"/tmp/e49b_{tag}")
    work.mkdir(parents=True, exist_ok=True)
    try:
        comp = {"chain_id": "Q", "seq": seq}   # 'Q' avoids clash with template chain 'A'
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
            json.dump([{"name": tag, "components": components,
                        "template_selection": ["A"]}], f, indent=2)
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
            raise RuntimeError(f"missing: zii={zf.exists()} disto={dp.exists()}")
        shutil.copy(zf, out_dir / 'zii.pt')
        shutil.copy(dp, out_dir / 'disto.pt')
        with open(out_dir / 'meta.json', 'w') as f:
            json.dump({'tag': tag, 'L': len(seq), 'template': str(tpl_cif)}, f)
        log.info(f'{tag} done (L={len(seq)})')
    except Exception as e:
        log.error(f'{tag} FAILED: {e}', exc_info=True)
    finally:
        os.environ.pop("RF3_DISTOGRAM_PATH", None)
        shutil.rmtree(work, ignore_errors=True)
        torch.cuda.empty_cache()


def main():
    ga88 = open(SEQ_DIR / "E2-06_GA88.fa").read().split('\n')[1].strip()
    gb88 = open(SEQ_DIR / "E2-06_GB88.fa").read().split('\n')[1].strip()
    ga_lines = open(A3M).readlines()
    gb_lines = ga_lines.copy()
    gb_lines[0] = ">GB88\n"
    gb_lines[1] = gb88 + "\n"
    tplGA = TPL_DIR / "E2-06_tplGA_2JWS.cif"
    tplGB = TPL_DIR / "E2-06_tplGB_2JWU.cif"
    engine = init_engine()
    run(engine, "GA88_tplOwn", ga88, None, tplGA, OUT_BASE / "GA88_tplOwn")
    run(engine, "GA88_tplOther", ga88, None, tplGB, OUT_BASE / "GA88_tplOther")
    run(engine, "GA88_msa_tplOther", ga88, ga_lines, tplGB, OUT_BASE / "GA88_msa_tplOther")
    run(engine, "GB88_tplOwn", gb88, None, tplGB, OUT_BASE / "GB88_tplOwn")
    run(engine, "GB88_tplOther", gb88, None, tplGA, OUT_BASE / "GB88_tplOther")
    run(engine, "GB88_msa_tplOther", gb88, gb_lines, tplGA, OUT_BASE / "GB88_msa_tplOther")
    log.info("E49b all done")


if __name__ == "__main__":
    main()
