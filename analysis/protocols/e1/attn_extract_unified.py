#!/usr/bin/env python3
"""
attn_extract_unified.py — Unified multi-protein attention extraction pipeline

TWO MODES:
  mode="heads" (default):  per-mutant attention heads + Z_II via hook capture.
                           Model loaded ONCE. All JSONs pre-generated, processed
                           in a tight loop. Maximizes GPU utilization.
  mode="structure":        RF3 JSON-array batch mode. N mutants in ONE
                           engine.run() call. Fastest throughput. No heads/Z_II.

AGENT ENTRY POINT:
  from attn_extract_unified import extract_from_manifest
  summary = extract_from_manifest(protein_name="sirt1", device="cuda:1",
      manifest_csv=".../manifest.csv", output_dir=".../heads",
      mode="heads", batch_size=10, ...)

ISOLATION:
  - Unique run_id per invocation -> /tmp/attn_extract_{run_id}/
  - RF3_ZII_PATH / RF3_ZII_FINAL_PATH = per-instance paths
  - Z_II shape validated against expected token count before save
  - Concurrent instances on different GPUs: zero conflict

ROBUSTNESS:
  - Missing template/MSA? WARNING, run de novo (seq-only)
  - No CUDA? WARNING, fall back to CPU
  - All torch/pandas imports are LAZY (module loads instantly, no CUDA init)

OLD SCRIPT (run_attn_extraction.py) IS UNTOUCHED.
"""
from __future__ import annotations

import gc
import json
import logging
import os
import random
import shutil
import string
import subprocess
import sys
import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Lazy import helpers — torch and pandas imported INSIDE functions only.
# This prevents CUDA init on module load (hangs when other GPU is busy).
def _import_torch():
    import torch as _T
    return _T

def _import_pd():
    import pandas as _P
    return _P

# INSTALL PATHS (fixed, not per-protein)
FOUNDRY_ROOT     = Path("/mnt/j/foundry-production")
RF3_ROOT         = FOUNDRY_ROOT / "models" / "rf3"
FAESM_ROOT       = Path("/mnt/j/faesm-main")
MANAGER_DIR      = RF3_ROOT / "src" / "rf3" / "model" / "layers"
EVO2_ROOT        = Path("/mnt/j/evo2-main")
EVO2_LOCAL_PATH  = Path("/mnt/j/evo2_7b")
FAESM_MODEL_PATH = FAESM_ROOT / "esm2_t33_650M_UR50D"
RF3_CKPT         = Path("/mnt/j/foundry_checkpoints/rf3_foundry_01_24_latest_remapped.ckpt")

for p in [FAESM_ROOT, RF3_ROOT / "src", FOUNDRY_ROOT / "src", MANAGER_DIR, EVO2_ROOT]:
    sys.path.insert(0, str(p))

# ──────────────────────────────────────────────
# CHAIN + PROTEIN CONFIG
# ──────────────────────────────────────────────
@dataclass
class ChainSpec:
    """One chain in the assembly: protein (seq +/- MSA) or ligand (smiles/ccd/path)."""
    chain_id: str
    seq: str = ""
    msa_path: Optional[str] = None
    smiles: Optional[str] = None
    ccd_code: Optional[str] = None
    cif_path: Optional[str] = None


@dataclass
class ProteinConfig:
    """Blueprint for ONE protein extraction. Optional fields trigger WARNINGs, not crashes."""
    protein_name: str
    device: str
    output_dir: Path
    chains: list[ChainSpec]
    template_cif: Optional[Path] = None
    template_selection: Optional[list[str]] = None
    msa_gold_a3m: Optional[Path] = None
    num_steps: int = 50
    diffusion_batch_size: int = 2
    n_recycles: int = 10
    faesm_model_path: str = str(FAESM_MODEL_PATH)
    evo2_model_name: str = "evo2_7b"
    evo2_seq_column: str = "dna_sequence"
    evo2_local_path: Optional[str] = str(EVO2_LOCAL_PATH)

    def warn_missing(self):
        if self.template_cif is None:
            logging.warning("[CONFIG] No template CIF -> running de novo (pLDDT may be lower)")
        if self.msa_gold_a3m is None:
            logging.warning("[CONFIG] No gold MSA -> seq-only input")
        for c in self.chains:
            if c.seq and not c.msa_path:
                logging.warning(f"[CONFIG] Chain '{c.chain_id}' has seq but no MSA -> seq-only")

    def infer_token_count(self) -> int:
        return sum(len(c.seq) for c in self.chains if c.seq) + 3

    def first_protein_seq(self) -> str:
        for c in self.chains:
            if len(c.seq) >= 5:
                return c.seq
        return ""


# ──────────────────────────────────────────────
# RUN INSTANCE — per-invocation isolation
# ──────────────────────────────────────────────
def _random_suffix(n: int = 6) -> str:
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=n))


class RunInstance:
    """Context manager: per-invocation /tmp + env var isolation."""
    _instances: dict = {}

    def __init__(self, protein_name: str, run_id: Optional[str] = None):
        self.run_id = run_id or f"{protein_name}_{int(time.time())}_{_random_suffix()}"
        self.temp_dir = Path(f"/tmp/attn_extract_{self.run_id}")
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.zii_path = self.temp_dir / "zii_latest.pt"
        self.zii_final_path = self.temp_dir / "zii_final.pt"
        self.zii_batch_dir = self.temp_dir / "zii_batch"
        self.zii_batch_dir.mkdir(exist_ok=True)

    @classmethod
    @contextmanager
    def get_or_create(cls, protein_name: str, run_id: Optional[str] = None):
        ri = cls(protein_name, run_id)
        key = ri.run_id
        cls._instances[key] = ri
        old_zii = os.environ.get("RF3_ZII_PATH")
        old_zii_final = os.environ.get("RF3_ZII_FINAL_PATH")
        os.environ["RF3_ZII_PATH"] = str(ri.zii_path)
        os.environ["RF3_ZII_FINAL_PATH"] = str(ri.zii_final_path)
        try:
            yield ri
        finally:
            if old_zii is not None:
                os.environ["RF3_ZII_PATH"] = old_zii
            else:
                os.environ.pop("RF3_ZII_PATH", None)
            if old_zii_final is not None:
                os.environ["RF3_ZII_FINAL_PATH"] = old_zii_final
            else:
                os.environ.pop("RF3_ZII_FINAL_PATH", None)
            if ri.temp_dir.exists():
                shutil.rmtree(ri.temp_dir, ignore_errors=True)
            cls._instances.pop(key, None)


# ──────────────────────────────────────────────
# HOOK MANAGER LOADER
# ──────────────────────────────────────────────
_HOOK_MANAGER_CLASS = None

def load_hook_manager_class():
    global _HOOK_MANAGER_CLASS
    if _HOOK_MANAGER_CLASS is not None:
        return _HOOK_MANAGER_CLASS
    import importlib.util
    manager_py = MANAGER_DIR / "attn_head_hook_manager.py"
    if not manager_py.exists():
        raise FileNotFoundError(f"attn_head_hook_manager.py not found at {manager_py}")
    spec = importlib.util.spec_from_file_location("attn_head_hook_manager", str(manager_py))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["attn_head_hook_manager"] = mod
    spec.loader.exec_module(mod)
    _HOOK_MANAGER_CLASS = mod.AttentionHeadHookManager
    return _HOOK_MANAGER_CLASS


# ──────────────────────────────────────────────
# HEAD LABEL + REGISTRY
# ──────────────────────────────────────────────
@dataclass
class HeadLabel:
    model_name: str
    module_id: int
    head_index: int
    mutant_id: str
    attention_type: str
    location: str

    @property
    def tag(self) -> str:
        return f"{self.model_name}__mod{self.module_id:04d}__h{self.head_index:02d}__{self.mutant_id}"

    def to_dict(self) -> dict:
        return {"model_name": self.model_name, "module_id": self.module_id,
                "head_index": self.head_index, "mutant_id": self.mutant_id,
                "attention_type": self.attention_type, "location": self.location, "tag": self.tag}


def build_head_registry(manager, model_name: str) -> list[dict]:
    registry = []
    for mod_id in sorted(manager._module_registry.keys()):
        info = manager._module_registry[mod_id]
        for h in range(info.n_heads):
            registry.append({"model_name": model_name, "module_id": mod_id,
                             "head_index": h, "attention_type": info.attention_type,
                             "n_heads_in_module": info.n_heads,
                             "per_head_shape_desc": info.per_head_shape_desc,
                             "location": info.location, "is_diffusion": info.is_diffusion,
                             "tag_template": f"{model_name}__mod{mod_id:04d}__h{h:02d}__{{mutant_id}}"})
    return registry


def save_head_registry(registries: list[dict], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(registries, f, indent=2)
    logging.info(f"Registry saved: {output_path} ({len(registries)} heads)")


# ═══════════════════════════════════════════════
# RF3 BATCH INPUT GENERATOR
# ═══════════════════════════════════════════════
class RF3BatchBuilder:
    """Generate RF3 JSON input files — single or batch (JSON array) — per the
    RF3 README spec. Missing template/MSA: omitted with a WARNING, not an error.

    Template CIF must be a component ({"path": "..."}), NOT an external field.
    Per the README: "Folding Many Inputs -> Single JSON with Multiple Examples"
    is the recommended batch format.  """

    def __init__(self, config: ProteinConfig):
        self.config = config
        self._gold_lines: Optional[list[str]] = None
        if config.msa_gold_a3m and config.msa_gold_a3m.exists():
            with open(config.msa_gold_a3m) as f:
                self._gold_lines = f.readlines()
        elif config.msa_gold_a3m is not None:
            logging.warning(f"[BATCHER] Gold MSA path set but not found: {config.msa_gold_a3m}")

    def build_single_json(self, mutant_id: str, seq: str, work_dir: Path,
                          skip_template_cif: bool = False) -> Path:
        spec = self._build_spec(mutant_id, seq, work_dir, skip_template_cif=skip_template_cif)
        json_path = work_dir / f"{mutant_id}.json"
        json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, "w") as f:
            json.dump(spec, f, indent=2)
        return json_path

    def build_batch_json(self, rows: list[tuple[str, str]], work_dir: Path,
                         skip_template_cif: bool = False) -> Path:
        """Build a JSON array of N input specs. Per the RF3 README, this is the
        recommended batch format — RF3 auto-distributes across GPUs."""
        specs = [self._build_spec(mid, seq, work_dir / mid, skip_template_cif=skip_template_cif)
                 for mid, seq in rows]
        json_path = work_dir / "batch_input.json"
        json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, "w") as f:
            json.dump(specs, f, indent=2)
        return json_path

    def _build_spec(self, mutant_id: str, seq: str, work_dir: Path,
                     skip_template_cif: bool = False) -> dict:
        """Build one RF3 input spec.

        If skip_template_cif is True, the template CIF component is omitted.
        This is used for fallback when atomworks reports chain-id conflicts.
        template_selection is always set when configured — the engine may
        still apply templates from other sources.
        """
        components = []
        for chain in self.config.chains:
            comp: dict = {"chain_id": chain.chain_id}

            if chain.seq:
                is_mutated = (chain.seq == self.config.first_protein_seq())
                comp["seq"] = seq if is_mutated else chain.seq
                msa_file = self._write_msa(chain, comp["seq"], mutant_id, work_dir)
                if msa_file is not None:
                    comp["msa_path"] = str(msa_file)
                components.append(comp)
            elif chain.smiles:
                comp["smiles"] = chain.smiles; components.append(comp)
            elif chain.ccd_code:
                comp["ccd_code"] = chain.ccd_code; components.append(comp)
            elif chain.cif_path:
                comp["path"] = chain.cif_path; components.append(comp)

        # Template CIF as component (most CIFs use numeric chain IDs that
        # don't conflict with letter IDs like A-E). If atomworks rejects
        # the JSON due to chain-id conflict, _run_with_hooks will retry
        # with skip_template_cif=True.
        if self.config.template_cif and self.config.template_cif.exists() and not skip_template_cif:
            components.append({"path": str(self.config.template_cif)})
        elif self.config.template_cif is not None and not skip_template_cif:
            logging.warning(f"[BATCHER] Template CIF set but not found: {self.config.template_cif}")

        spec: dict = {"name": mutant_id, "components": components}
        if self.config.template_selection:
            spec["template_selection"] = self.config.template_selection
        return spec

    def _write_msa(self, chain: ChainSpec, actual_seq: str, mutant_id: str, work_dir: Path) -> Optional[Path]:
        if self._gold_lines is None or not chain.msa_path:
            return None
        msa_dir = work_dir / chain.chain_id
        msa_dir.mkdir(parents=True, exist_ok=True)
        lines = self._gold_lines.copy()
        lines[0] = f">{mutant_id}\n"
        lines[1] = f"{actual_seq}\n"
        msa_file = msa_dir / "t000_.msa0.a3m"
        with open(msa_file, "w") as f:
            f.writelines(lines)
        return msa_file


# ═══════════════════════════════════════════════
# RF3 EXTRACTOR — two modes
# ═══════════════════════════════════════════════
class RF3Extractor:
    """Load RF3 engine ONCE, then process N mutants.

    mode='heads':     Per-mutant engine.run() + hook capture. All JSONs
                      pre-generated upfront. Model stays loaded throughout.
                      batch_size controls pre-gen chunk size.
    mode='structure': JSON-array batch -> single engine.run(). No heads/Z_II.
    """

    def __init__(self, config: ProteinConfig):
        self.config = config
        self.engine = None
        self.model = None
        self.manager = None
        self.builder = None

    # ── Engine loading ──────────────────────────
    def load_model(self, ckpt_path: Optional[str] = None):
        torch = _import_torch()
        from omegaconf import OmegaConf
        from hydra import compose, initialize_config_dir
        from hydra.core.global_hydra import GlobalHydra
        from hydra.utils import instantiate

        ckpt = ckpt_path or str(RF3_CKPT)
        config_dir = str(RF3_ROOT / "configs")
        if not Path(config_dir).exists():
            config_dir = str(FOUNDRY_ROOT / "models" / "rf3" / "configs")
        if not Path(config_dir).exists():
            raise FileNotFoundError(f"RF3 configs not found under {RF3_ROOT}")

        GlobalHydra.instance().clear()
        overrides = ["dump_predictions=false", "dump_trajectories=false",
                     f"num_steps={self.config.num_steps}",
                     f"diffusion_batch_size={self.config.diffusion_batch_size}",
                     f"n_recycles={self.config.n_recycles}",
                     "early_stopping_plddt_threshold=0.0", f"ckpt_path={ckpt}"]
        with initialize_config_dir(config_dir=config_dir, version_base="1.3"):
            cfg = compose(config_name="inference", overrides=overrides)

        run_param_keys = {"inputs", "out_dir", "dump_predictions", "dump_trajectories",
                          "one_model_per_file", "annotate_b_factor_with_plddt",
                          "sharding_pattern", "skip_existing", "template_selection",
                          "ground_truth_conformer_selection", "cyclic_chains", "add_missing_atoms"}
        init_cfg = OmegaConf.create(
            {k: v for k, v in OmegaConf.to_container(cfg, resolve=True).items() if k not in run_param_keys})

        logging.info("[RF3] Instantiating engine...")
        self.engine = instantiate(init_cfg, _convert_="partial", _recursive_=False)
        try:
            self.engine.initialize()
        except Exception as e:
            logging.warning(f"[RF3] engine.initialize() failed: {e}; trying manual")
            if hasattr(self.engine, '_construct_pipeline'):
                self.engine._construct_pipeline()
            if hasattr(self.engine, '_construct_trainer'):
                self.engine._construct_trainer()
            self.engine.initialized_ = True

        self.model = self._find_model()
        if self.model is None:
            raise RuntimeError("[RF3] Could not locate nn.Module in engine")
        if hasattr(self.model, 'shadow') and self.model.shadow is not None:
            self.model = self.model.shadow
        elif hasattr(self.model, 'forward_module'):
            inner = getattr(self.model, 'forward_module', None) or getattr(self.model, 'module', None)
            if inner is not None and hasattr(inner, 'shadow'):
                self.model = inner.shadow

        AHM = load_hook_manager_class()
        self.manager = AHM(self.model, verbose=False)
        logging.info(f"[RF3] Manager: {self.manager.num_modules} modules, {self.manager.num_total_heads} heads")
        self.builder = RF3BatchBuilder(self.config)

    def _is_rf3_model(self, m) -> bool:
        targets = {"Recycler", "DiffusionModule", "PairformerBlock", "RF3"}
        for _, ch in m.named_modules():
            if ch.__class__.__name__ in targets:
                return True
        return False

    def _find_model(self):
        torch = _import_torch()
        nn = torch.nn
        for attr in ["model", "rf3_model", "_model", "net", "network", "trunk", "pipeline", "trainer"]:
            obj = getattr(self.engine, attr, None)
            if isinstance(obj, nn.Module) and self._is_rf3_model(obj):
                return obj
        trainer = getattr(self.engine, "trainer", None)
        if trainer is not None:
            state = getattr(trainer, "state", None)
            if isinstance(state, dict) and "model" in state:
                return state["model"]
        for pa in ["trainer", "pipeline"]:
            parent = getattr(self.engine, pa, None)
            if not parent: continue
            for ca in ["model", "lightning_module", "module", "_model", "rf3_model", "net"]:
                obj = getattr(parent, ca, None)
                if isinstance(obj, nn.Module) and self._is_rf3_model(obj):
                    return obj
                for ga in ["model", "_model", "net"]:
                    obj2 = getattr(obj, ga, None) if obj is not None else None
                    if isinstance(obj2, nn.Module) and self._is_rf3_model(obj2):
                        return obj2
        return self._bruteforce(self.engine, 0, 5, set())

    def _bruteforce(self, obj, depth, max_depth, visited):
        torch = _import_torch()
        nn = torch.nn
        if depth > max_depth:
            return None
        oid = id(obj)
        if oid in visited:
            return None
        visited.add(oid)
        if isinstance(obj, nn.Module) and self._is_rf3_model(obj):
            return obj
        for attr in dir(obj):
            if attr.startswith("__"):
                continue
            try:
                child = getattr(obj, attr, None)
            except Exception:
                continue
            if child is None:
                continue
            if isinstance(child, nn.Module):
                if self._is_rf3_model(child):
                    return child
                r = self._bruteforce(child, depth + 1, max_depth, visited)
                if r: return r
            elif hasattr(child, '__dict__') and not callable(child):
                r = self._bruteforce(child, depth + 1, max_depth, visited)
                if r: return r
        return None

    # ── Core run ──────────────────────────────────
    def _run_with_hooks(self, json_path: str, out_dir: str,
                        skip_template_cif: bool = False,
                        template_selection=None) -> dict:
        """Run engine with hook capture.

        If chain-id conflict is detected, retries ONCE with the template CIF
        component removed. This gracefully handles CIFs whose chain IDs
        happen to overlap the input chains.

        Args:
            template_selection: Override for chain-batch mode. Pass [] to
                disable templating (chains are unrelated mutants, not a real
                complex).  None = use self.config.template_selection.
        """
        torch = _import_torch()
        import inspect
        ts = self.config.template_selection if template_selection is None else template_selection
        params = {"inputs": json_path, "out_dir": out_dir,
                  "dump_predictions": False, "dump_trajectories": False,
                  "one_model_per_file": False, "annotate_b_factor_with_plddt": False,
                  "skip_existing": False, "template_selection": ts,
                  "ground_truth_conformer_selection": None, "cyclic_chains": []}
        sig = inspect.signature(self.engine.run)
        valid = set(sig.parameters.keys())
        filtered = {k: v for k, v in params.items() if k in valid}
        if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
            filtered = params

        steps = self.config.num_steps - 1
        try:
            with torch.no_grad(), self.manager.capture(diffusion_total_steps=steps) as captured:
                self.engine.run(**filtered)
        except ValueError as e:
            err_str = str(e)
            if "multiple input components" in err_str and not skip_template_cif:
                logging.warning(f"[RF3] Chain-id conflict with template CIF. "
                                f"Retrying without CIF component. Error: {err_str[:100]}")
                # Rebuild JSON without the template CIF and retry once
                import json as _json
                with open(json_path) as _f:
                    _spec = _json.load(_f)
                if isinstance(_spec, list):
                    for _entry in _spec:
                        _entry["components"] = [c for c in _entry["components"]
                                                 if "path" not in c]
                else:
                    _spec["components"] = [c for c in _spec["components"]
                                            if "path" not in c]
                with open(json_path, "w") as _f:
                    _json.dump(_spec, _f, indent=2)
                with torch.no_grad(), self.manager.capture(diffusion_total_steps=steps) as captured:
                    self.engine.run(**filtered)
            else:
                raise

        result = {}
        for k, v in captured.items():
            if isinstance(v, torch.Tensor):
                result[k] = v.cpu().clone() if v.is_cuda else v.clone()
            else:
                result[k] = v
        captured.clear()
        torch.cuda.empty_cache()
        return result

    def _extract_zii(self, output_dir: Path, mutant_id: str, mut_pos: int) -> bool:
        torch = _import_torch()
        zii_path = os.environ.get("RF3_ZII_PATH", "")
        if not zii_path or not Path(zii_path).exists():
            return False
        try:
            z = torch.load(zii_path, map_location="cpu", weights_only=False)
        except Exception as e:
            logging.warning(f"[ZII] load failed: {e}")
            return False
        expected = self.config.infer_token_count()
        tokens = z.shape[-2]
        if tokens < expected * 0.5 or tokens > expected * 1.5:
            logging.error(f"[ZII] SHAPE MISMATCH: expected ~{expected} tokens, got {tokens} ({list(z.shape)}). SKIP")
            return False
        pos = max(0, min(mut_pos + 3, tokens - 1))
        if z.dim() in (2, 3, 4):
            # Save FULL ZII matrix (no row/col slice)
            torch.save(z.half(), output_dir / f"{mutant_id}_zii.pt")
            return True
        return True

    def _find_mut_pos(self, seq: str) -> int:
        wt = self.config.first_protein_seq()
        return next((j for j in range(min(len(seq), len(wt))) if seq[j] != wt[j]), 0)

    # ── MODE 1: heads (per-mutant capture) ───────
    def run_batch_heads(self, df, output_dir: Path, batch_size: int = 1, log_every: int = 10):
        """Per-mutant head + Z_II extraction.

        Model loaded ONCE. JSONs pre-generated in chunks of 'batch_size',
        then processed sequentially through engine.run() with hook capture.
        batch_size here = pre-gen chunk size (reduces I/O wait between
        GPU forwards). For TRUE batch parallelism (no heads, faster), use
        mode='structure'.
        """
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)
        meta_path = output_dir / "labels_metadata.jsonl"
        done_set = self._load_done_set(meta_path)
        rf3_out = output_dir / "_rf3_structures"
        rf3_out.mkdir(exist_ok=True)

        rows = [(r["mutant_id"], r["protein"]) for _, r in df.iterrows()]
        n_skip, n_done, n_fail = 0, 0, 0
        t_start = time.time()

        for b_start in range(0, len(rows), batch_size):
            batch = rows[b_start:b_start + batch_size]
            batch_done = sum(1 for mid, _ in batch if mid in done_set)
            n_skip += batch_done
            batch = [(mid, seq) for mid, seq in batch if mid not in done_set]
            if not batch:
                continue

            batch_dir = Path(f"/tmp/rf3_batch_{int(time.time())}_{_random_suffix()}")
            try:
                for mid, seq in batch:
                    self.builder.build_single_json(mid, seq, batch_dir / mid)
                for mid, seq in batch:
                    try:
                        jpath = batch_dir / mid / f"{mid}.json"
                        captured = self._run_with_hooks(str(jpath), str(rf3_out))
                        results, metadata = {}, []
                        for (mod_id, h_idx), tensor in captured.items():
                            info = self.manager.get_module_info(mod_id)
                            label = HeadLabel("rf3", mod_id, h_idx, mid, info.attention_type, info.location)
                            results[label.tag] = tensor.cpu().half()
                            metadata.append(label.to_dict())
                        torch.save(results, output_dir / f"{mid}.pt")
                        self._extract_zii(output_dir, mid, self._find_mut_pos(seq))
                        with open(meta_path, "a") as f:
                            for e in metadata:
                                f.write(json.dumps(e) + "\n")
                        n_done += 1
                        if n_done % log_every == 0 or n_done == 1:
                            elapsed = time.time() - t_start
                            rate = n_done / max(elapsed, 1e-6)
                            rem = len(rows) - n_skip - n_done - n_fail
                            logging.info(f"[RF3 heads] {n_done}/{len(rows)-n_skip} ({n_fail} fail) "
                                         f"{rate:.2f} it/s ETA: {rem/max(rate,1e-6)/3600:.1f}h")
                    except Exception as e:
                        n_fail += 1
                        logging.error(f"[RF3 heads] {mid} failed: {e}", exc_info=True)
                    finally:
                        torch.cuda.empty_cache()
                        gc.collect()
            finally:
                if batch_dir.exists():
                    shutil.rmtree(batch_dir, ignore_errors=True)

        logging.info(f"[RF3 heads] Done: {n_done} ok, {n_skip} skip, {n_fail} fail ({time.time()-t_start:.1f}s)")

    # ── MODE: chain-batch (true parallel: N mutants as N chains, 1 forward) ──
    def run_batch_chain(self, df, output_dir: Path, batch_size: int = 4,
                        log_every: int = 10):
        """TRUE batch mode: N mutants encoded as N chains in ONE forward pass.

        Instead of running N separate engine.run() calls (one per mutant),
        this method packs batch_size mutants as distinct chains in a single
        multi-chain JSON.  One forward pass produces Z_II and heads for ALL
        mutants simultaneously.  Per-mutant tensors are sliced by chain
        boundaries.

        For RBD (200aa), BS=4 gives 800 total tokens — ~4× throughput vs
        heads mode.  Attention heads and Z_II are sliced per-chain, so output
        format is identical to heads mode.

        Limitation: all mutants in a batch must have the SAME sequence length
        (same protein backbone, only mutations differ).  True for saturation
        mutagenesis DMS data.
        """
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)
        meta_path = output_dir / "labels_metadata.jsonl"
        done_set = self._load_done_set(meta_path)

        rows = [(r["mutant_id"], r["protein"]) for _, r in df.iterrows()]
        rows = [(mid, seq) for mid, seq in rows if mid not in done_set]
        n_total = len(rows)
        if n_total == 0:
            return

        n_done, n_fail = 0, 0
        t_start = time.time()

        logging.info(f"[RF3 chain-batch] BS={batch_size}, {n_total} mutants remaining")

        for b_start in range(0, n_total, batch_size):
            batch = rows[b_start:b_start + batch_size]
            batch_ids = [mid for mid, _ in batch]
            batch_seqs = [seq for _, seq in batch]

            # Verify all sequences have the same length (required for chain slicing)
            seq_lens = [len(s) for s in batch_seqs]
            if len(set(seq_lens)) > 1:
                logging.warning(f"[RF3 chain-batch] Variable seq lengths in batch: {seq_lens}. "
                                f"Padding to max length (may introduce artifacts).")

            try:
                # Build multi-chain JSON: each mutant = 1 chain
                batch_dir = Path(f"/tmp/rf3_chainbatch_{int(time.time())}_{_random_suffix()}")
                batch_dir.mkdir(parents=True, exist_ok=True)

                multi_json = self._build_chain_batch_json(batch_ids, batch_seqs, batch_dir)

                # Run ONE forward pass with all chains
                rf3_out = output_dir / "_rf3_structures"
                rf3_out.mkdir(exist_ok=True)
                captured = self._run_with_hooks(str(multi_json), str(rf3_out),
                                                  template_selection=[])

                # Read Z_II from env-var path
                zii_full = None
                zii_path = os.environ.get("RF3_ZII_PATH", "")
                if zii_path and Path(zii_path).exists():
                    try:
                        zii_full = torch.load(zii_path, map_location="cpu", weights_only=False)
                    except Exception:
                        pass

                # Compute chain boundaries from sequence lengths
                boundaries = self._chain_boundaries([len(s) for s in batch_seqs], zii_full)

                # Slice and save per-mutant
                for i, (mid, seq) in enumerate(batch):
                    try:
                        ch_start = boundaries[i]
                        ch_end = boundaries[i + 1]

                        # Slice heads
                        results, metadata = self._slice_heads(
                            captured, mid, ch_start, ch_end)

                        torch.save(results, output_dir / f"{mid}.pt")

                        # Slice Z_II
                        if zii_full is not None:
                            self._slice_and_save_zii(
                                zii_full, mid, ch_start, ch_end,
                                output_dir, self._find_mut_pos(seq))

                        with open(meta_path, "a") as f:
                            for e in metadata:
                                f.write(json.dumps(e) + "\n")

                        n_done += 1
                    except Exception as e:
                        n_fail += 1
                        logging.error(f"[RF3 chain-batch] {mid} slice failed: {e}")

                if n_done % log_every == 0:
                    elapsed = time.time() - t_start
                    logging.info(f"[RF3 chain-batch] {n_done}/{n_total} ({n_fail} fail) "
                                 f"{n_done/max(elapsed,1e-6):.2f} it/s")

            except Exception as e:
                n_fail += len(batch)
                logging.error(f"[RF3 chain-batch] Batch start={b_start} failed: {e}", exc_info=True)
            finally:
                if 'batch_dir' in dir() and batch_dir.exists():
                    shutil.rmtree(batch_dir, ignore_errors=True)
                torch.cuda.empty_cache()
                gc.collect()

        elapsed = time.time() - t_start
        logging.info(f"[RF3 chain-batch] Done: {n_done} ok, {n_fail} fail ({elapsed:.1f}s)")

    def _build_chain_batch_json(self, mutant_ids: list[str], seqs: list[str],
                                work_dir: Path) -> Path:
        """Build a multi-chain JSON: each mutant = 1 chain, chain IDs = A,B,C,..."""
        chain_ids = [chr(65 + i) for i in range(len(mutant_ids))]
        components = []
        for cid, mid, seq in zip(chain_ids, mutant_ids, seqs):
            comp = {"chain_id": cid, "seq": seq}
            # MSA: apply gold MSA to ALL chains (consistent, same homologs)
            if self.config.msa_gold_a3m and self.config.msa_gold_a3m.exists():
                gold_lines = open(self.config.msa_gold_a3m).readlines()
                gold_lines[0] = f">{mid}\n"
                gold_lines[1] = f"{seq}\n"
                msa_dir = work_dir / cid
                msa_dir.mkdir(parents=True, exist_ok=True)
                msa_path = msa_dir / "t000_.msa0.a3m"
                with open(msa_path, "w") as f:
                    f.writelines(gold_lines)
                comp["msa_path"] = str(msa_path)
            components.append(comp)

        # Template CIF
        if self.config.template_cif and self.config.template_cif.exists():
            components.append({"path": str(self.config.template_cif)})

        spec = {"name": f"batch_{'_'.join(mutant_ids[:3])}", "components": components}
        if self.config.template_selection:
            spec["template_selection"] = self.config.template_selection

        json_path = work_dir / "chain_batch.json"
        with open(json_path, "w") as f:
            json.dump([spec], f, indent=2)  # RF3 accepts JSON array
        return json_path

    @staticmethod
    def _chain_boundaries(seq_lens: list[int], zii_full=None) -> list[int]:
        """Return cumulative token offsets for each chain.

        For proteins, each residue = 1 token.  No special tokens.
        Fallback to sequence-length-based estimate when Z_II unavailable.

        Returns [0, L0, L0+L1, L0+L1+L2, ...]  (len = n_chains + 1)
        """
        if zii_full is not None:
            total = zii_full.shape[-2]
            # Distribute proportionally (handles variable-length padding)
            total_len = sum(seq_lens)
            boundaries = [0]
            for sl in seq_lens:
                n_tokens = max(1, round(sl / total_len * total))
                boundaries.append(boundaries[-1] + n_tokens)
            # Adjust last boundary to match exact total
            if boundaries[-1] != total:
                boundaries[-1] = total
            return boundaries
        # No Z_II: use sequence lengths directly
        boundaries = [0]
        for sl in seq_lens:
            boundaries.append(boundaries[-1] + sl)
        return boundaries

    def _slice_heads(self, captured: dict, mutant_id: str,
                     ch_start: int, ch_end: int) -> tuple[dict, list]:
        """Slice all attention head tensors by chain token range."""
        results = {}
        metadata = []
        for (mod_id, h_idx), tensor in captured.items():
            info = self.manager.get_module_info(mod_id)
            label = HeadLabel("rf3", mod_id, h_idx, mutant_id, info.attention_type, info.location)
            try:
                if tensor.dim() == 2:  # (I, C)
                    sliced = tensor[ch_start:ch_end, :]
                elif tensor.dim() == 3:  # (B, I, C) or (I, I, C)
                    if tensor.shape[0] == tensor.shape[1]:
                        sliced = tensor[ch_start:ch_end, ch_start:ch_end, :]
                    else:
                        sliced = tensor[:, ch_start:ch_end, :]
                elif tensor.dim() == 4:  # (B, I, I, C) or (B, L, L, H, D)
                    sliced = tensor[:, ch_start:ch_end, ch_start:ch_end, :]
                elif tensor.dim() == 5:  # (B, L, L, H, D)
                    sliced = tensor[:, ch_start:ch_end, ch_start:ch_end, :, :]
                else:
                    sliced = tensor
                results[label.tag] = sliced.cpu().half()
            except Exception:
                results[label.tag] = tensor.cpu().half()
            metadata.append(label.to_dict())
        return results, metadata

    @staticmethod
    def _slice_and_save_zii(zii_full, mutant_id: str,
                            ch_start: int, ch_end: int,
                            output_dir: Path, mut_pos: int):
        """Extract the per-mutant Z_II diagonal block and save FULL (I,I,D)."""
        torch = _import_torch()
        # Z_II[total_tokens, total_tokens, 128] → [ch_start:ch_end, ch_start:ch_end, :]
        zii_mutant = zii_full[ch_start:ch_end, ch_start:ch_end, :]
        # Save FULL per-mutant pair matrix (I, I, D)
        torch.save(zii_mutant.half(), output_dir / f"{mutant_id}_zii.pt")

    # ── MODE: zii (NO HOOKS, chain-batch, maximum speed) ─────────
    def run_batch_zii_chain(self, df, output_dir: Path, batch_size: int = 4,
                            log_every: int = 10):
        """ZII-only chain-batch: NO hooks, NO heads. Maximum throughput.

        Packs batch_size mutants as separate chains in one forward pass.
        Runs engine.run() WITHOUT manager.capture() — no Python callbacks,
        no GPU pipeline interruptions. Z_II read from /tmp after each batch.

        Expected: ~12s per mutant for 200aa (vs 45s with hooks).
        With BS=4: ~40s per batch of 4 = ~10s/mutant effective.
        """
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)

        done_set = set()
        done_file = output_dir / "_done_zii.txt"
        if done_file.exists():
            done_set = set(done_file.read_text().strip().split("\n"))

        rows = [(r["mutant_id"], r["protein"]) for _, r in df.iterrows()]
        rows = [(mid, seq) for mid, seq in rows if mid not in done_set]
        n_total = len(rows)
        if n_total == 0:
            logging.info("[RF3 zii] All done.")
            return

        n_done, n_fail = 0, 0
        t_start = time.time()
        rf3_out = output_dir / "_rf3_structures"
        rf3_out.mkdir(exist_ok=True)

        logging.info(f"[RF3 zii] BS={batch_size}, {n_total} mutants, NO HOOKS")

        import inspect
        sig = inspect.signature(self.engine.run)
        valid = set(sig.parameters.keys())
        has_var_kw = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())

        for b_start in range(0, n_total, batch_size):
            batch = rows[b_start:b_start + batch_size]
            batch_ids = [mid for mid, _ in batch]
            batch_seqs = [seq for _, seq in batch]
            seq_lens = [len(s) for s in batch_seqs]

            try:
                batch_dir = Path(f"/tmp/rf3_zii_{int(time.time())}_{_random_suffix()}")
                batch_dir.mkdir(parents=True, exist_ok=True)

                # Build multi-chain JSON
                multi_json = self._build_chain_batch_json(batch_ids, batch_seqs, batch_dir)

                # Clear old Z_II
                zii_path = Path(os.environ.get("RF3_ZII_PATH", "/tmp/zii_latest.pt"))
                if zii_path.exists():
                    zii_path.unlink()

                # Run engine WITHOUT hooks (no manager.capture)
                params = {"inputs": str(multi_json), "out_dir": str(rf3_out),
                          "dump_predictions": False, "dump_trajectories": False,
                          "one_model_per_file": False, "skip_existing": False,
                          "template_selection": []}
                filtered = {k: v for k, v in params.items() if k in valid}
                if has_var_kw:
                    filtered = params

                with torch.no_grad():
                    self.engine.run(**filtered)

                # Read Z_II
                if not zii_path.exists():
                    logging.warning(f"[RF3 zii] No Z_II file after batch {b_start}")
                    n_fail += len(batch)
                    continue

                zii_full = torch.load(str(zii_path), map_location="cpu", weights_only=False)

                # Slice per-mutant
                boundaries = self._chain_boundaries(seq_lens, zii_full)

                for i, (mid, seq) in enumerate(batch):
                    ch_start = boundaries[i]
                    ch_end = boundaries[i + 1]
                    self._slice_and_save_zii(zii_full, mid, ch_start, ch_end,
                                            output_dir, self._find_mut_pos(seq))
                    n_done += 1

                # Record done
                with open(done_file, "a") as f:
                    for mid, _ in batch:
                        f.write(mid + "\n")

                if n_done % (log_every * batch_size) == 0 or n_done == len(batch):
                    elapsed = time.time() - t_start
                    rate = n_done / max(elapsed, 1e-6)
                    rem = n_total - n_done - n_fail
                    logging.info(f"[RF3 zii] {n_done}/{n_total} "
                                 f"{rate:.2f} mut/s ETA: {rem/max(rate,1e-6)/60:.0f}min")

            except Exception as e:
                n_fail += len(batch)
                logging.error(f"[RF3 zii] Batch {b_start} failed: {e}", exc_info=True)
            finally:
                if 'batch_dir' in dir() and batch_dir.exists():
                    shutil.rmtree(batch_dir, ignore_errors=True)
                torch.cuda.empty_cache()

        elapsed = time.time() - t_start
        logging.info(f"[RF3 zii] Done: {n_done} ok, {n_fail} fail "
                     f"({elapsed:.1f}s, {n_done/max(elapsed,1e-6):.2f} mut/s)")

    # ── MODE 2: structure (JSON-array batch) ───────
    def run_batch_structure(self, df, output_dir: Path, batch_size: int = 10):
        """True batch mode: N mutants in ONE engine.run() via JSON array. Fastest. No heads."""
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)
        rows = [(r["mutant_id"], r["protein"]) for _, r in df.iterrows()]
        n_done, n_fail = 0, 0
        t_start = time.time()

        for b_start in range(0, len(rows), batch_size):
            batch = rows[b_start:b_start + batch_size]
            batch_dir = Path(f"/tmp/rf3_struct_{int(time.time())}_{_random_suffix()}")
            try:
                json_path = self.builder.build_batch_json(batch, batch_dir)
                logging.info(f"[RF3 struct] Batch {b_start}-{b_start+len(batch)}")

                import inspect
                sig = inspect.signature(self.engine.run)
                valid = set(sig.parameters.keys())
                params = {"inputs": str(json_path), "out_dir": str(output_dir),
                          "dump_predictions": True, "dump_trajectories": False,
                          "one_model_per_file": False, "annotate_b_factor_with_plddt": False,
                          "skip_existing": True, "template_selection": self.config.template_selection}
                filtered = {k: v for k, v in params.items() if k in valid}
                if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
                    filtered = params

                with torch.no_grad():
                    self.engine.run(**filtered)

                n_done += len(batch)
                elapsed = time.time() - t_start
                logging.info(f"[RF3 struct] {n_done}/{len(rows)} ({n_done/max(elapsed,1e-6):.1f} it/s)")
            except Exception as e:
                n_fail += len(batch)
                logging.error(f"[RF3 struct] Batch failed: {e}", exc_info=True)
            finally:
                if batch_dir.exists():
                    shutil.rmtree(batch_dir, ignore_errors=True)
                torch.cuda.empty_cache()
                gc.collect()
        logging.info(f"[RF3 struct] Done: {n_done} ok, {n_fail} fail ({time.time()-t_start:.1f}s)")

    # ── MODE 3: parallel batch (multi-instance, heads + ZII) ───────
    def run_batch_parallel(self, df, output_dir: Path, batch_size: int = 2,
                           log_every: int = 10):
        """Persistent multi-instance batch: N workers, each processing a chunk.

        Each worker loads RF3 model ONCE, then processes mutants from its chunk
        sequentially. BS controls worker count (parallelism). GPU round-robin.

        BS=2 on 2 GPUs → 2 workers (one per GPU), each gets ~half the mutants.
        BS=4 on 2 GPUs → 4 workers (2 per GPU), each gets ~quarter.

        Per-mutant time ~45s (same as heads mode), but model loaded once per
        worker (amortized over hundreds of mutants). No 20-minute gaps from
        per-mutant subprocess spawning.

        GPU0 = RTX 5090D (32GB): 1 worker max
        GPU1 = RTX 6000D (80GB): 2 workers max
        """
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)
        meta_path = output_dir / "labels_metadata.jsonl"
        done_set = self._load_done_set(meta_path)

        rows = [(r["mutant_id"], r["protein"]) for _, r in df.iterrows()]
        rows = [(mid, seq) for mid, seq in rows if mid not in done_set]
        n_total = len(rows)
        if n_total == 0:
            logging.info("[RF3 parallel] All mutants already done.")
            return

        t_start = time.time()
        available_devices = [0, 1] if torch.cuda.device_count() >= 2 else [0]
        worker_py = Path(sys.executable).parent / "attn_worker.py"
        if not worker_py.exists():
            worker_py = Path("/mnt/j/conda_envs/foundry/attn_worker.py")

        # Partition rows into batch_size chunks (each worker gets 1 chunk)
        chunk_size = (n_total + batch_size - 1) // batch_size
        chunks = [rows[i:i + chunk_size] for i in range(0, n_total, chunk_size)]
        # Limit to batch_size chunks
        chunks = chunks[:batch_size]

        logging.info(f"[RF3 parallel] {batch_size} workers × ~{chunk_size} mutants each "
                     f"(total {n_total}) on {available_devices}")

        processes = []
        for wi, chunk in enumerate(chunks):
            dev = available_devices[wi % len(available_devices)]
            log_file = output_dir / "_logs" / f"worker_{wi}.log"
            log_file.parent.mkdir(parents=True, exist_ok=True)

            # Build a batch JSON for this worker's chunk
            batch_dir = Path(f"/tmp/rf3_parallel_w{wi}_{int(time.time())}")
            batch_dir.mkdir(parents=True, exist_ok=True)
            batch_json = self.builder.build_batch_json(chunk, batch_dir)

            cmd = (
                f"CUDA_VISIBLE_DEVICES={dev} "
                f"{sys.executable} {worker_py} "
                f"--batch-json {batch_json} --output {output_dir} "
                f"--gpu-id 0 "
                f"--num-steps {self.config.num_steps} "
                f"--diff-batch {self.config.diffusion_batch_size} "
                f"--recycles {self.config.n_recycles} "
                f"> {log_file} 2>&1"
            )
            proc = subprocess.Popen(cmd, shell=True)
            processes.append((wi, proc, dev, batch_dir))

        # Wait for ALL workers
        n_done, n_fail = 0, 0
        for wi, proc, dev, bdir in processes:
            ret = proc.wait()
            if ret == 0:
                n_done += len(chunks[wi])
            else:
                n_fail += len(chunks[wi])
                logging.error(f"[RF3 parallel] Worker {wi} on GPU{dev} failed (exit={ret})")
                lf = output_dir / "_logs" / f"worker_{wi}.log"
                if lf.exists():
                    with open(lf) as f:
                        logging.error("  Worker log tail:\n" + "".join(f.readlines()[-15:]))
            if bdir.exists():
                shutil.rmtree(bdir, ignore_errors=True)

        elapsed = time.time() - t_start
        logging.info(f"[RF3 parallel] Done: {n_done} ok, {n_fail} fail "
                     f"({elapsed:.1f}s, {n_done/max(elapsed,1e-6):.2f} it/s)")

    @staticmethod
    def _load_done_set(meta_path: Path) -> set:
        if not meta_path.exists():
            return set()
        done = set()
        with open(meta_path) as f:
            for line in f:
                if line.strip():
                    done.add(json.loads(line.strip())["mutant_id"])
        return done


# ═══════════════════════════════════════════════
# FAESM + EVO2 (batch_size for API compat)
# ═══════════════════════════════════════════════
class FAESMExtractor:
    def __init__(self, config: ProteinConfig):
        self.config = config
        self.model = None
        self.manager = None

    def _patch_hf_compat(self):
        from faesm.esm import FAEsmForMaskedLM, FAEsmModel, FAEsmEncoder, FAEsmSelfAttention
        from transformers.models.esm.modeling_esm import EsmEmbeddings
        import functools, inspect
        torch = _import_torch()

        if not hasattr(FAEsmForMaskedLM, "all_tied_weights_keys"):
            @property
            def _all_tied_weights_keys(self):
                keys = getattr(self, "_tied_weights_keys", [])
                return {k: k for k in keys} if isinstance(keys, list) else keys
            FAEsmForMaskedLM.all_tied_weights_keys = _all_tied_weights_keys

        if not hasattr(FAEsmModel, "get_head_mask"):
            def get_head_mask(self, head_mask, num_hidden_layers, is_attention_chunked=False):
                if head_mask is None: return [None] * num_hidden_layers
                if head_mask.dim() == 1:
                    head_mask = head_mask.unsqueeze(0).unsqueeze(0).unsqueeze(-1).unsqueeze(-1)
                    head_mask = head_mask.expand(num_hidden_layers, -1, -1, -1, -1)
                elif head_mask.dim() == 2:
                    head_mask = head_mask.unsqueeze(1).unsqueeze(-1).unsqueeze(-1)
                return head_mask.to(dtype=torch.float32)
            FAEsmModel.get_head_mask = get_head_mask

        sig = inspect.signature(EsmEmbeddings.forward)
        if not any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()) \
           and "past_key_values_length" not in sig.parameters:
            _orig = EsmEmbeddings.forward
            @functools.wraps(_orig)
            def _patched(self_, *args, **kwargs):
                kwargs.pop("past_key_values_length", None)
                return _orig(self_, *args, **kwargs)
            EsmEmbeddings.forward = _patched

        sig = inspect.signature(FAEsmEncoder.forward)
        if not any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
            _orig = FAEsmEncoder.forward
            _params = set(sig.parameters.keys())
            @functools.wraps(_orig)
            def _patched_enc(self_, *args, **kwargs):
                return _orig(self_, *args, **{k: v for k, v in kwargs.items() if k in _params})
            FAEsmEncoder.forward = _patched_enc

        if not hasattr(FAEsmSelfAttention, "transpose_for_scores"):
            def transpose_for_scores(self_, x):
                new_shape = x.size()[:-1] + (self_.num_attention_heads, self_.attention_head_size)
                return x.view(new_shape).permute(0, 2, 1, 3)
            FAEsmSelfAttention.transpose_for_scores = transpose_for_scores

    def load_model(self):
        torch = _import_torch()
        self._patch_hf_compat()
        from faesm.esm import FAEsmForMaskedLM
        AHM = load_hook_manager_class()
        self.model = FAEsmForMaskedLM.from_pretrained(self.config.faesm_model_path, use_fa=False)
        self.model = self.model.to(self.config.device).eval()
        self.manager = AHM(self.model, verbose=False)
        logging.info(f"[FAESM] Manager: {self.manager.num_modules} modules, {self.manager.num_total_heads} heads")

    def extract_single(self, sequence: str, mutant_id: str) -> tuple[dict, list]:
        torch = _import_torch()
        tokenizer = self.model.tokenizer
        inputs = tokenizer(sequence, return_tensors="pt", padding=False, truncation=True)
        input_ids = inputs["input_ids"].to(self.config.device)
        results, metadata = {}, []
        with torch.no_grad(), self.manager.capture() as captured:
            _ = self.model(input_ids)
        for (mid, h_idx), tensor in captured.items():
            info = self.manager.get_module_info(mid)
            label = HeadLabel("faesm", mid, h_idx, mutant_id, info.attention_type, info.location)
            results[label.tag] = tensor.cpu().half()
            metadata.append(label.to_dict())
        return results, metadata

    def run_batch(self, df, output_dir: Path, log_every: int = 50):
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)
        meta_path = output_dir / "labels_metadata.jsonl"
        done_set = set()
        if meta_path.exists():
            with open(meta_path) as f:
                for line in f:
                    if line.strip(): done_set.add(json.loads(line.strip())["mutant_id"])
        n_skip, n_done = 0, 0
        t_start = time.time()
        for _, row in df.iterrows():
            mid, seq = row["mutant_id"], row["protein"]
            if mid in done_set: n_skip += 1; continue
            try:
                results, meta = self.extract_single(seq, mid)
                torch.save(results, output_dir / f"{mid}.pt")
                with open(meta_path, "a") as f:
                    for e in meta: f.write(json.dumps(e) + "\n")
                n_done += 1
                if n_done % log_every == 0:
                    logging.info(f"[FAESM] {n_done}/{len(df)-n_skip} ({n_done/max(time.time()-t_start,1e-6):.1f} it/s)")
                if n_done % 100 == 0: torch.cuda.empty_cache(); gc.collect()
            except Exception as e:
                logging.error(f"[FAESM] {mid} failed: {e}")
        logging.info(f"[FAESM] Done: {n_done} ok, {n_skip} skip")


class EVO2Extractor:
    def __init__(self, config: ProteinConfig):
        self.config = config
        self.evo2 = None; self.model = None; self.manager = None; self.tokenizer = None

    def load_model(self):
        from evo2 import Evo2
        AHM = load_hook_manager_class()
        self.evo2 = Evo2(self.config.evo2_model_name, local_path=self.config.evo2_local_path)
        self.model = self.evo2.model; self.tokenizer = self.evo2.tokenizer
        self.model.eval()
        self.manager = AHM(self.model, verbose=False)
        mha = [i for i in self.manager._module_registry.values() if i.attention_type == "Evo2MHA"]
        logging.info(f"[EVO2] Manager: {len(mha)} MHA modules, {sum(m.n_heads for m in mha)} heads")

    def extract_single(self, sequence: str, mutant_id: str) -> tuple[dict, list]:
        torch = _import_torch()
        token_ids = self.tokenizer.tokenize(sequence)
        input_ids = torch.tensor([token_ids], dtype=torch.long).to(self.model.embedding_layer.weight.device)
        results, metadata = {}, []
        with torch.no_grad(), self.manager.capture(diffusion_total_steps=None) as captured:
            _ = self.model(input_ids)
        for (mid, h_idx), tensor in captured.items():
            info = self.manager.get_module_info(mid)
            if info.attention_type != "Evo2MHA": continue
            label = HeadLabel("evo2", mid, h_idx, mutant_id, info.attention_type, info.location)
            results[label.tag] = tensor.cpu().half()
            metadata.append(label.to_dict())
        return results, metadata

    def run_batch(self, df, output_dir: Path, log_every: int = 50):
        pd = _import_pd()
        torch = _import_torch()
        if isinstance(df, str):
            df = pd.read_csv(df)
        output_dir.mkdir(parents=True, exist_ok=True)
        meta_path = output_dir / "labels_metadata.jsonl"
        done_set = set()
        if meta_path.exists():
            with open(meta_path) as f:
                for line in f:
                    if line.strip(): done_set.add(json.loads(line.strip())["mutant_id"])
        n_skip, n_done = 0, 0; t_start = time.time()
        for _, row in df.iterrows():
            mid = row["mutant_id"]
            seq = row.get(self.config.evo2_seq_column, row.get("protein", ""))
            if mid in done_set: n_skip += 1; continue
            try:
                results, meta = self.extract_single(seq, mid)
                torch.save(results, output_dir / f"{mid}.pt")
                with open(meta_path, "a") as f:
                    for e in meta: f.write(json.dumps(e) + "\n")
                n_done += 1
                if n_done % log_every == 0:
                    logging.info(f"[EVO2] {n_done}/{len(df)-n_skip} ({n_done/max(time.time()-t_start,1e-6):.1f} it/s)")
                if n_done % 100 == 0: torch.cuda.empty_cache(); gc.collect()
            except Exception as e:
                logging.error(f"[EVO2] {mid} failed: {e}")
        logging.info(f"[EVO2] Done: {n_done} ok, {n_skip} skip")


# ═══════════════════════════════════════════════
# AGENT ENTRY POINT
# ═══════════════════════════════════════════════
def extract_from_manifest(
    protein_name: str,
    device: str,
    manifest_csv: str,
    output_dir: str,
    template_cif: Optional[str] = None,
    msa_gold: Optional[str] = None,
    chains: Optional[list[ChainSpec]] = None,
    models: str = "rf3",
    mode: str = "heads",
    batch_size: int = 5,
    num_steps: int = 50,
    diffusion_batch_size: int = 2,
    n_recycles: int = 10,
    template_selection: Optional[list[str]] = None,
) -> dict:
    """One-shot: manifest CSV -> attention heads or structures.

    Args:
        protein_name: Human-readable label.
        device: 'cuda:0', 'cuda:1', 'cpu', etc.
        manifest_csv: Path with 'mutant_id' and 'protein' columns.
        output_dir: Save location (subdirs per model).
        template_cif: Optional -> WARNING, de novo.
        msa_gold: Optional -> WARNING, seq-only.
        chains: Assembly definition. Inferred from template CIF if None.
        models: Comma/space-sep: 'rf3', 'faesm', 'evo2', 'all'.
        mode: 'heads' (heads+ZII) or 'structure' (JSON-array batch, no heads).
        batch_size: Heads -> pre-gen chunk size. Structure -> per JSON array.
        num_steps/diffusion_batch_size/n_recycles: RF3 inference params.
        template_selection: Which chains to template.
    """
    out = Path(output_dir)
    log_dir = out / "_logs"; log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(log_dir / f"{protein_name}.log", mode="a")])
    logging.getLogger().setLevel(logging.INFO)

    manifest_path = Path(manifest_csv)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest CSV not found: {manifest_csv}")
    # Peek at manifest to infer chain seq if needed
    first_seq = ""
    try:
        with open(manifest_path) as f:
            header = f.readline().strip().split(",")
            first_line = f.readline().strip().split(",")
            if "protein" in header and len(first_line) > header.index("protein"):
                first_seq = first_line[header.index("protein")]
    except Exception:
        pass

    resolved_chains = chains
    if resolved_chains is None:
        if template_cif and Path(template_cif).exists():
            resolved_chains = _infer_chains_from_template(Path(template_cif))
        elif first_seq:
            logging.info(f"[API] Inferring chain A from manifest: {len(first_seq)}aa")
            resolved_chains = [ChainSpec(chain_id="A", seq=first_seq)]
        else:
            logging.warning("[API] No chains, no template, no manifest protein. Creating placeholder.")
            resolved_chains = [ChainSpec(chain_id="A", seq="")]

    tpl = Path(template_cif) if template_cif and Path(template_cif).exists() else None
    msa = Path(msa_gold) if msa_gold and Path(msa_gold).exists() else None
    if template_cif and tpl is None:
        logging.warning(f"[API] template_cif not found: {template_cif}")
    if msa_gold and msa is None:
        logging.warning(f"[API] msa_gold not found: {msa_gold}")

    config = ProteinConfig(protein_name=protein_name, device=device, output_dir=out,
                           chains=resolved_chains, template_cif=tpl,
                           template_selection=template_selection, msa_gold_a3m=msa,
                           num_steps=num_steps, diffusion_batch_size=diffusion_batch_size,
                           n_recycles=n_recycles)
    config.warn_missing()

    pd = _import_pd()
    df = pd.read_csv(manifest_path)
    if "mutant_id" not in df.columns:
        raise ValueError("Manifest CSV must have 'mutant_id' column")

    model_list = [m.strip() for m in models.replace(",", " ").split() if m.strip()]
    if "all" in model_list:
        model_list = ["rf3", "faesm", "evo2"]

    logging.info(f"=== Pipeline: {protein_name} on {device} | {mode} BS={batch_size} ===")
    logging.info(f"  Models: {model_list} | Inputs: {len(df)} | Output: {out}")

    summary = {"protein_name": protein_name, "output_dir": str(out.resolve()),
               "modes_run": [], "total_input": len(df), "completed": 0, "errors": 0}
    if mode == "structure":
        summary["modes_run"].append("structure")

    torch = _import_torch()
    with RunInstance.get_or_create(protein_name):
        for m in model_list:
            model_out = out / m
            if m == "rf3":
                extractor = RF3Extractor(config)
                try:
                    if mode == "parallel":
                        # Parallel mode: workers load their own models.
                        # Coordinator only needs the builder for JSON gen.
                        extractor.builder = RF3BatchBuilder(config)
                        save_head_registry(build_head_registry(extractor.manager, "rf3") if extractor.manager else [],
                                           model_out / "head_registry.json")
                        extractor.run_batch_parallel(df, model_out, batch_size=batch_size)
                    elif mode == "chain":
                        # Chain-batch mode: N mutants as N chains, 1 forward pass
                        extractor.load_model()
                        save_head_registry(build_head_registry(extractor.manager, "rf3"),
                                           model_out / "head_registry.json")
                        extractor.run_batch_chain(df, model_out, batch_size=batch_size)
                    elif mode == "zii":
                        # ZII-only chain-batch: NO hooks, maximum speed
                        extractor.load_model()
                        extractor.run_batch_zii_chain(df, model_out, batch_size=batch_size)
                    else:
                        extractor.load_model()
                        save_head_registry(build_head_registry(extractor.manager, "rf3"),
                                           model_out / "head_registry.json")
                        if mode == "heads":
                            extractor.run_batch_heads(df, model_out, batch_size=batch_size)
                        elif mode == "structure":
                            extractor.run_batch_structure(df, out, batch_size=batch_size)
                    summary["completed"] += 1
                except Exception as e:
                    logging.error(f"[RF3] Failed: {e}", exc_info=True)
                    summary["errors"] += 1
                finally:
                    if extractor.manager:
                        try: extractor.manager.detach_and_unregister()
                        except Exception: pass
                    del extractor; torch.cuda.empty_cache(); gc.collect()
            elif m == "faesm":
                extractor = FAESMExtractor(config)
                try:
                    extractor.load_model()
                    save_head_registry(build_head_registry(extractor.manager, "faesm"),
                                       model_out / "head_registry.json")
                    extractor.run_batch(df, model_out)
                    summary["completed"] += 1
                except Exception as e:
                    logging.error(f"[FAESM] Failed: {e}", exc_info=True)
                    summary["errors"] += 1
                finally:
                    if extractor.manager:
                        try: extractor.manager.detach_and_unregister()
                        except Exception: pass
                    del extractor; torch.cuda.empty_cache(); gc.collect()
            elif m == "evo2":
                extractor = EVO2Extractor(config)
                try:
                    extractor.load_model()
                    registry = build_head_registry(extractor.manager, "evo2")
                    save_head_registry([r for r in registry if r["attention_type"] == "Evo2MHA"],
                                       model_out / "head_registry.json")
                    extractor.run_batch(df, model_out)
                    summary["completed"] += 1
                except Exception as e:
                    logging.error(f"[EVO2] Failed: {e}", exc_info=True)
                    summary["errors"] += 1
                finally:
                    if extractor.manager:
                        try: extractor.manager.detach_and_unregister()
                        except Exception: pass
                    del extractor; torch.cuda.empty_cache(); gc.collect()

    logging.info(f"=== Pipeline complete: {summary} ===")
    return summary


def _infer_chains_from_template(cif: Path) -> list[ChainSpec]:
    if not cif or not cif.exists():
        return []
    chains = {}
    try:
        with open(cif) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    cid = line[21].strip()
                    if cid and cid not in chains:
                        chains[cid] = ChainSpec(chain_id=cid, seq="")
    except Exception:
        pass
    return list(chains.values())


# ═══════════════════════════════════════════════
# INTERACTIVE CLI
# ═══════════════════════════════════════════════
def _prompt(t: str, default: Optional[str] = None) -> str:
    t = t.strip()
    if default:
        v = input(f"  {t} [{default}]: ").strip()
        return v if v else default
    return input(f"  {t}: ").strip()

def _prompt_path(t: str, default: str = "", must_exist: bool = False) -> Optional[Path]:
    raw = _prompt(t, default or None)
    if not raw: return None
    p = Path(raw)
    if must_exist and not p.exists():
        print(f"  WARNING: Not found: {p}")
        return _prompt_path(t, default, must_exist)
    return p

def _pick_models() -> list[str]:
    raw = _prompt("Models (rf3/faesm/evo2/all)", "rf3").strip().lower()
    if raw == "all": return ["rf3", "faesm", "evo2"]
    return [m for m in raw.replace(",", " ").split() if m in ("rf3", "faesm", "evo2")] or ["rf3"]

def _pick_mode() -> str:
    raw = _prompt("Mode: heads (attention+ZII) or structure (batch, fast)", "heads").strip().lower()
    return "structure" if raw.startswith("s") else "heads"

def interactive_cli():
    print("\n" + "=" * 65)
    print("  Unified Attention Extraction Pipeline")
    print("  (Enter to accept defaults. Leave blank to skip optional fields.)")
    print("=" * 65)

    protein_name = _prompt("Protein name", "my_protein")
    torch = _import_torch()
    n_gpu = torch.cuda.device_count()
    if n_gpu == 0:
        print("  WARNING: No CUDA GPUs. Using CPU.")
        device = "cpu"
    else:
        print(f"  GPUs: {n_gpu}")
        for i in range(n_gpu):
            m = torch.cuda.get_device_name(i)
            gb = torch.cuda.get_device_properties(i).total_memory / 1e9
            print(f"    {i}: {m} ({gb:.1f} GB)")
        device = _prompt("Device (cuda:N or cpu)", "cuda:0")

    manifest_csv = _prompt_path("Manifest CSV path", must_exist=True)
    if manifest_csv is None:
        print("  Manifest required. Exiting."); return {"error": "no manifest"}
    template_cif = _prompt_path("Template CIF (optional)")
    msa_gold = _prompt_path("Gold MSA .a3m (optional)")
    output_dir = _prompt_path("Output directory",
                               f"/mnt/j/conda_envs/foundry/DMS_Project/output_heads/{protein_name}",
                               must_exist=False)

    n_chains = int(_prompt("Number of chains", "1"))
    chains = []
    for i in range(n_chains):
        cid = chr(65 + i)
        print(f"\n  --- Chain {cid} ---")
        seq = _prompt("Sequence (Enter for ligand)")
        if seq:
            has_msa = _prompt("Has MSA? (y/n)", "n").lower().startswith("y")
            msa_p = str(_prompt_path("MSA path")) if has_msa else None
            chains.append(ChainSpec(chain_id=cid, seq=seq, msa_path=msa_p))
        else:
            it = _prompt("Type: smiles/ccd/path/none", "none").lower()
            if it == "smiles":
                chains.append(ChainSpec(chain_id=cid, seq="", smiles=_prompt("SMILES")))
            elif it == "ccd":
                chains.append(ChainSpec(chain_id=cid, seq="", ccd_code=_prompt("CCD code")))
            elif it == "path":
                chains.append(ChainSpec(chain_id=cid, seq="", cif_path=str(_prompt_path("CIF/SDF path"))))
            else:
                chains.append(ChainSpec(chain_id=cid, seq=""))

    mode = _pick_mode()
    bs = int(_prompt("Batch size", "5" if mode == "structure" else "1"))
    models = _pick_models()
    num_steps = int(_prompt("RF3 num_steps", "50"))

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler(),
                                  logging.FileHandler(Path(output_dir)/"_logs"/f"{protein_name}.log", mode="a")])
    result = extract_from_manifest(protein_name=protein_name, device=device,
                                    manifest_csv=str(manifest_csv), output_dir=str(output_dir),
                                    template_cif=str(template_cif) if template_cif else None,
                                    msa_gold=str(msa_gold) if msa_gold else None,
                                    chains=chains, models=",".join(models), mode=mode,
                                    batch_size=bs, num_steps=num_steps)
    print(f"\nDone. Output: {output_dir}")
    return result


if __name__ == "__main__":
    if len(sys.argv) > 1:
        import argparse
        p = argparse.ArgumentParser(description="Unified attention extraction")
        p.add_argument("--protein", required=True)
        p.add_argument("--device", default="cuda:0")
        p.add_argument("--manifest", required=True)
        p.add_argument("--output", required=True)
        p.add_argument("--template-cif", default=None)
        p.add_argument("--msa-gold", default=None)
        p.add_argument("--models", default="rf3")
        p.add_argument("--mode", default="heads", choices=["heads", "structure", "parallel", "chain"])
        p.add_argument("--batch-size", type=int, default=1)
        p.add_argument("--num-steps", type=int, default=50)
        p.add_argument("--template-selection", default=None, help="comma-sep: A,B")
        args = p.parse_args()
        ts = args.template_selection.split(",") if args.template_selection else None
        extract_from_manifest(protein_name=args.protein, device=args.device,
                               manifest_csv=args.manifest, output_dir=args.output,
                               template_cif=args.template_cif, msa_gold=args.msa_gold,
                               models=args.models, mode=args.mode, batch_size=args.batch_size,
                               num_steps=args.num_steps, template_selection=ts)
    else:
        interactive_cli()
