"""Shared path configuration for the RF3 MSA/template audit repository.

All lightweight-reproduction paths resolve relative to the repository root
so that analysis scripts run on any machine with the repo checked out.
Heavy inference assets (RF3 code root, checkpoints, project feature roots)
are NOT part of the repo and are configured through environment variables;
see analysis/README.md for the full list.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = DATA_DIR / "results"
INPUTS_DIR = DATA_DIR / "inputs"
STRUCTURES_DIR = DATA_DIR / "structures"
FIGURES_DIR = REPO_ROOT / "figures"

# Lightweight-reproduction inputs/outputs
E2_MSA_DIR = INPUTS_DIR / "e2_targets" / "msas"
E49_TPL_DIR = INPUTS_DIR / "e49_templates"
E428_STRUCTURES_DIR = STRUCTURES_DIR / "e428"
E428_RESULTS = RESULTS_DIR / "e428" / "results.json"

# Heavy inference assets (external, configured via environment)
RF3_ROOT = os.environ.get("RF3_ROOT", "")          # RF3/foundry model code root
RF3_CHECKPOINT = os.environ.get("RF3_CHECKPOINT", "")  # RF3 remapped checkpoint
RF3_CONFIG_DIR = os.environ.get("RF3_CONFIG_DIR", "")  # RF3 hydra configs (default: RF3_ROOT/models/rf3/configs)


def rf3_config_dir():
    if RF3_CONFIG_DIR:
        return RF3_CONFIG_DIR
    if RF3_ROOT:
        return str(Path(RF3_ROOT) / "models" / "rf3" / "configs")
    return ""


def require_rf3():
    """Raise a clear error if RF3 inference assets are not configured."""
    missing = [v for v, val in (("RF3_ROOT", RF3_ROOT), ("RF3_CHECKPOINT", RF3_CHECKPOINT)) if not val]
    if missing:
        raise RuntimeError(
            "RF3 inference assets not configured; set environment variables: "
            + ", ".join(missing)
            + " (optional: RF3_CONFIG_DIR). See analysis/README.md."
        )
    return Path(RF3_ROOT), Path(RF3_CHECKPOINT)


def env_path(var, default):
    """Path from environment variable `var`, falling back to `default`."""
    return Path(os.environ[var]) if var in os.environ else Path(default)
