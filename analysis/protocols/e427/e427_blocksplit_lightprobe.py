#!/usr/bin/env python3
"""E427: robustness of the MSA/template functional effects (M1/M2 sidepapers) to
(A) contiguous-position-block cross-validation and (B) low-capacity probes.

Motivation: M3 (PONE-D-26-44955) reviewer #1 — random-position holdout lets
neighbouring positions leak, and a 3-layer transformer probe on small N may
overfit. The same two criticisms apply to M1's canonical Delta (sp_msa - sp_no)
and M2's RBD template arm. This experiment re-tests both under:

  A1  RBD full N=3998, SinglePredictorV2, arms no/msa, CONTIGUOUS 5-fold position
      blocks (np.array_split of sorted mutated positions; full CV, ~20% val/fold)
  A2  RBD 500-subset (every 8th mutant, metadata file order), arms
      no/msa/tpl/msatpl (M2 sec.2.3), same contiguous CV
  A3  CASP3 full set, SinglePredictorV2, arms no/msa, same contiguous CV
  B_* light probes on CANONICAL RANDOM splits (seeds 100-104, inits 7/107/207) so
      they pair directly with the archived SP numbers: flattened Linear and a
      per-residue 2-layer MLP (no inter-position mixing), arms no/msa, systems
      rbd / casp3 / ci (ci = N=351, the overfitting worst case)

Everything else follows the canonical protocol (e2_four_condition.py /
e4_*_msa_analysis.py): 300 ep, BS=64, AdamW lr 1e-4 wd 0.05, cosine, MSE,
best-epoch val Spearman, raw labels, overall per-position per-channel std
(axis=0; optimistic convention, shared by all cells within each comparison).

GPU is selected BY NAME (5090D) — enumeration has flipped on this box before.
Run experiments SEQUENTIALLY (5090D serial rule). Resume-safe: completed
(cond, probe, split, init) cells are skipped on restart.

Output: /mnt/k/output_heads/e427_blocksplit_lightprobe/results_{EXP}.json
Feature cache: .../cache/{system}_{arm}_{subset}.npy (row-sliced fp32 [N,L,128]).
"""
import argparse, csv, json, os, re, time
import numpy as np
import torch, torch.nn as nn
from pathlib import Path
from scipy import stats

torch.set_float32_matmul_precision('high')
import sys
sys.path.insert(0, '/mnt/j/conda_envs/foundry/DMS_Project')
from fusion_v3 import SinglePredictorV2, AttnPool

PROJ = Path('/mnt/j/conda_envs/foundry/DMS_Project')
HEADS = Path('/mnt/k/output_heads')
OUT = HEADS / 'e427_blocksplit_lightprobe'
CACHE = OUT / 'cache'
CACHE.mkdir(parents=True, exist_ok=True)

BS, LR, WD = 64, 1e-4, 0.05
t0 = time.time()


def log(m):
    print(f'[{time.time()-t0:7.1f}s] {m}', flush=True)


def pick_device():
    for i in range(torch.cuda.device_count()):
        if '5090' in torch.cuda.get_device_name(i):
            log(f'GPU: cuda:{i} = {torch.cuda.get_device_name(i)}')
            return torch.device(f'cuda:{i}')
    raise RuntimeError('5090D not found; refusing to run on an unnamed device')


DEVICE = pick_device()


class FlatLinear(nn.Module):
    def __init__(self, L, d=128):
        super().__init__()
        self.lin = nn.Linear(L * d, 1)

    def forward(self, x):
        return self.lin(x.reshape(x.shape[0], -1)).squeeze(-1)


class ResidueMLP(nn.Module):
    """2-layer per-residue MLP + attention pooling; no inter-position mixing."""

    def __init__(self, d=128, hidden=256, out=64):
        super().__init__()
        self.mlp = nn.Sequential(nn.Linear(d, hidden), nn.GELU(), nn.Dropout(0.1),
                                 nn.Linear(hidden, out), nn.GELU())
        self.pool = AttnPool(out)
        self.head = nn.Linear(out, 1)

    def forward(self, x):
        return self.head(self.pool(self.mlp(x))).squeeze(-1)


# ── metadata ──

def get_meta(system):
    if system == 'rbd':
        meta = [json.loads(l) for l in open(HEADS / 'rbd' / 'metadata.jsonl') if l.strip()]
        y = np.array([m['bind_avg'] for m in meta], np.float32)
        pos = np.array([m['site_rbd'] - 1 for m in meta])
        return meta, y, pos, 201
    if system == 'casp3':
        meta = [json.loads(l) for l in open(HEADS / 'casp3' / 'manifest.jsonl') if l.strip()]
        y = np.array([m['label'] for m in meta], np.float32)
        pos = np.array([m['position'] - 1 for m in meta])
        return meta, y, pos, 488
    if system == 'ci':
        rows = list(csv.DictReader(open(PROJ / 'data' / 'DMS_ProteinGym_substitutions' /
                                        'RPC1_LAMBD_Li_2019_high-expression.csv')))
        meta = [{'mutant_id': f'ci_{i:03d}'} for i in range(len(rows))]
        y = np.array([float(r['DMS_score']) for r in rows], np.float32)
        pos = np.array([int(re.match(r'[A-Z](\d+)[A-Z]', r['mutant']).group(1)) - 1
                        for r in rows])
        return meta, y, pos, 237
    raise ValueError(system)


def feature_path(system, arm, i, m):
    b = HEADS / system
    if system == 'rbd':
        if arm == 'no':
            return b / 'zii' / f"{m['mutant_id']}_zii.pt", None
        return b / {'msa': 'zii_e1_full15', 'tpl': 'zii_e49_tpl',
                    'msatpl': 'zii_e49_msatpl'}[arm] / f'mutant_{i:04d}_zii.pt', None
    if system == 'casp3':
        d = b / ('zii' if arm == 'no' else 'zii_msa_sp69')
        return d / f"{m['mutant_id']}_zii.pt", m['position'] - 1
    if system == 'ci':
        d = b / ('zii_clean' if arm == 'no' else 'zii_msa')
        return d / f'ci_{i:03d}_zii.pt', None
    raise ValueError((system, arm))


def load_features(system, arm, meta, subset, L):
    tag = 'full' if subset is None else f'sub{len(subset)}'
    cache = CACHE / f'{system}_{arm}_{tag}.npy'
    if cache.exists():
        log(f'  cache hit {cache.name}')
        return np.load(cache)
    idx = range(len(meta)) if subset is None else subset
    X = np.empty((len(idx), L, 128), np.float32)
    for j, i in enumerate(idx):
        fp, row = feature_path(system, arm, i, meta[i])
        z = torch.load(fp, map_location='cpu', weights_only=True)
        if isinstance(z, dict):
            z = list(z.values())[0]
        X[j] = (z[0] if row is None else z[row]).float().numpy()
        if j % 500 == 499:
            log(f'  {system}/{arm} {j + 1}/{len(idx)}')
    np.save(cache, X)
    return X


def std_per_channel(X):
    """Canonical overall per-position, per-channel z-score (axis=0 vector check)."""
    X = torch.from_numpy(X)
    for i in range(X.shape[1]):
        mu = X[:, i].mean(0)
        sd = X[:, i].std(0).clamp_min(1e-8)
        X[:, i] = (X[:, i] - mu) / sd
    return X


# ── splits ──

def random_split(pos, n, seed):
    """Canonical cross-position 70/30 (identical to e2/e4 cross_split)."""
    rng = np.random.RandomState(seed)
    ap = sorted(np.unique(pos))
    rng.shuffle(ap)
    vp = set(ap[:max(1, int(len(ap) * 0.3))])
    return (np.array([i for i in range(n) if pos[i] not in vp]),
            np.array([i for i in range(n) if pos[i] in vp]))


def contiguous_folds(pos, n, k=5):
    """Full K-fold CV over contiguous blocks of sorted mutated positions."""
    ap = np.array(sorted(np.unique(pos)))
    blocks = np.array_split(ap, k)
    folds = []
    for b in blocks:
        vp = set(b.tolist())
        folds.append((np.array([i for i in range(n) if pos[i] not in vp]),
                      np.array([i for i in range(n) if pos[i] in vp])))
    return folds


# ── training ──

def train_one(model, X, y, tr, va, init, epochs):
    torch.manual_seed(init)
    model = model.to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    tr_t = torch.from_numpy(tr).to(DEVICE)
    va_t = torch.from_numpy(va).to(DEVICE)
    yt = torch.from_numpy(y[tr]).to(DEVICE)
    Xt, Xv = X[tr_t], X[va_t]
    yv = y[va]
    best, best_ep = -1.0, -1
    tstart = time.time()
    for ep in range(epochs):
        model.train()
        idx = torch.randperm(len(tr), device=DEVICE)
        for bi in range(0, len(tr), BS):
            b = idx[bi:bi + BS]
            loss = nn.MSELoss()(model(Xt[b]), yt[b])
            opt.zero_grad(); loss.backward(); opt.step()
        sched.step()
        model.eval()
        with torch.no_grad():
            pv = torch.cat([model(Xv[i:i + BS]) for i in range(0, len(va), BS)]).cpu().numpy()
        r, _ = stats.spearmanr(pv, yv)
        r = -1.0 if np.isnan(r) else float(r)
        if r > best:
            best, best_ep = r, ep
    del model, Xt, Xv
    torch.cuda.empty_cache()
    return best, best_ep, time.time() - tstart


def summarize(runs, exp):
    summary = {}
    cells = sorted({(r['probe'], r['cond']) for r in runs})
    for probe, cond in cells:
        vals = [r['val'] for r in runs if r['probe'] == probe and r['cond'] == cond]
        summary[f'{probe}:{cond}'] = {'mean': float(np.mean(vals)), 'std': float(np.std(vals)),
                                      'n': len(vals)}
    paired = {}
    for probe in sorted({r['probe'] for r in runs}):
        pr = [r for r in runs if r['probe'] == probe]
        for a, b in ([('msa', 'no')] if not exp.startswith('A2') else
                     [('msa', 'no'), ('tpl', 'no'), ('msatpl', 'no'), ('msatpl', 'tpl')]):
            ds = []
            for sd in sorted({r['split'] for r in pr}):
                va = [r['val'] for r in pr if r['cond'] == a and r['split'] == sd]
                vb = [r['val'] for r in pr if r['cond'] == b and r['split'] == sd]
                if va and vb:
                    ds.append(float(np.mean(va) - np.mean(vb)))
            if ds:
                paired[f'{probe}:{a}-{b}'] = {'per_split': ds, 'mean': float(np.mean(ds))}
    return summary, paired


EXPS = {
    'A1': dict(system='rbd', arms=['no', 'msa'], probes=['sp'], subset=None, split='block'),
    'A2': dict(system='rbd', arms=['no', 'msa', 'tpl', 'msatpl'], probes=['sp'],
               subset='sub500', split='block'),
    'A3': dict(system='casp3', arms=['no', 'msa'], probes=['sp'], subset=None, split='block'),
    'B_rbd': dict(system='rbd', arms=['no', 'msa'], probes=['linear', 'mlp'],
                  subset=None, split='random'),
    'B_casp3': dict(system='casp3', arms=['no', 'msa'], probes=['linear', 'mlp'],
                    subset=None, split='random'),
    'B_ci': dict(system='ci', arms=['no', 'msa'], probes=['linear', 'mlp'],
                 subset=None, split='random'),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp', required=True, choices=sorted(EXPS))
    ap.add_argument('--epochs', type=int, default=300)
    ap.add_argument('--inits', default='7,107,207')
    ap.add_argument('--max-runs', type=int, default=0, help='smoke-test cap; 0 = all')
    ap.add_argument('--probes', default='', help='override probe list, e.g. linear,mlp')
    ap.add_argument('--tag', default='', help='suffix for results file (smoke tests)')
    args = ap.parse_args()
    cfg = EXPS[args.exp]
    if args.probes:
        cfg = dict(cfg, probes=args.probes.split(','))
    system = cfg['system']
    meta, y, pos, L = get_meta(system)
    subset = None
    if cfg['subset'] == 'sub500':
        subset = list(range(0, len(meta), 8))
        meta = [meta[i] for i in subset]
        y, pos = y[subset], pos[subset]
        subset = list(range(len(meta)))  # already sliced; loader uses positional idx
    # loader needs ORIGINAL indices for file naming; rebuild explicit index list
    if cfg['subset'] == 'sub500':
        orig_idx = list(range(0, len(get_meta(system)[0]), 8))
    else:
        orig_idx = None
    n = len(meta)
    log(f'{args.exp}: {system} N={n} L={L} arms={cfg["arms"]} probes={cfg["probes"]} '
        f'split={cfg["split"]} epochs={args.epochs}')

    feats = {}
    for arm in cfg['arms']:
        X = load_features(system, arm, meta if orig_idx is None else get_meta(system)[0],
                          orig_idx, L)
        feats[arm] = std_per_channel(X).to(DEVICE)
        log(f'{arm} ready {tuple(feats[arm].shape)}')

    if cfg['split'] == 'block':
        splits = contiguous_folds(pos, n)
        split_ids = list(range(len(splits)))
    else:
        split_ids = [100, 101, 102, 103, 104]
        splits = [random_split(pos, n, sd) for sd in split_ids]

    res_path = OUT / f'results_{args.exp}{"_" + args.tag if args.tag else ""}.json'
    runs = []
    if res_path.exists():
        runs = json.load(open(res_path))['runs']
    done = {(r['cond'], r['probe'], r['split'], r['init']) for r in runs}
    inits = [int(x) for x in args.inits.split(',')]
    nrun = 0
    for probe in cfg['probes']:
        for arm in cfg['arms']:
            for sd, (tr, va) in zip(split_ids, splits):
                for init in inits:
                    if (arm, probe, sd, init) in done:
                        continue
                    if probe == 'sp':
                        model = SinglePredictorV2(128)
                    elif probe == 'linear':
                        model = FlatLinear(L)
                    else:
                        model = ResidueMLP()
                    r, ep, sec = train_one(model, feats[arm], y, tr, va, init, args.epochs)
                    runs.append({'cond': arm, 'probe': probe, 'split': sd, 'init': init,
                                 'val': r, 'best_ep': ep, 'sec': sec})
                    nrun += 1
                    log(f'{args.exp} {probe}/{arm} split{sd} init{init}: '
                        f'{r:.4f} (best_ep {ep}, {sec:.0f}s)')
                    summary, paired = summarize(runs, args.exp)
                    with open(res_path, 'w') as f:
                        json.dump({'config': {k: v for k, v in cfg.items()},
                                   'epochs': args.epochs, 'runs': runs,
                                   'summary': summary, 'paired': paired}, f, indent=1)
                    if args.max_runs and nrun >= args.max_runs:
                        log(f'smoke cap reached; saved {res_path}')
                        return
    summary, paired = summarize(runs, args.exp)
    for k, v in summary.items():
        log(f'{k}: mu={v["mean"]:.4f} sd={v["std"]:.4f} n={v["n"]}')
    for k, v in paired.items():
        log(f'paired {k}: {["%+.3f" % d for d in v["per_split"]]} mu={v["mean"]:+.4f}')
    log(f'saved {res_path}')


if __name__ == '__main__':
    main()
