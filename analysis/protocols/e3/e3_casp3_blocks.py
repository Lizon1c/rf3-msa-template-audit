#!/usr/bin/env python3
"""E3: CASP3 docking-component test — block decomposition of MSA effect on Z_II.

The "free docking" hypothesis: for the CASP3 homodimer (2x244, 488 tokens), if MSA
acts partly as inter-chain docking pinning, its effect on the pair representation
should be CONCENTRATED in the cross-chain (A-B) block relative to the within-chain
(A-A, B-B) blocks.

Per mutant (stride 4 subset of the common set of casp3/zii/ [noMSA] and
casp3/zii_msa_sp69/ [MSA], full [488,488,128] pair):
  blocks AA = [0:244, 0:244], BB = [244:488, 244:488], AB = [0:244, 244:488]
  per block: cos(flat(Z_msa), flat(Z_no))  (direction preservation)
             ||Z_msa|| / ||Z_no||          (amplitude change)
Aggregate: median/mean per block + paired Wilcoxon (AB vs AA).
Output: /mnt/k/output_heads/casp3/e3_blocks/results.json + stdout table.
"""
import json, os, time
import numpy as np
import torch
from pathlib import Path
from scipy import stats

torch.set_num_threads(16)
t0 = time.time()
BASE = Path('/mnt/k/output_heads/casp3')
OUT = BASE / 'e3_blocks' / os.environ.get('E3_MSA_DIR', 'zii_msa_sp69')
OUT.mkdir(parents=True, exist_ok=True)
L = 244
STRIDE = 4


def log(m):
    print(f'[{time.time()-t0:7.1f}s] {m}', flush=True)


meta = [json.loads(l) for l in open(BASE / 'manifest.jsonl') if l.strip()]
mids = [m['mutant_id'] for m in meta]
sub = [i for i in range(0, len(mids), STRIDE)]
log(f'{len(mids)} mutants total, subset {len(sub)}')

cos = torch.nn.functional.cosine_similarity
rows = []
for k, i in enumerate(sub):
    mid = mids[i]
    p_no = BASE / 'zii' / f'{mid}_zii.pt'
    p_ms = BASE / os.environ.get('E3_MSA_DIR', 'zii_msa_sp69') / f'{mid}_zii.pt'
    if not p_no.exists() or not p_ms.exists():
        continue
    z_no = torch.load(p_no, map_location='cpu', weights_only=True).float()
    z_ms = torch.load(p_ms, map_location='cpu', weights_only=True).float()
    blocks = {
        'AA': (slice(0, L), slice(0, L)),
        'BB': (slice(L, 2 * L), slice(L, 2 * L)),
        'AB': (slice(0, L), slice(L, 2 * L)),
    }
    rec = {'mutant_id': mid}
    for name, (r, c) in blocks.items():
        a = z_no[r, c, :].reshape(-1)
        b = z_ms[r, c, :].reshape(-1)
        rec[f'cos_{name}'] = float(cos(a, b, dim=0))
        rec[f'normratio_{name}'] = float(b.norm() / a.norm().clamp_min(1e-12))
    rows.append(rec)
    if (k + 1) % 50 == 0:
        log(f'{k+1}/{len(sub)}')

log(f'computed {len(rows)} mutants')
res = {'n': len(rows), 'per_mutant': rows}
summary = {}
for name in ['AA', 'BB', 'AB']:
    c = np.array([r[f'cos_{name}'] for r in rows])
    n = np.array([r[f'normratio_{name}'] for r in rows])
    summary[f'cos_{name}'] = {'median': float(np.median(c)), 'mean': float(c.mean()),
                              'q10': float(np.percentile(c, 10)), 'q90': float(np.percentile(c, 90))}
    summary[f'normratio_{name}'] = {'median': float(np.median(n)), 'mean': float(n.mean())}
cos_aa = np.array([r['cos_AA'] for r in rows])
cos_ab = np.array([r['cos_AB'] for r in rows])
w, p = stats.wilcoxon(cos_ab - cos_aa)
summary['wilcoxon_cosAB_minus_cosAA'] = {'W': float(w), 'p': p,
                                         'median_diff': float(np.median(cos_ab - cos_aa))}
res['summary'] = summary
with open(OUT / 'results.json', 'w') as f:
    json.dump(res, f)

print('\n════════ E3 BLOCK TABLE ════════')
for name in ['AA', 'BB', 'AB']:
    s = summary[f'cos_{name}']
    n_ = summary[f'normratio_{name}']
    print(f"{name}: cos med={s['median']:.4f} (q10 {s['q10']:.3f}, q90 {s['q90']:.3f}) | "
          f"normratio med={n_['median']:.4f}")
w_ = summary['wilcoxon_cosAB_minus_cosAA']
print(f"AB - AA cos: median {w_['median_diff']:+.4f}, Wilcoxon p={w_['p']:.2e}")
log('saved')
