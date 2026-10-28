#!/usr/bin/env python3
"""E49-E2 analyzer: template-state steering vs MSA-state steering in RF3.

Extends e2_analyze.py (conditions nomsa/msa) to the E49 template arms
(tplA/tplB/msa_tplA). Per target x condition x state(A/B): NLL of the
experimental Cα distance map under the predicted distogram; preference =
NLL(B) - NLL(A) (positive => state A matches better).

Key contrasts:
  steer_tpl  = pref_tplB - pref_tplA   (template state steering: should be > 0
               if the model follows the supplied template — tplA pulls toward A,
               tplB pulls toward B)
  steer_msa  = pref_msa - pref_nomsa   (archived E2 shift, recomputed here)
  conflict   = pref_msa_tplA - pref_tplA (does adding MSA on top of tplA move
               the preference? Kovalevskiy: deep/strong MSA dominates templates)
Seq source: a3m row0 (no /tmp/e2/seqs dependency).
Output: /mnt/k/output_heads/e2/e49_results.json + stdout table.
"""
import json
import re
import numpy as np
import torch
from pathlib import Path

PROJ = Path('/mnt/j/conda_envs/foundry/DMS_Project')
TGT_JSON = PROJ / 'inputs' / 'e2_targets' / 'e2_targets.json'
STRUCT_DIR = PROJ / 'inputs' / 'e2_targets' / 'structs'
MSA_DIR = PROJ / 'inputs' / 'e2_targets' / 'msas'
BASE = Path('/mnt/k/output_heads/e2')
CONDS = ['nomsa', 'msa', 'tplA', 'tplB', 'msa_tplA']
EDGES = np.linspace(2.0, 22.0, 64)
BINW = EDGES[1] - EDGES[0]
CENTERS = np.concatenate([[EDGES[0] - BINW / 2],
                          (EDGES[:-1] + EDGES[1:]) / 2,
                          [EDGES[-1] + BINW / 2]])


def dist_metrics(disto_logits, dmap, off=0):
    p = torch.softmax(disto_logits.float(), dim=-1).numpy()
    M = dmap.shape[0]
    if off + M > p.shape[0]:
        M = p.shape[0] - off
        dmap = dmap[:M, :M]
    p = p[off:off + M, off:off + M, :]
    d = dmap
    valid = np.isfinite(d)
    if valid.sum() < 50:
        return None
    idx = np.clip(np.digitize(d[valid], EDGES), 0, 64)
    pv = p[valid]
    probs = pv[np.arange(len(idx)), idx]
    nll = float(-np.log(np.clip(probs, 1e-9, 1)).mean())
    exp_d = (p * CENTERS).sum(-1)
    mae = float(np.abs(exp_d[valid] - d[valid]).mean())
    ent = float(-(p * np.log(np.clip(p, 1e-9, 1))).sum(-1)[valid].mean())
    return {'nll': nll, 'mae': mae, 'entropy': ent, 'n_pairs': int(valid.sum())}


def fold_seq_of(tid, t):
    a3m = MSA_DIR / f'{tid}.a3m'
    lines = open(a3m).read().split('\n')
    seq = lines[1].strip()
    dom_lo = 0
    if len(seq) > 600:
        m_ = re.match(r'\s*(\d+)\s*-\s*(\d+)', t.get('construct', {}).get('boundaries', ''))
        if m_:
            dom_lo = int(m_.group(1)) - 1
            seq = seq[dom_lo:int(m_.group(2))]
    return seq, dom_lo


def main():
    targets = {t['id']: t for t in json.load(open(TGT_JSON))['targets']}
    results = []
    for tid in ['E2-01', 'E2-02', 'E2-03', 'E2-04', 'E2-10']:
        t = targets[tid]
        row = {'id': tid, 'name': t['name'], 'crystallized': t.get('crystallized_state'),
               'category': t.get('category'),
               'state_a_label': t.get('state_a', {}).get('label', ''),
               'state_b_label': t.get('state_b', {}).get('label', '')}
        fold_seq, dom_lo = fold_seq_of(tid, t)
        smeta_all = None
        smp = STRUCT_DIR / f'{tid}.json'
        if smp.exists():
            smeta_all = json.load(open(smp))
        for cond in CONDS:
            dp = BASE / tid / cond / 'disto.pt'
            if not dp.exists():
                row[cond] = None
                continue
            logits = torch.load(dp, map_location='cpu', weights_only=True)
            per_state = {}
            for sname, sk in [('A', 'stateA'), ('B', 'stateB')]:
                if smeta_all and sk in smeta_all.get('states', {}):
                    sm = smeta_all['states'][sk]
                    dmap = np.load(STRUCT_DIR / sm['npy'])
                    off = fold_seq.find(sm.get('sequence', '')[:40])
                    if off < 0:
                        rr = sm.get('residue_range_used')
                        off = max(0, (rr[0] - 1) - dom_lo) if rr else 0
                    per_state[sname] = dist_metrics(logits, dmap, off)
                else:
                    per_state[sname] = None
            row[cond] = per_state
        for cond in CONDS:
            if row.get(cond) and row[cond].get('A') and row[cond].get('B'):
                row[f'pref_{cond}'] = row[cond]['B']['nll'] - row[cond]['A']['nll']
        if row.get('pref_tplA') is not None and row.get('pref_tplB') is not None:
            row['steer_tpl'] = row['pref_tplB'] - row['pref_tplA']
        if row.get('pref_nomsa') is not None and row.get('pref_msa') is not None:
            row['steer_msa'] = row['pref_msa'] - row['pref_nomsa']
        if row.get('pref_tplA') is not None and row.get('pref_msa_tplA') is not None:
            row['conflict'] = row['pref_msa_tplA'] - row['pref_tplA']
        results.append(row)

    with open(BASE / 'e49_results.json', 'w') as f:
        json.dump(results, f, indent=1)

    hdr = f"{'target':<8}{'pref_no':>8}{'pref_msa':>9}{'pref_tplA':>10}{'pref_tplB':>10}{'pref_msa_tplA':>13}{'steer_tpl':>10}{'steer_msa':>10}{'conflict':>9}"
    print('\n' + hdr)
    for r in results:
        def g(k):
            v = r.get(k)
            return f"{v:>8.3f}" if v is not None else f"{'—':>8}"
        print(f"{r['id']:<8}{g('pref_nomsa'):>8}{g('pref_msa'):>9}{g('pref_tplA'):>10}{g('pref_tplB'):>10}{g('pref_msa_tplA'):>13}{g('steer_tpl'):>10}{g('steer_msa'):>10}{g('conflict'):>9}")
    print('\nInterpretation guide: pref>0 favors state A. steer_tpl>0 = model follows')
    print('the supplied template; conflict vs steer_msa sign = who dominates (MSA or tpl).')
    print('saved', BASE / 'e49_results.json')


if __name__ == '__main__':
    main()
