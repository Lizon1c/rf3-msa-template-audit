#!/usr/bin/env python3
"""E2 range-matched reanalysis: score both conformational states on a SHARED
residue universe (2026-08-02).

Bug fixed (flagged in README §6 / findings §6 / paper ledger L18):
e2_analyze.py scored each state's NLL over its own npy span, so the two sides
of pref = NLL(B) - NLL(A) covered different residue sets —
  XCL1 (E2-04): stateA 72 res (auth 22-93) vs stateB 39 res (22-60)
  RfaH (E2-02): stateA 141 res (auth 2-156, NTD+CTD) vs stateB 66 res (97-162,
  CTD only)
plus a second, unreported bug: RfaH stateA (2OUG) lacks the disordered loop
UniProt 101-114, so its npy rows are COMPRESSED and the contiguous-slice
assumption (p[off:off+M]) misaligns rows >= 99 by 14 residues.
Milder span mismatches also exist on KaiB (98v106), Mad2 (187v205),
LacY (417v391), Abl1 (272v263).

Method: per target, per state, align the npy row sequence (structs/{id}.json
'sequence') to the folded sequence (a3m row0, domain-sliced like
e49_analyze_e2.py) with difflib matching blocks (>=5 residues, kills tag
artifacts like RfaH stateB's GAMG). Shared universe = fold-sequence positions
covered by BOTH states. Metrics (NLL/MAE/entropy) recomputed per cond on the
pair set finite in BOTH experimental maps — identical pairs on both sides of
every preference. No new RF3 forwards: archived full-length disto.pt files
suffice (template-arm distos carry query tokens first, templates appended).

Conds: nomsa/msa for all 13 targets; tplA/tplB/msa_tplA for the 5 E49
targets. Output: /mnt/k/output_heads/e2/e2_rangematched_results.json +
comparison table vs the archived (buggy) e2_results.json / e49_results.json.
"""
import json
import re
import difflib
import numpy as np
import torch
from pathlib import Path

PROJ = Path('/mnt/j/conda_envs/foundry/DMS_Project')
TGT_JSON = PROJ / 'inputs' / 'e2_targets' / 'e2_targets.json'
STRUCT_DIR = PROJ / 'inputs' / 'e2_targets' / 'structs'
MSA_DIR = PROJ / 'inputs' / 'e2_targets' / 'msas'
BASE = Path('/mnt/k/output_heads/e2')
OUT = BASE / 'e2_rangematched_results.json'
E49_TIDS = ['E2-01', 'E2-02', 'E2-03', 'E2-04', 'E2-10']
EDGES = np.linspace(2.0, 22.0, 64)
BINW = EDGES[1] - EDGES[0]
CENTERS = np.concatenate([[EDGES[0] - BINW / 2],
                          (EDGES[:-1] + EDGES[1:]) / 2,
                          [EDGES[-1] + BINW / 2]])


def fold_seq_of(tid, t):
    lines = open(MSA_DIR / f'{tid}.a3m').read().split('\n')
    seq = lines[1].strip()
    dom_lo = 0
    if len(seq) > 600:
        m_ = re.match(r'\s*(\d+)\s*-\s*(\d+)', t.get('construct', {}).get('boundaries', ''))
        if m_:
            dom_lo = int(m_.group(1)) - 1
            seq = seq[dom_lo:int(m_.group(2))]
    return seq, dom_lo


def row_map(state_seq, fold_seq, min_block=5):
    """map[j] = fold-seq index of npy row j (None if unaligned)."""
    sm = difflib.SequenceMatcher(None, state_seq, fold_seq, autojunk=False)
    mp = [None] * len(state_seq)
    for a, b, n in sm.get_matching_blocks():
        if n < min_block:
            continue
        for k in range(n):
            mp[a + k] = b + k
    return mp


def state_load(tid, sk, smeta, fold_seq):
    sm = smeta['states'][sk]
    dmap = np.load(STRUCT_DIR / sm['npy']).astype(np.float64)
    mp = row_map(sm.get('sequence', ''), fold_seq)
    return dmap, mp


def metrics_on(probs, d, valid):
    idx = np.clip(np.digitize(d[valid], EDGES), 0, 64)
    pv = probs[valid]
    pr = pv[np.arange(len(idx)), idx]
    nll = float(-np.log(np.clip(pr, 1e-9, 1)).mean())
    exp_d = (probs * CENTERS).sum(-1)
    mae = float(np.abs(exp_d[valid] - d[valid]).mean())
    ent = float(-(probs * np.log(np.clip(probs, 1e-9, 1))).sum(-1)[valid].mean())
    return {'nll': nll, 'mae': mae, 'entropy': ent, 'n_pairs': int(valid.sum())}


def main():
    targets = {t['id']: t for t in json.load(open(TGT_JSON))['targets']}
    old_e2 = {r['id']: r for r in json.load(open(BASE / 'e2_results.json'))}
    old_e49 = {r['id']: r for r in json.load(open(BASE / 'e49_results.json'))} \
        if (BASE / 'e49_results.json').exists() else {}
    results = []
    for tid in sorted(targets):
        t = targets[tid]
        smp = STRUCT_DIR / f'{tid}.json'
        if not smp.exists():
            continue
        smeta = json.load(open(smp))
        if not all(sk in smeta.get('states', {}) for sk in ('stateA', 'stateB')):
            continue
        fold_seq, _ = fold_seq_of(tid, t)
        dA, mpA = state_load(tid, 'stateA', smeta, fold_seq)
        dB, mpB = state_load(tid, 'stateB', smeta, fold_seq)
        posA = {p for p in mpA if p is not None}
        posB = {p for p in mpB if p is not None}
        shared = sorted(posA & posB)
        if len(shared) < 8:
            print(f'{tid}: shared universe too small ({len(shared)} res) — skipped')
            continue
        rowsA = [j for j, p in enumerate(mpA) if p in set(shared)]
        rowsB = [j for j, p in enumerate(mpB) if p in set(shared)]
        assert [mpA[j] for j in rowsA] == shared and [mpB[j] for j in rowsB] == shared
        dAs = dA[np.ix_(rowsA, rowsA)]
        dBs = dB[np.ix_(rowsB, rowsB)]
        both = np.isfinite(dAs) & np.isfinite(dBs)
        if both.sum() < 50:
            print(f'{tid}: shared pair set too small ({both.sum()}) — skipped')
            continue
        row = {'id': tid, 'name': t['name'], 'crystallized': t.get('crystallized_state'),
               'category': t.get('category'),
               'state_a_label': t.get('state_a', {}).get('label', ''),
               'state_b_label': t.get('state_b', {}).get('label', ''),
               'shared_res': len(shared), 'shared_pairs': int(both.sum()),
               'spanA_res': dA.shape[0], 'spanB_res': dB.shape[0]}
        conds = ['nomsa', 'msa'] + (['tplA', 'tplB', 'msa_tplA'] if tid in E49_TIDS else [])
        for cond in conds:
            dp = BASE / tid / cond / 'disto.pt'
            if not dp.exists():
                row[cond] = None
                continue
            logits = torch.load(dp, map_location='cpu', weights_only=True)
            p = torch.softmax(logits.float(), dim=-1).numpy()
            pS = p[np.ix_(shared, shared, np.arange(p.shape[2]))]
            per = {'A': metrics_on(pS, dAs, both), 'B': metrics_on(pS, dBs, both)}
            row[cond] = per
            row[f'pref_{cond}'] = per['B']['nll'] - per['A']['nll']
        if row.get('pref_nomsa') is not None and row.get('pref_msa') is not None:
            row['shift'] = row['pref_msa'] - row['pref_nomsa']
        if row.get('pref_tplA') is not None and row.get('pref_tplB') is not None:
            row['steer_tpl'] = row['pref_tplB'] - row['pref_tplA']
            row['conflict'] = row['pref_msa_tplA'] - row['pref_tplA']
        results.append(row)

    with open(OUT, 'w') as f:
        json.dump(results, f, indent=1)

    def verdict(r):
        cr = (r['crystallized'] or '').lower()
        al = (r.get('state_a_label') or '').lower()
        bl = (r.get('state_b_label') or '').lower()
        side = 'A' if (cr and (cr in al or al in cr)) or not (cr and (cr in bl or bl in cr)) else 'B'
        if r.get('shift') is None:
            return side, '?'
        toward = (side == 'A' and r['shift'] > 0) or (side == 'B' and r['shift'] < 0)
        return side, 'TOWARD' if toward else 'away'

    print(f"\n{'tgt':<7}{'res':>5}{'pairs':>7} {'old_pref_no':>11}{'new_pref_no':>11}"
          f"{'old_pref_msa':>12}{'new_pref_msa':>12}{'old_shift':>10}{'new_shift':>10}  verdict(old->new)")
    n_flip = 0
    for r in results:
        o = old_e2.get(r['id'], {})
        side, vn = verdict(r)
        vo = '?'
        if o.get('shift') is not None:
            cr = (r['crystallized'] or '').lower()
            bl = (r.get('state_b_label') or '').lower()
            toward_o = (side == 'A' and o['shift'] > 0) or (side == 'B' and o['shift'] < 0)
            vo = 'TOWARD' if toward_o else 'away'
        flip = '  <<< FLIP' if (vo != '?' and vn != '?' and vo != vn) else ''
        if flip:
            n_flip += 1
        fmt = lambda v: f"{v:>10.3f}" if v is not None else f"{'—':>10}"
        print(f"{r['id']:<7}{r['shared_res']:>5}{r['shared_pairs']:>7} "
              f"{fmt(o.get('pref_nomsa')):>11}{fmt(r.get('pref_nomsa')):>11}"
              f"{fmt(o.get('pref_msa')):>12}{fmt(r.get('pref_msa')):>12}"
              f"{fmt(o.get('shift')):>10}{fmt(r.get('shift')):>10}  {vo}->{vn}{flip}")
    print(f'\nVerdict flips: {n_flip}')
    e49 = [r for r in results if r.get('steer_tpl') is not None]
    if e49:
        print(f"\nE49 arms (range-matched): {'tgt':<7}{'steer_tpl':>10}{'steer_msa':>10}{'conflict':>9}  vs old")
        for r in e49:
            o = old_e49.get(r['id'], {})
            fmt = lambda v: f"{v:>9.3f}" if v is not None else f"{'—':>9}"
            print(f"{r['id']:<7}{fmt(r['steer_tpl']):>10}{fmt(r.get('shift')):>10}{fmt(r['conflict']):>9}"
                  f"  old: steer_tpl={fmt(o.get('steer_tpl'))} steer_msa={fmt(o.get('steer_msa'))} conflict={fmt(o.get('conflict'))}")
    print('\nsaved', OUT)


if __name__ == '__main__':
    main()
