#!/usr/bin/env python3
"""E429: same-init, same-split, same-cell training noise floor for M1 sec.2.5.

M1's +/-0.02 floor was calibrated on UNPAIRED cross-process/cross-GPU cell
reruns (E14). The gradient-rebalancing claims in sec.2.5 (+0.013/+0.014) are
SPLIT-PAIRED intervention deltas, so the relevant noise is the same-init
same-split rerun spread of ONE cell, which has never been measured.

Cell: RBD x_msa_faesm (CrossAttnOrthoConcatFusion(Z_msa, FAESM)), split seed
100, init 7, canonical 300-epoch protocol (identical to AB-V3/AB-V4's base).
Design: process A and process B each run K identical replicates IN ONE
PROCESS (torch.manual_seed(init) per rep => algorithmically deterministic;
residual variation = GPU nondeterminism). Report per-process sd and the
pooled same-init spread. Usage:
  CUDA_VISIBLE_DEVICES=0 python -B e429_paired_noise_floor.py --proc A --reps 8
Output: /mnt/k/output_heads/e429_paired_noise_floor/results_{proc}.json
"""
import argparse, json, time
from pathlib import Path
import numpy as np

OUT = Path('/mnt/k/output_heads/e429_paired_noise_floor')
OUT.mkdir(exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--proc', required=True)
    ap.add_argument('--reps', type=int, default=8)
    args = ap.parse_args()
    import e2_four_condition as e2
    from fusion_v3 import CrossAttnOrthoConcatFusion
    t0 = time.time()

    def log(m):
        print(f'[{time.time()-t0:7.1f}s] {m}', flush=True)

    Z_msa = e2.std_per_residue(e2.load_z_msa()).to(e2.DEVICE)
    log('z_msa ready')
    F = e2.std_per_residue(e2.load_faesm()).to(e2.DEVICE)
    log('faesm ready')
    tr, va = e2.cross_split(100)
    res_path = OUT / f'results_{args.proc}.json'
    runs = []
    if res_path.exists():
        runs = json.load(open(res_path))['runs']
    for rep in range(len(runs), args.reps):
        t1 = time.time()
        r = e2.train_one(CrossAttnOrthoConcatFusion(), [Z_msa, F], tr, va, 7)
        runs.append({'rep': rep, 'val': float(r), 'sec': time.time() - t1})
        log(f'proc{args.proc} rep{rep}: {r:.4f} ({time.time()-t1:.0f}s)')
        with open(res_path, 'w') as f:
            json.dump({'runs': runs,
                       'mean': float(np.mean([x['val'] for x in runs])),
                       'sd': float(np.std([x['val'] for x in runs], ddof=1)) if len(runs) > 1 else None,
                       'range': float(np.ptp([x['val'] for x in runs]))}, f, indent=1)
    vals = [x['val'] for x in runs]
    log(f'proc{args.proc}: n={len(vals)} mean={np.mean(vals):.4f} '
        f'sd={np.std(vals, ddof=1):.4f} range={np.ptp(vals):.4f}')


if __name__ == '__main__':
    main()
