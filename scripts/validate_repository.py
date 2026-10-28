#!/usr/bin/env python3
"""Repository integrity check for the RF3 MSA/template audit (analysis utility).

Verifies that the canonical results and archived structures shipped in the
repo match the audited experiment dimensions:

  E2          13 entries in e2/e2_rangematched_results.json
  E428        480 rows in e428/results.json, 30 cells x 16 structures each,
              and 480 CIFs under structures/e428
  E49-RBD     60 runs in e49_rbd_supervised/results_e49.json (4 arms)
  E427        6 result files in e427/
  E429        2 protocols x 8 runs (results_A.json, results_B.json)
  Figures     Fig1/Fig2 present in figures/ (pdf + png)

Usage:  python scripts/validate_repository.py [--repo ROOT]
Exits 0 with a PASS summary if every check holds, 1 otherwise.
No GPU, training, or network access required.
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


@check("E2 entries == 13")
def e2_entries(repo):
    rows = json.loads((repo / "data/results/e2/e2_rangematched_results.json").read_text())
    assert len(rows) == 13, f"got {len(rows)}"
    assert len({r["id"] for r in rows}) == 13, "duplicate ids"
    return f"13 ids: {rows[0]['id']}..{rows[-1]['id']}"


@check("E49b GA/GB results present")
def e49b(repo):
    d = json.loads((repo / "data/results/e2/e49b_results.json").read_text())
    assert "GA88" in d and "GB88" in d, f"got {sorted(d)}"
    for k in ("GA88", "GB88"):
        assert "msa_tplOther" in d[k], f"{k}: got {sorted(d[k])}"
    extra = sorted(set(d) - {"GA88", "GB88"})
    return "GA88 + GB88 msa_tplOther values present" + (f" (extra keys: {extra})" if extra else "")


@check("E428 results == 480 rows, 30 cells x 16")
def e428_results(repo):
    rows = json.loads((repo / "data/results/e428/results.json").read_text())
    assert len(rows) == 480, f"got {len(rows)} rows"
    cells = defaultdict(list)
    for r in rows:
        cells[(r["tid"], r["arm"])].append(r)
    assert len(cells) == 30, f"got {len(cells)} cells"
    bad = {k: len(v) for k, v in cells.items() if len(v) != 16}
    assert not bad, f"cells without 16 rows: {bad}"
    return "30 cells (4x5 E2 + 2x5 GA/GB) x 16 structures"


@check("E428 structures == 480 CIFs (30x16)")
def e428_structures(repo):
    base = repo / "data/structures/e428"
    tids = sorted(p.name for p in base.iterdir() if p.is_dir())
    assert len(tids) == 6, f"got {tids}"
    total = 0
    per_target = {}
    for tid in tids:
        n = sum(1 for _ in base.glob(f"{tid}/*/*/*.cif"))
        per_target[tid] = n
        total += n
    assert total == 480, f"got {total}: {per_target}"
    assert all(n == 80 for n in per_target.values()), per_target
    return "6 targets x 80 CIFs = 480"


@check("E49-RBD runs == 60 (4 arms)")
def e49_rbd(repo):
    d = json.loads((repo / "data/results/e49_rbd_supervised/results_e49.json").read_text())
    runs = d["runs"]
    assert len(runs) == 60, f"got {len(runs)}"
    arms = defaultdict(int)
    for r in runs:
        arms[r["cond"]] += 1
    assert len(arms) == 4 and all(n == 15 for n in arms.values()), dict(arms)
    return "4 arms x 15 runs = 60"


@check("E427 result files == 6")
def e427(repo):
    files = sorted((repo / "data/results/e427").glob("*.json"))
    assert len(files) == 6, f"got {[f.name for f in files]}"
    names = [f.stem for f in files]
    return ", ".join(names)


@check("E429 2 protocols x 8 runs")
def e429(repo):
    out = []
    for proc in ("A", "B"):
        d = json.loads((repo / f"data/results/e429/results_{proc}.json").read_text())
        assert len(d["runs"]) == 8, f"{proc}: got {len(d['runs'])}"
        out.append(f"{proc}: 8 runs mean={d['mean']:.4f}")
    return "; ".join(out)


@check("Figures Fig1 + Fig2 (pdf/png)")
def figures(repo):
    figs = []
    for stem in ("Fig1_conformational_selection", "Fig2_template_conflict"):
        for ext in (".pdf", ".png"):
            p = repo / "figures" / (stem + ext)
            assert p.is_file() and p.stat().st_size > 0, f"missing/empty: {p}"
        figs.append(stem)
    return " + ".join(figs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parents[1]),
                    help="repository root (default: %(default)s)")
    args = ap.parse_args()
    repo = Path(args.repo)
    failed = 0
    for name, fn in CHECKS:
        try:
            detail = fn(repo)
            print(f"PASS  {name:<38} {detail}")
        except Exception as e:
            failed += 1
            print(f"FAIL  {name:<38} {e}")
    print(f"\n{'PASS' if failed == 0 else 'FAIL'} — {len(CHECKS) - failed}/{len(CHECKS)} checks ok"
          f"  ({'all canonical dimensions verified' if failed == 0 else 'see failures above'})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
