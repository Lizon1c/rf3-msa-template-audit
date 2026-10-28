from pathlib import Path
from collections import defaultdict
import argparse
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
OUT = REPO_ROOT / "figures"


def resolve_inputs(data_dir):
    data_dir = Path(data_dir)
    return {
        "e2": data_dir / "results" / "e2" / "e2_rangematched_results.json",
        "e49b": data_dir / "results" / "e2" / "e49b_results.json",
        "e428": data_dir / "results" / "e428" / "results.json",
        "e49_rbd": data_dir / "results" / "e49_rbd_supervised" / "results_e49.json",
        "e427": data_dir / "results" / "e427" / "results_A2.json",
    }


INPUT_PATHS = resolve_inputs(DATA_DIR)
E2_PATH = INPUT_PATHS["e2"]
E49B_PATH = INPUT_PATHS["e49b"]
E428_PATH = INPUT_PATHS["e428"]
E49_RBD_PATH = INPUT_PATHS["e49_rbd"]
E427_PATH = INPUT_PATHS["e427"]

BLUE = "#0077BB"
CYAN = "#33BBEE"
TEAL = "#009988"
ORANGE = "#EE7733"
RED = "#CC3311"
GREY = "#BBBBBB"
DARK_GREY = "#666666"
BLACK = "#000000"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


def load_inputs(paths=None):
    paths = paths or INPUT_PATHS
    e2 = json.loads(paths["e2"].read_text())
    e49b = json.loads(paths["e49b"].read_text())
    e428_rows = json.loads(paths["e428"].read_text())
    e49_rbd = json.loads(paths["e49_rbd"].read_text())
    e427 = json.loads(paths["e427"].read_text())
    assert len(e2) == 13
    assert len(e428_rows) == 480
    assert len(e49_rbd["runs"]) == 60
    e428 = defaultdict(list)
    for row in e428_rows:
        e428[(row["tid"], row["arm"])].append(row)
    assert len(e428) == 30
    assert all(len(rows) == 16 for rows in e428.values())
    return e2, e49b, e428, e49_rbd, e427


def coordinate_cell(e428, tid, arm):
    rows = e428[(tid, arm)]
    return {
        "a": float(np.mean([r["rmsdA"] for r in rows])),
        "b": float(np.mean([r["rmsdB"] for r in rows])),
        "delta": float(np.mean([r["rmsdB"] - r["rmsdA"] for r in rows])),
        "pct_a": 100 * sum(r["closer_to"] == "A" for r in rows) / len(rows),
    }


def panel_label(ax, label, x=-0.16, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def signed_bar_colors(values):
    return [BLUE if value >= 0 else ORANGE for value in values]


def save_figure(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", metadata={"Title": stem, "Author": "Cheng Ding"})
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    tif_path = OUT / f"{stem}.tif"
    fig.savefig(tif_path, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tif_path) as image:
        rgb = Image.new("RGB", image.size, "white")
        rgb.paste(image, mask=image.getchannel("A"))
        rgb.save(tif_path, dpi=(600, 600), compression="tiff_lzw")
    plt.close(fig)


def make_fig1(e2, e428):
    names = {
        "E2-01": "KaiB",
        "E2-02": "RfaH",
        "E2-03": "Mad2 (cryst.=B)",
        "E2-04": "XCL1",
        "E2-05": "selecase",
        "E2-06": "GA88/GB88",
        "E2-07": "β2AR",
        "E2-08": "Abl1",
        "E2-09": "EGFR",
        "E2-10": "LacY",
        "E2-11": "ubiquitin",
        "E2-12": "GB1",
        "E2-13": "thioredoxin",
    }
    oriented = []
    for row in e2:
        sign = -1 if row["id"] == "E2-03" else 1
        oriented.append(sign * row["shift"])
    categories = []
    for row in e2:
        if row["id"] in {"E2-01", "E2-02"}:
            categories.append("Away from crystal")
        elif row["id"] in {"E2-03", "E2-04"}:
            categories.append("Toward crystal")
        elif row["id"] == "E2-06":
            categories.append("Own-fold sharpening")
        else:
            categories.append("Within noise")
    colors = {
        "Away from crystal": BLUE,
        "Toward crystal": ORANGE,
        "Own-fold sharpening": TEAL,
        "Within noise": GREY,
    }

    fig = plt.figure(figsize=(7.2, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.48, 1.05], hspace=0.50)
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(e2))
    ax.axvspan(-0.27, 0.27, color="#F0F0F0", zorder=0)
    ax.axvspan(-0.14, 0.14, color="#E2E2E2", zorder=1)
    ax.barh(y, oriented, color=[colors[c] for c in categories], height=0.66,
            edgecolor="none", zorder=2)
    ax.axvline(0, color=BLACK, linewidth=0.7, zorder=3)
    ax.set_yticks(y, [names[row["id"]] for row in e2])
    ax.invert_yaxis()
    ax.set_xlim(-3.65, 3.65)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel("MSA-induced distogram shift (NLL units)\noriented so positive points toward the crystallized state")
    ax.set_title("Two-state distogram panel (13 entries)", loc="left", pad=32)
    for yi, value in zip(y, oriented):
        if value >= 0:
            ax.text(value + 0.08, yi, f"{value:+.2f}", va="center", ha="left", fontsize=7)
        else:
            ax.text(value - 0.08, yi, f"{value:+.2f}", va="center", ha="right", fontsize=7)
    ax.legend(handles=[Patch(facecolor=colors[key], label=key) for key in colors],
              loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=4,
              frameon=False, handlelength=1.2, columnspacing=1.2)
    panel_label(ax, "A", x=-0.16, y=1.06)

    ax = fig.add_subplot(gs[1])
    systems = [
        ("KaiB", "E2-01"),
        ("RfaH", "E2-02"),
        ("Mad2 (cryst.=B)", "E2-03"),
        ("XCL1", "E2-04"),
        ("GA88", "GA88"),
        ("GB88", "GB88"),
    ]
    y = np.arange(len(systems))
    for yi, (label, tid) in zip(y, systems):
        no = coordinate_cell(e428, tid, "nomsa")
        msa = coordinate_cell(e428, tid, "msa")
        ax.annotate("", xy=(msa["delta"], yi), xytext=(no["delta"], yi),
                    arrowprops={"arrowstyle": "->", "color": "#999999", "lw": 1.1})
        ax.scatter(no["delta"], yi, s=38, marker="o", color=DARK_GREY,
                   edgecolor=BLACK, linewidth=0.4, zorder=3)
        ax.scatter(msa["delta"], yi, s=42, marker="s", color=BLUE,
                   edgecolor=BLACK, linewidth=0.4, zorder=4)
        ax.text(no["delta"], yi - 0.19, f"{no['delta']:+.2f}", ha="center",
                va="bottom", fontsize=6.7, color=DARK_GREY)
        ax.text(msa["delta"], yi + 0.20, f"{msa['delta']:+.2f}; {msa['pct_a']:.0f}% A",
                ha="center", va="top", fontsize=6.7, color=BLUE)
    ax.axvline(0, color=BLACK, linewidth=0.7)
    ax.set_yticks(y, [label for label, _ in systems])
    ax.set_ylim(len(systems) - 0.35, -0.55)
    ax.set_xlim(-14.8, 14.8)
    ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    ax.set_xlabel("Coordinate preference, ΔRMSD = RMSD$_B$ − RMSD$_A$ (Å)\npositive = closer to state A")
    ax.set_title("Explicit coordinate sampling: no MSA → MSA", loc="left", pad=22)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DARK_GREY,
               markeredgecolor=BLACK, markersize=5.5, label="no MSA"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=BLUE,
               markeredgecolor=BLACK, markersize=5.5, label="MSA"),
    ], loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False)
    panel_label(ax, "B", x=-0.16, y=1.06)
    save_figure(fig, "Fig1_conformational_selection")


def make_fig2(e2, e49b, e428, e49_rbd, e427):
    e2_by_id = {row["id"]: row for row in e2}
    comparisons = [
        ("Mad2\ntplA", e2_by_id["E2-03"]["pref_tplA"], coordinate_cell(e428, "E2-03", "tplA")["delta"]),
        ("XCL1\ntplA", e2_by_id["E2-04"]["pref_tplA"], coordinate_cell(e428, "E2-04", "tplA")["delta"]),
        ("RfaH\ntplB", e2_by_id["E2-02"]["pref_tplB"], coordinate_cell(e428, "E2-02", "tplB")["delta"]),
        ("GA88\nMSA+tplOther", e49b["GA88"]["msa_tplOther"], coordinate_cell(e428, "GA88", "msa_tplOther")["delta"]),
        ("GB88\nMSA+tplOther", e49b["GB88"]["msa_tplOther"], coordinate_cell(e428, "GB88", "msa_tplOther")["delta"]),
    ]
    labels = [row[0] for row in comparisons]
    distogram = [row[1] for row in comparisons]
    coordinates = [row[2] for row in comparisons]

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.62, wspace=0.38)
    y = np.arange(len(labels))

    ax = fig.add_subplot(gs[0, 0])
    ax.barh(y, distogram, color=signed_bar_colors(distogram), height=0.62)
    ax.axvline(0, color=BLACK, linewidth=0.7)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(-5.7, 5.7)
    ax.set_xticks([-5, -2.5, 0, 2.5, 5])
    ax.set_xlabel("Distogram preference (NLL units)")
    ax.set_title("Internal distogram", loc="left", pad=8)
    for yi, value in zip(y, distogram):
        ax.text(value + (0.13 if value >= 0 else -0.13), yi, f"{value:+.2f}",
                ha="left" if value >= 0 else "right", va="center", fontsize=7)
    panel_label(ax, "A", x=-0.22, y=1.08)

    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, coordinates, color=signed_bar_colors(coordinates), height=0.62)
    ax.axvline(0, color=BLACK, linewidth=0.7)
    ax.set_yticks(y, [])
    ax.invert_yaxis()
    ax.set_xlim(-15.8, 15.8)
    ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    ax.set_xlabel("Coordinate preference, ΔRMSD (Å)")
    ax.set_title("Sampled 3D coordinates", loc="left", pad=8)
    for yi, value in zip(y, coordinates):
        ax.text(value + (0.38 if value >= 0 else -0.38), yi, f"{value:+.2f}",
                ha="left" if value >= 0 else "right", va="center", fontsize=7)
    panel_label(ax, "B", x=-0.12, y=1.08)

    ax = fig.add_subplot(gs[1, 0])
    arm_info = [
        ("no", "no MSA", DARK_GREY),
        ("zii_e1_full15", "MSA", BLUE),
        ("zii_e49_tpl", "template", ORANGE),
        ("zii_e49_msatpl", "MSA+\ntemplate", TEAL),
    ]
    rng = np.random.default_rng(427)
    for xi, (cond, label, color) in enumerate(arm_info):
        vals = np.array([r["val"] for r in e49_rbd["runs"] if r["cond"] == cond])
        jitter = rng.uniform(-0.10, 0.10, len(vals))
        ax.scatter(np.full(len(vals), xi) + jitter, vals, s=13, color=color,
                   alpha=0.58, edgecolor="none", zorder=2)
        mean = vals.mean()
        se = vals.std(ddof=1) / math.sqrt(len(vals))
        ax.errorbar(xi, mean, yerr=se, fmt="D", markersize=5.5, color=BLACK,
                    capsize=3, linewidth=0.9, zorder=4)
        ax.text(xi + 0.18, mean, f"{mean:.3f}", ha="left", va="center",
                fontsize=7.2, fontweight="bold",
                bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    ax.set_xticks(range(len(arm_info)), [x[1] for x in arm_info])
    ax.set_xlim(-0.4, 3.65)
    ax.set_ylim(0.15, 0.82)
    ax.set_ylabel("DMS Spearman ρ")
    ax.set_xlabel("RBD input arm")
    ax.set_title("500-mutant subset; 15 runs per arm", loc="left", pad=8)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    panel_label(ax, "C", x=-0.22, y=1.08)

    ax = fig.add_subplot(gs[1, 1])
    contrast_info = [
        ("sp:msa-no", "MSA − no", BLUE),
        ("sp:tpl-no", "template − no", ORANGE),
        ("sp:msatpl-no", "MSA+tpl − no", TEAL),
        ("sp:msatpl-tpl", "MSA+tpl − tpl", DARK_GREY),
    ]
    y = np.arange(len(contrast_info))
    for yi, (key, label, color) in zip(y, contrast_info):
        block = e427["paired"][key]
        vals = np.array(block["per_split"])
        ax.scatter(vals, np.full(len(vals), yi + 0.12), s=18, color=color, alpha=0.65,
                   edgecolor="none", zorder=2)
        mean = float(block["mean"])
        se = vals.std(ddof=1) / math.sqrt(len(vals))
        ax.errorbar(mean, yi - 0.04, xerr=se, fmt="D", markersize=5.5, color=BLACK,
                    capsize=3, linewidth=0.9, zorder=4)
        ax.text(mean, yi - 0.24, f"{mean:+.3f}", ha="center", va="bottom",
                fontsize=7.2, fontweight="bold")
    ax.axvline(0, color=BLACK, linewidth=0.7)
    ax.set_yticks(y, [x[1] for x in contrast_info])
    ax.set_ylim(len(contrast_info) - 0.4, -0.55)
    ax.set_xlim(-0.34, 0.46)
    ax.set_xticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlabel("Paired Δ Spearman ρ")
    ax.set_title("E427 block-CV; five folds", loc="left", pad=8)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    panel_label(ax, "D", x=-0.22, y=1.08)

    save_figure(fig, "Fig2_template_conflict")


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=str(DATA_DIR),
                    help="repo data directory containing results/ (default: %(default)s)")
    ap.add_argument("--output-dir", default=str(OUT),
                    help="figure output directory (default: %(default)s)")
    args = ap.parse_args()
    OUT = Path(args.output_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    paths = resolve_inputs(args.data_dir)
    e2, e49b, e428, e49_rbd, e427 = load_inputs(paths)
    make_fig1(e2, e428)
    make_fig2(e2, e49b, e428, e49_rbd, e427)
    for path in sorted(OUT.glob("Fig[12]_*")):
        if path.suffix in {".pdf", ".png", ".tif"}:
            print(f"{path.name}\t{path.stat().st_size}")


if __name__ == "__main__":
    main()