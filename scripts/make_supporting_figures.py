#!/usr/bin/env python3
"""Generate the two supporting figures from canonical archived results."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from make_figures import (
    BLACK,
    BLUE,
    CYAN,
    DARK_GREY,
    GREY,
    ORANGE,
    TEAL,
    panel_label,
    save_figure,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS = REPO_ROOT / "data" / "results"


def load_json(relative_path):
    return json.loads((RESULTS / relative_path).read_text())


def mean_se(values):
    values = np.asarray(values, dtype=float)
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(len(values)))


def add_bar_labels(ax, bars, fmt="{:.3f}", dy=0.02):
    y0, y1 = ax.get_ylim()
    for bar in bars:
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + dy * (y1 - y0),
            fmt.format(bar.get_height()),
            ha="center",
            va="bottom",
            fontsize=7,
        )


def add_run_panel(ax, runs_by_arm, labels, colors, ylim):
    rng = np.random.default_rng(427)
    for index, arm in enumerate(labels):
        values = np.array([row["val"] for row in runs_by_arm[arm]], dtype=float)
        jitter = rng.uniform(-0.09, 0.09, len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s=10,
            color=colors[arm],
            alpha=0.52,
            edgecolor="none",
            zorder=2,
        )
        mean, se = mean_se(values)
        ax.errorbar(
            index,
            mean,
            yerr=se,
            fmt="D",
            markersize=4.8,
            color=BLACK,
            capsize=2.5,
            linewidth=0.8,
            zorder=4,
        )
        label_offset = 0.040 if index % 2 == 0 else -0.050
        ax.text(
            index,
            mean + label_offset,
            f"{mean:.3f}",
            ha="center",
            va="bottom" if label_offset > 0 else "top",
            fontsize=6.5,
            bbox={"boxstyle": "round,pad=0.08", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        )
    ax.set_xticks(range(len(labels)), labels, rotation=32, ha="right")
    ax.set_ylim(*ylim)
    ax.set_ylabel("DMS Spearman ρ")
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def make_s1_fig():
    rbd = load_json("e1_geometry/results_fullset.json")["arms"]
    casp3 = load_json("casp3/e4_msa_analysis.json")["geo"]
    ci = load_json("ci/e4_msa_analysis.json")["geo"]
    no_values = [rbd["noMSA_ref"]["eff90"], casp3["eff90_no"], ci["eff90_no"]]
    msa_values = [rbd["zii_e1_full15"]["eff90"], casp3["eff90_msa"], ci["eff90_msa"]]
    cka_values = [rbd["zii_e1_full15"]["cka_vs_noMSA"], casp3["cka"], ci["cka"]]
    dcos_values = [
        rbd["zii_e1_full15"]["delta_cos_vs_noMSA_median"],
        casp3["dcos_median"],
        ci["dcos_median"],
    ]
    labels = ["RBD", "CASP3", "CI"]

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75))
    ax = axes[0]
    x = np.arange(3)
    width = 0.34
    bars_no = ax.bar(x - width / 2, no_values, width, color=DARK_GREY, label="no MSA")
    bars_msa = ax.bar(x + width / 2, msa_values, width, color=BLUE, label="MSA")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 390)
    ax.set_ylabel("Effective rank (90% variance)")
    ax.set_title("Effective rank", loc="left", pad=8)
    ax.legend(frameon=False, loc="upper right")
    for bars in (bars_no, bars_msa):
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 8,
                f"{int(bar.get_height())}",
                ha="center",
                va="bottom",
                fontsize=7,
            )
    for index, (before, after) in enumerate(zip(no_values, msa_values)):
        ax.text(
            index,
            max(before, after) + 42,
            f"{100 * (after / before - 1):+.0f}%",
            ha="center",
            va="bottom",
            fontsize=7.5,
            fontweight="bold",
        )
    panel_label(ax, "A", x=-0.18, y=1.08)

    for ax, values, title, ylabel, panel, label_x in [
        (axes[1], cka_values, "Linear CKA", "Linear CKA", "B", -0.18),
        (axes[2], dcos_values, "Mutant directions", "Median dCos", "C", -0.18),
    ]:
        bars = ax.bar(x, values, width=0.58, color=[BLUE, TEAL, ORANGE])
        ax.set_xticks(x, labels)
        ax.set_ylim(0, 1.12)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", pad=8)
        add_bar_labels(ax, bars, dy=0.018)
        panel_label(ax, panel, x=label_x, y=1.08)

    fig.subplots_adjust(wspace=0.40)
    save_figure(fig, "S1_Fig")


def make_s2_fig():
    geometry = load_json("e1_geometry/results_results.json")["arms"]
    e7_geometry = load_json("e1_geometry/results_e7.json")["arms"]
    e1_runs = load_json("e1_supervised/results_arms.json")["runs"]
    endpoint_runs = load_json("e1_supervised/results_no_vs_full15.json")["runs"]
    e7_runs = load_json("e1_supervised/results_e7.json")["runs"]
    runs_by_arm = {}
    for arm in ["noMSA", "depth01", "random", "rowshuffle", "consensus", "colshuffle", "depth04", "depth08", "full15"]:
        if arm == "noMSA":
            runs_by_arm[arm] = [row for row in endpoint_runs if row["cond"] == "no"]
        elif arm == "full15":
            runs_by_arm[arm] = [row for row in endpoint_runs if row["cond"] == "full15"]
        else:
            runs_by_arm[arm] = [row for row in e1_runs if row["cond"] == f"zii_e1_{arm}"]
    for arm in ["k1", "k2", "k3"]:
        runs_by_arm[arm] = [row for row in e7_runs if row["cond"] == f"zii_e7_{arm}"]

    e1_arms = ["noMSA", "depth01", "random", "rowshuffle", "consensus", "colshuffle", "depth04", "depth08", "full15"]
    geometry_keys = {arm: ("noMSA_ref" if arm == "noMSA" else f"zii_e1_{arm}") for arm in e1_arms}
    e1_colors = {
        "noMSA": DARK_GREY,
        "depth01": GREY,
        "random": ORANGE,
        "rowshuffle": ORANGE,
        "consensus": CYAN,
        "colshuffle": TEAL,
        "depth04": BLUE,
        "depth08": BLUE,
        "full15": BLUE,
    }
    e7_arms = ["noMSA", "k1", "k2", "k3", "full15"]
    e7_colors = {"noMSA": DARK_GREY, "k1": CYAN, "k2": CYAN, "k3": CYAN, "full15": BLUE}

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
    ax = axes[0, 0]
    x = np.arange(len(e1_arms))
    ranks = [geometry[geometry_keys[arm]]["eff90"] for arm in e1_arms]
    bars = ax.bar(x, ranks, color=[e1_colors[arm] for arm in e1_arms], width=0.72)
    ax.set_xticks(x, e1_arms, rotation=32, ha="right")
    ax.set_ylim(0, 213)
    ax.set_ylabel("Effective rank (90% variance)")
    ax.set_title("E1: alignment dose × content", loc="left", pad=8)
    add_bar_labels(ax, bars, fmt="{:.0f}", dy=0.012)
    panel_label(ax, "A", x=-0.20, y=1.08)

    ax = axes[0, 1]
    add_run_panel(ax, runs_by_arm, e1_arms, e1_colors, (0.30, 0.79))
    ax.set_title("E1: mutation-effect prediction", loc="left", pad=8)
    panel_label(ax, "B", x=-0.18, y=1.08)

    ax = axes[1, 0]
    x = np.arange(len(e7_arms))
    e7_keys = {arm: ("noMSA_ref" if arm == "noMSA" else ("zii_e1_full15" if arm == "full15" else f"zii_e7_{arm}")) for arm in e7_arms}
    ranks = [e7_geometry[e7_keys[arm]]["eff90"] for arm in e7_arms]
    bars = ax.bar(x, ranks, color=[e7_colors[arm] for arm in e7_arms], width=0.62)
    ax.set_xticks(x, ["no MSA", "K=1", "K=2", "K=3", "K=10"])
    ax.set_ylim(0, 213)
    ax.set_ylabel("Effective rank (90% variance)")
    ax.set_title("E7: MSA recycle gating", loc="left", pad=8)
    add_bar_labels(ax, bars, fmt="{:.0f}", dy=0.012)
    panel_label(ax, "C", x=-0.20, y=1.08)

    ax = axes[1, 1]
    add_run_panel(ax, runs_by_arm, e7_arms, e7_colors, (0.30, 0.79))
    ax.set_xticks(range(len(e7_arms)), ["no MSA", "K=1", "K=2", "K=3", "K=10"], rotation=0, ha="center")
    ax.set_title("E7: mutation-effect prediction", loc="left", pad=8)
    panel_label(ax, "D", x=-0.18, y=1.08)

    fig.subplots_adjust(wspace=0.33, hspace=0.60)
    save_figure(fig, "S2_Fig")


def main():
    make_s1_fig()
    make_s2_fig()
    for name in ("S1_Fig", "S2_Fig"):
        for suffix in (".pdf", ".png", ".tif"):
            path = REPO_ROOT / "figures" / f"{name}{suffix}"
            print(f"{path.relative_to(REPO_ROOT)}\t{path.stat().st_size}")


if __name__ == "__main__":
    main()
