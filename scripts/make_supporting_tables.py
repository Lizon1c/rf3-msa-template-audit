#!/usr/bin/env python3
"""Build S1--S4 Table workbooks directly from canonical JSON results."""
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import xlsxwriter

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS = REPO_ROOT / "data" / "results"
OUT = REPO_ROOT / "manuscript" / "supporting_information"
TRAJECTORY_KEYS = {"traj_tr", "traj1", "traj2", "grad_traj"}
PREFERRED_COLUMNS = [
    "source",
    "target",
    "protocol",
    "condition",
    "cond",
    "probe",
    "arm",
    "normalization",
    "split",
    "fold",
    "init",
    "rep",
    "replicate",
    "cif",
    "val",
]


def load_json(relative_path):
    return json.loads((RESULTS / relative_path).read_text())


def flatten_scalars(value, prefix=""):
    row = {}
    for key, item in value.items():
        name = f"{prefix}.{key}" if prefix else key
        if key in TRAJECTORY_KEYS:
            continue
        if isinstance(item, dict):
            row.update(flatten_scalars(item, name))
        elif not isinstance(item, list):
            row[name] = item
    return row


def run_rows(relative_path, target, protocol, conditions=None):
    data = load_json(relative_path)
    rows = []
    for archived in data.get("runs", []):
        if conditions is not None and archived.get("cond") not in conditions:
            continue
        row = {
            "source": relative_path,
            "target": target,
            "protocol": protocol,
        }
        row.update(flatten_scalars(archived))
        rows.append(row)
    return rows


def summary_rows(relative_path, target, protocol, sections=("summary", "paired"), conditions=None):
    data = load_json(relative_path)
    rows = []

    def visit(value, path):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in TRAJECTORY_KEYS:
                    continue
                visit(item, path + [key])
        elif isinstance(value, list):
            if all(not isinstance(item, (dict, list)) for item in value):
                for index, item in enumerate(value, start=1):
                    rows.append({
                        "source": relative_path,
                        "target": target,
                        "protocol": protocol,
                        "metric": ".".join(path),
                        "index": index,
                        "value": item,
                    })
        else:
            rows.append({
                "source": relative_path,
                "target": target,
                "protocol": protocol,
                "metric": ".".join(path),
                "index": None,
                "value": value,
            })

    for section in sections:
        content = data.get(section, {})
        if conditions is not None:
            content = {key: value for key, value in content.items() if key in conditions}
        visit(content, [section])
    for key, value in data.items():
        if key not in {"runs", *sections, *TRAJECTORY_KEYS} and not isinstance(value, (dict, list)):
            rows.append({
                "source": relative_path,
                "target": target,
                "protocol": protocol,
                "metric": key,
                "index": None,
                "value": value,
            })
    return rows


def config_rows(relative_path, target, protocol):
    data = load_json(relative_path)
    rows = []
    for key, value in data.get("config", {}).items():
        rows.append({
            "source": relative_path,
            "target": target,
            "protocol": protocol,
            "field": key,
            "value": ", ".join(map(str, value)) if isinstance(value, list) else value,
        })
    if "epochs" in data:
        rows.append({
            "source": relative_path,
            "target": target,
            "protocol": protocol,
            "field": "epochs",
            "value": data["epochs"],
        })
    return rows


class TableWorkbook:
    def __init__(self, path, title, description):
        self.path = path
        self.workbook = xlsxwriter.Workbook(path, {"nan_inf_to_errors": True})
        self.workbook.set_properties({
            "title": title,
            "author": "Cheng Ding",
            "comments": "Generated from canonical archived JSON results; no rounded manuscript values used as inputs.",
        })
        self.header = self.workbook.add_format({
            "bold": True,
            "font_color": "#FFFFFF",
            "bg_color": "#164477",
            "border": 0,
            "text_wrap": True,
            "valign": "top",
        })
        self.text = self.workbook.add_format({"valign": "top"})
        self.wrap = self.workbook.add_format({"text_wrap": True, "valign": "top"})
        self.number = self.workbook.add_format({"num_format": "0.0000000000", "valign": "top"})
        self.integer = self.workbook.add_format({"num_format": "0", "valign": "top"})
        self.title = self.workbook.add_format({"bold": True, "font_size": 14, "font_color": "#164477"})
        self.note = self.workbook.add_format({"text_wrap": True, "valign": "top", "font_color": "#444444"})
        self.add_readme(title, description)

    def add_readme(self, title, description):
        sheet = self.workbook.add_worksheet("Read me")
        sheet.set_column("A:A", 24)
        sheet.set_column("B:B", 110)
        sheet.write(0, 0, title, self.title)
        sheet.merge_range(0, 0, 0, 1, title, self.title)
        sheet.write(2, 0, "Contents", self.header)
        sheet.write(2, 1, description, self.wrap)
        notes = [
            ("Source", "Generated by scripts/make_supporting_tables.py from data/results/*.json in the public repository."),
            ("Precision", "Numeric cells retain the archived floating-point values; the number format only controls display."),
            ("Comparisons", "Compare absolute Spearman rho values only within the same feature universe, split scheme, normalization, readout, and checkpoint-selection protocol."),
            ("Trajectory fields", "Long per-epoch training trajectories are not duplicated in these tables; run-level endpoint and selection metrics are retained. Canonical JSON files remain the full archive."),
            ("Manuscript", "An input-channel audit of RoseTTAFold3: MSA writing, mutation-effect prediction, and conformational-state selection."),
            ("Repository", "https://github.com/Lizon1c/rf3-msa-template-audit"),
            ("Data DOI", "https://doi.org/10.5281/zenodo.23859152"),
        ]
        for row_index, (key, value) in enumerate(notes, start=4):
            sheet.write(row_index, 0, key, self.header)
            sheet.write(row_index, 1, value, self.wrap)
        sheet.set_row(2, 48)
        sheet.freeze_panes(3, 0)

    def add_sheet(self, name, rows, columns=None):
        if not rows:
            raise ValueError(f"No rows supplied for {self.path.name}:{name}")
        sheet = self.workbook.add_worksheet(name)
        discovered = []
        for row in rows:
            for key in row:
                if key not in discovered:
                    discovered.append(key)
        if columns is None:
            columns = [key for key in PREFERRED_COLUMNS if key in discovered]
            columns += [key for key in discovered if key not in columns]
        else:
            columns += [key for key in discovered if key not in columns]
        for column_index, key in enumerate(columns):
            sheet.write(0, column_index, key, self.header)
        for row_index, row in enumerate(rows, start=1):
            for column_index, key in enumerate(columns):
                value = row.get(key)
                if value is None:
                    sheet.write_blank(row_index, column_index, None, self.text)
                elif isinstance(value, bool):
                    sheet.write_boolean(row_index, column_index, value, self.text)
                elif isinstance(value, int):
                    sheet.write_number(row_index, column_index, value, self.integer)
                elif isinstance(value, float):
                    sheet.write_number(row_index, column_index, value, self.number)
                else:
                    sheet.write_string(row_index, column_index, str(value), self.text)
        sheet.freeze_panes(1, 0)
        sheet.autofilter(0, 0, len(rows), len(columns) - 1)
        for index, key in enumerate(columns):
            width = min(42, max(11, len(key) + 2))
            if key in {"source", "protocol", "metric"}:
                width = 42
            elif key in {"target", "condition", "cond", "cif"}:
                width = 24
            sheet.set_column(index, index, width)

    def close(self):
        self.workbook.close()
        print(f"{self.path.relative_to(REPO_ROOT)}\t{self.path.stat().st_size}")


def make_s1():
    book = TableWorkbook(
        OUT / "S1_Table.xlsx",
        "S1 Table. DMS de-biasing, block-CV, and lightweight probes",
        "Canonical, dual-validation, normalization, E427 block-CV, and lightweight-probe summaries and run-level values. RBD template-arm block-CV is isolated in its own worksheet.",
    )
    canonical_specs = [
        ("rbd_canonical/four_condition_v3.json", "RBD", "Canonical full-set", {"sp_no", "sp_msa"}),
        ("rbd_canonical/four_condition_A.json", "RBD", "Canonical full-set independent replicate", {"sp_no", "sp_msa"}),
        ("casp3/e4_msa_analysis.json", "CASP3", "Canonical full-set", {"sp_no", "sp_msa"}),
        ("ci/e4_msa_analysis.json", "CI", "Canonical full-set; not fully de-biased", {"sp_no", "sp_msa"}),
    ]
    rows = []
    summaries = []
    for path, target, protocol, conditions in canonical_specs:
        rows += run_rows(path, target, protocol, conditions)
        summaries += summary_rows(path, target, protocol, conditions=conditions)
    book.add_sheet("Canonical runs", rows)
    book.add_sheet("Canonical summaries", summaries)

    dual_specs = [
        ("ab_validation/v1_rbd_full.json", "RBD"),
        ("ab_validation/v1_casp3_main.json", "CASP3"),
    ]
    rows = []
    summaries = []
    for path, target in dual_specs:
        rows += run_rows(path, target, "Dual-validation 60/20/20 position split", {"sp_no", "sp_msa"})
        summaries += summary_rows(path, target, "Dual-validation 60/20/20 position split", conditions={"sp_no", "sp_msa"})
    book.add_sheet("Dual-validation runs", rows)
    book.add_sheet("Dual-validation summaries", summaries)

    path = "ab_validation/v2_std_full_sp.json"
    book.add_sheet("Train-only standardization", run_rows(path, "RBD", "Full-set overall versus strict train-only standardization"))
    book.add_sheet("Standardization summaries", summary_rows(path, "RBD", "Full-set overall versus strict train-only standardization"))
    path = "ab_validation/v2_std_subset.json"
    book.add_sheet("Subset normalization", run_rows(path, "RBD 500-mutant subset", "Subset normalization sensitivity analysis"))

    block_paths = [
        ("e427/results_A1.json", "RBD full-set"),
        ("e427/results_A2.json", "RBD 500-mutant subset"),
        ("e427/results_A3.json", "CASP3 full-set"),
    ]
    rows = []
    summaries = []
    configs = []
    for path, target in block_paths:
        protocol = "E427 five-fold contiguous sequence-block CV"
        rows += run_rows(path, target, protocol)
        summaries += summary_rows(path, target, protocol)
        configs += config_rows(path, target, protocol)
    book.add_sheet("E427 block-CV runs", rows)
    book.add_sheet("E427 block-CV summaries", summaries)
    book.add_sheet("E427 configurations", configs)

    light_paths = [
        ("e427/results_B_rbd.json", "RBD full-set"),
        ("e427/results_B_casp3.json", "CASP3 full-set"),
        ("e427/results_B_ci.json", "CI full-set"),
    ]
    rows = []
    summaries = []
    for path, target in light_paths:
        protocol = "E427 lightweight probes on canonical random position splits"
        rows += run_rows(path, target, protocol)
        summaries += summary_rows(path, target, protocol)
    book.add_sheet("Light-probe runs", rows)
    book.add_sheet("Light-probe summaries", summaries)
    book.close()


def make_s2():
    book = TableWorkbook(
        OUT / "S2_Table.xlsx",
        "S2 Table. Two-state distogram preferences and template steering",
        "Thirteen range-matched two-state entries, five template-steering targets, and the separate GA88/GB88 E49b scoring flow. Nested NLL, MAE, entropy, and pair-count fields are flattened into named columns.",
    )
    path = "e2/e2_rangematched_results.json"
    data = load_json(path)
    input_root = REPO_ROOT / "data" / "inputs" / "e2_targets"
    definitions = {
        row["id"]: row
        for row in json.loads((input_root / "e2_targets.json").read_text())["targets"]
    }
    rows = []
    for archived in data:
        row = {"source": path, "protocol": "Range-matched distogram preference"}
        row.update(flatten_scalars(archived))
        definition = definitions[archived["id"]]
        row["input_definition"] = "data/inputs/e2_targets/e2_targets.json"
        row["msa_rows_including_query"] = sum(
            line.startswith(">")
            for line in (input_root / "msas" / f"{archived['id']}.a3m").read_text().splitlines()
        )
        for state, prefix in (("state_a", "state_A"), ("state_b", "state_B")):
            for key in ("pdb", "chain", "method", "resolution"):
                row[f"{prefix}_{key}"] = definition[state][key]
        rows.append(row)
    book.add_sheet("13-entry panel", rows)
    template_rows = [row for row in rows if "pref_tplA" in row]
    if len(template_rows) != 5:
        raise ValueError(f"Expected five distogram template targets, found {len(template_rows)}")
    book.add_sheet("Template steering", template_rows)

    path = "e2/e49b_results.json"
    data = load_json(path)
    rows = []
    for target, archived in data.items():
        if not isinstance(archived, dict):
            continue
        row = {"source": path, "target": target, "protocol": "E49b separate residue-universe scoring flow"}
        row.update(flatten_scalars(archived))
        rows.append(row)
    book.add_sheet("GA88-GB88 E49b", rows)
    book.close()


def make_s3():
    book = TableWorkbook(
        OUT / "S3_Table.xlsx",
        "S3 Table. Protocol E428 explicit 3D coordinate sampling",
        "Thirty target-by-arm cell summaries and all 480 sample-level RMSD records. Positive delta RMSD = RMSD_B - RMSD_A indicates coordinates closer to State A.",
    )
    path = "e428/results.json"
    archived_rows = load_json(path)
    if len(archived_rows) != 480:
        raise ValueError(f"Expected 480 structures, found {len(archived_rows)}")
    cells = defaultdict(list)
    sample_rows = []
    for archived in archived_rows:
        row = {"source": path, "protocol": "E428 explicit coordinate sampling"}
        row.update(flatten_scalars(archived))
        row["delta_RMSD"] = archived["rmsdB"] - archived["rmsdA"]
        row["closer_to_A_recomputed"] = archived["rmsdA"] < archived["rmsdB"]
        if "closer_to" in archived:
            expected = "A" if row["closer_to_A_recomputed"] else "B"
            if archived["closer_to"] != expected:
                raise ValueError(f"Closer-state mismatch in {archived}")
        sample_rows.append(row)
        cells[(archived["tid"], archived["arm"])].append(archived)
    if len(cells) != 30 or any(len(rows) != 16 for rows in cells.values()):
        raise ValueError("E428 is not a complete 30-cell by 16-structure archive")

    cell_rows = []
    for (target, arm), rows in cells.items():
        a = [row["rmsdA"] for row in rows]
        b = [row["rmsdB"] for row in rows]
        delta = [row["rmsdB"] - row["rmsdA"] for row in rows]
        closer_a = sum(row["rmsdA"] < row["rmsdB"] for row in rows)
        cell_rows.append({
            "source": path,
            "target": target,
            "arm": arm,
            "n_structures": len(rows),
            "RMSD_A_mean": statistics.mean(a),
            "RMSD_B_mean": statistics.mean(b),
            "delta_RMSD_mean": statistics.mean(delta),
            "delta_RMSD_sd": statistics.stdev(delta),
            "n_closer_to_A": closer_a,
            "fraction_closer_to_A": closer_a / len(rows),
            "percent_closer_to_A": 100 * closer_a / len(rows),
        })
    book.add_sheet("30-cell summaries", cell_rows)
    book.add_sheet("480 structure-level RMSD", sample_rows)
    book.close()


def make_s4():
    book = TableWorkbook(
        OUT / "S4_Table.xlsx",
        "S4 Table. Fusion, interventions, distillation, noise, and interface statistics",
        "Three-stream fusion, RBD learning-rate and input-noise interventions, E429 same-initialization replicates, D1 reconstruction, and CASP3 paired/unpaired interface statistics.",
    )
    fusion_specs = [
        ("rbd_canonical/four_condition_v3.json", "RBD", "Canonical two-stream conditions"),
        ("rbd_canonical/matrix_complete_v3.json", "RBD", "Three-stream and ESM conditions"),
        ("casp3/e4_msa_analysis.json", "CASP3", "Canonical single-stream pair conditions"),
        ("casp3/e5_matrix.json", "CASP3", "Two-stream and three-stream fusion"),
        ("ci/e4_msa_analysis.json", "CI", "Canonical single-stream pair conditions; not fully de-biased"),
        ("ci/e5_matrix.json", "CI", "Two-stream and three-stream fusion; not fully de-biased"),
    ]
    rows = []
    summaries = []
    for path, target, protocol in fusion_specs:
        rows += run_rows(path, target, protocol)
        summaries += summary_rows(path, target, protocol)
    book.add_sheet("Fusion runs", rows)
    book.add_sheet("Fusion summaries", summaries)

    for path, sheet, protocol in [
        ("ab_validation/v4_pace_fenc.json", "Learning-rate intervention", "RBD FAESM encoder learning-rate scaling"),
        ("ab_validation/v3_esm_noise.json", "Input-noise intervention", "RBD FAESM Gaussian input-noise intervention"),
        ("rbd_canonical/results_gradprobe.json", "Gradient diagnostics", "RBD late-training encoder gradient diagnostics"),
    ]:
        book.add_sheet(sheet, run_rows(path, "RBD", protocol))
        intervention_summaries = summary_rows(path, "RBD", protocol)
        if intervention_summaries:
            book.add_sheet(f"{sheet[:22]} summary", intervention_summaries)

    rows = []
    summaries = []
    for process in ("A", "B"):
        path = f"e429/results_{process}.json"
        protocol = f"E429 same split and initialization; process {process}"
        rows += run_rows(path, "RBD x_msa_faesm", protocol)
        summaries += summary_rows(path, "RBD x_msa_faesm", protocol)
    book.add_sheet("E429 16 replicates", rows)
    book.add_sheet("E429 summaries", summaries)

    d1_paths = [
        "d1/results.json",
        "d1/deep_both.json",
        "d1/deep_z.json",
        "d1/deep_f.json",
        "d1/xattn_cka.json",
        "d1/xattn_cka_deep.json",
        "d1/onehot.json",
    ]
    rows = []
    summaries = []
    for path in d1_paths:
        rows += run_rows(path, "RBD", "D1 pair-representation reconstruction")
        summaries += summary_rows(path, "RBD", "D1 pair-representation reconstruction")
    book.add_sheet("D1 reconstruction runs", rows)
    if summaries:
        book.add_sheet("D1 reconstruction summaries", summaries)
    book.add_sheet("D1 functional recovery", run_rows("d1/functional.json", "RBD", "D1 downstream DMS functional recovery"))

    rows = []
    summaries = []
    for condition in ("paired", "unpaired"):
        path = f"casp3_interface/{condition}_results.json"
        data = load_json(path)
        protocol = f"CASP3 homodimer interface; {condition} MSA"
        for archived in data["per_mutant"]:
            row = {"source": path, "target": "CASP3", "protocol": protocol}
            row.update(flatten_scalars(archived))
            rows.append(row)
        summaries += summary_rows(path, "CASP3", protocol)
    book.add_sheet("CASP3 interface per mutant", rows)
    book.add_sheet("CASP3 interface summaries", summaries)
    book.close()


def verify_headline_inputs():
    checks = [
        ("e427/results_A1.json", "sp:msa-no", 0.256),
        ("e427/results_A3.json", "sp:msa-no", 0.264),
        ("e427/results_A2.json", "sp:tpl-no", -0.129),
        ("e427/results_A2.json", "sp:msatpl-tpl", 0.003),
        ("e427/results_B_rbd.json", "linear:msa-no", 0.030),
        ("e427/results_B_rbd.json", "mlp:msa-no", 0.157),
    ]
    for path, contrast, expected in checks:
        value = load_json(path)["paired"][contrast]["mean"]
        if not math.isclose(value, expected, abs_tol=0.0005):
            raise ValueError(f"Headline mismatch for {path}:{contrast}: {value} versus {expected}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    verify_headline_inputs()
    make_s1()
    make_s2()
    make_s3()
    make_s4()


if __name__ == "__main__":
    main()
