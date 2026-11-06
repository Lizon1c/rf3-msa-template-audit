#!/usr/bin/env python3
"""Validate the assembled PLOS submission files and their numerical contents."""
import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile

from PIL import Image
from pypdf import PdfReader

from assemble_submission import SOURCE_FILES, UPLOAD_FILES

REPO_ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = REPO_ROOT / "submission"
RESULTS = REPO_ROOT / "data" / "results"
CHECKS = []
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def check(name):
    def decorate(function):
        CHECKS.append((name, function))
        return function
    return decorate


def pdf_text(name):
    return "\n".join(page.extract_text() or "" for page in PdfReader(SUBMISSION / name).pages)


def load_json(relative_path):
    return json.loads((RESULTS / relative_path).read_text())


def xlsx_sheet_names(path):
    with ZipFile(path) as archive:
        root = ET.fromstring(archive.read("xl/workbook.xml"))
        return [sheet.attrib["name"] for sheet in root.find(f"{{{MAIN_NS}}}sheets")]


def xlsx_rows(path, sheet_name):
    with ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in shared_root.findall(f"{{{MAIN_NS}}}si")]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relation_id = None
        for sheet in workbook.find(f"{{{MAIN_NS}}}sheets"):
            if sheet.attrib["name"] == sheet_name:
                relation_id = sheet.attrib[f"{{{REL_NS}}}id"]
                break
        if relation_id is None:
            raise ValueError(f"Missing worksheet {sheet_name} in {path.name}")
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = None
        for relation in relations:
            if relation.attrib["Id"] == relation_id:
                target = relation.attrib["Target"]
                break
        if target is None:
            raise ValueError(f"Missing worksheet relationship {relation_id}")
        target_path = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
        worksheet = ET.fromstring(archive.read(target_path))
        raw_rows = []
        for row in worksheet.findall(f".//{{{MAIN_NS}}}sheetData/{{{MAIN_NS}}}row"):
            values = []
            for cell in row.findall(f"{{{MAIN_NS}}}c"):
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{MAIN_NS}}}v")
                if cell_type == "inlineStr":
                    inline = cell.find(f"{{{MAIN_NS}}}is")
                    value = "".join(inline.itertext()) if inline is not None else ""
                elif value_node is None:
                    value = ""
                elif cell_type == "s":
                    value = shared[int(value_node.text)]
                else:
                    value = value_node.text
                values.append(value)
            raw_rows.append(values)
    if not raw_rows:
        raise ValueError(f"Empty worksheet {path.name}:{sheet_name}")
    headers = raw_rows[0]
    rows = []
    for raw in raw_rows[1:]:
        rows.append({header: raw[index] if index < len(raw) else "" for index, header in enumerate(headers)})
    return rows


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@check("11 upload files and auxiliary guides present")
def required_files():
    required = set(UPLOAD_FILES) | {"README.md", "Submission_Metadata.md"}
    missing = [name for name in required if not (SUBMISSION / name).is_file()]
    assert not missing, missing
    oversized = [name for name in UPLOAD_FILES if (SUBMISSION / name).stat().st_size >= 10_000_000]
    assert not oversized, f"Files over 10 MB: {oversized}"
    return f"{len(UPLOAD_FILES)} upload files; all below 10 MB"


@check("LaTeX and vector source copies match originals")
def source_copies():
    for destination, source in SOURCE_FILES.items():
        assert sha256(SUBMISSION / destination) == sha256(REPO_ROOT / source), destination
    return f"{len(SOURCE_FILES)} byte-identical source files"


@check("Manuscript PDF content and page count")
def manuscript_pdf():
    reader = PdfReader(SUBMISSION / "Manuscript.pdf")
    assert len(reader.pages) == 14, len(reader.pages)
    raw_text = "\n".join(page.extract_text() or "" for page in reader.pages)
    text = " ".join(raw_text.split())
    compact_text = "".join(raw_text.split())
    for required in (
        "An input-channel audit of RoseTTAFold3",
        "Supporting Information Captions",
        "Corley N, Krishna R, Bauer MS, Thompson TR, Ahern W, Kazman MW",
        "Accelerating biomolecular modeling with AtomWorks and RF3",
    ):
        assert required in text, required
    assert any(marker in text for marker in ("S4 Table.", "S4 T able."))
    assert "10.5281/zenodo.23859152" in compact_text
    for obsolete in (
        "Corley N, Butterfoss GL",
        "AtomWorks and RosettaFold-3",
        "No external funding was received for this work",
        "Supporting Information Tables",
        "fig1_geometry_collapse.png",
    ):
        assert obsolete not in text, obsolete
    return "14 pages; corrected references; captions after bibliography; no embedded SI tables"


@check("Cover Letter PDF content")
def cover_letter_pdf():
    reader = PdfReader(SUBMISSION / "Cover_Letter.pdf")
    assert len(reader.pages) == 1, len(reader.pages)
    raw_text = "\n".join(page.extract_text() or "" for page in reader.pages)
    text = " ".join(raw_text.split())
    compact_text = "".join(raw_text.split())
    for required in ("Research Article", "Biological Macromolecules", "480 structures", "rf3-msa-template-audit"):
        assert required in text, required
    assert "10.5281/zenodo.23859152" in compact_text
    return "1 page; journal fit, significance, declarations, and repositories present"


@check("S1 Text PDF content")
def supporting_text_pdf():
    reader = PdfReader(SUBMISSION / "S1_Text.pdf")
    assert len(reader.pages) == 9, len(reader.pages)
    text = " ".join("\n".join(page.extract_text() or "" for page in reader.pages).split())
    for required in ("S1 Text", "Protocol map", "E429 same-initialization noise calibration", "Analysis-artifact register"):
        assert required in text, required
    return "9 pages; extended methods and artifact register present"


@check("Main and supporting TIFF specifications")
def tiff_specifications():
    details = []
    for name in ("Fig1.tif", "Fig2.tif", "S1_Fig.tif", "S2_Fig.tif"):
        with Image.open(SUBMISSION / name) as image:
            dpi = image.info.get("dpi", (0, 0))
            assert image.mode == "RGB", (name, image.mode)
            assert all(abs(value - 600) < 0.1 for value in dpi), (name, dpi)
            assert image.info.get("compression") == "tiff_lzw", (name, image.info.get("compression"))
            details.append(f"{name} {image.size[0]}x{image.size[1]}")
    return "; ".join(details)


@check("Workbook containers and worksheet names")
def workbook_containers():
    expected = {
        "S1_Table.xlsx": 13,
        "S2_Table.xlsx": 4,
        "S3_Table.xlsx": 3,
        "S4_Table.xlsx": 14,
    }
    for name, count in expected.items():
        path = SUBMISSION / name
        with ZipFile(path) as archive:
            assert archive.testzip() is None, name
        sheets = xlsx_sheet_names(path)
        assert len(sheets) == count, (name, sheets)
        assert sheets[0] == "Read me", (name, sheets)
    return "4 valid XLSX containers; 34 worksheets"


@check("S1 run counts and E427 archived values")
def s1_values():
    expected_counts = {
        "Canonical runs": 120,
        "Dual-validation runs": 60,
        "Train-only standardization": 60,
        "Subset normalization": 240,
        "E427 block-CV runs": 120,
        "Light-probe runs": 180,
    }
    for sheet, count in expected_counts.items():
        assert len(xlsx_rows(SUBMISSION / "S1_Table.xlsx", sheet)) == count, sheet
    for sheet, directory in (("E427 block-CV runs", "e427"), ("Light-probe runs", "e427")):
        rows = xlsx_rows(SUBMISSION / "S1_Table.xlsx", sheet)
        grouped = defaultdict(list)
        for row in rows:
            grouped[row["source"]].append(row)
        for source, archived_rows in grouped.items():
            archived = load_json(source)["runs"]
            assert len(archived_rows) == len(archived), source
            for row, original in zip(archived_rows, archived):
                assert math.isclose(float(row["val"]), original["val"], abs_tol=1e-12), source
    return "780 run-level rows; all 300 E427 values match canonical JSON"


@check("S2 panel, PDB, and MSA-row metadata")
def s2_values():
    rows = xlsx_rows(SUBMISSION / "S2_Table.xlsx", "13-entry panel")
    archived = load_json("e2/e2_rangematched_results.json")
    assert len(rows) == len(archived) == 13
    for row, original in zip(rows, archived):
        assert row["id"] == original["id"]
        assert math.isclose(float(row["pref_nomsa"]), original["pref_nomsa"], abs_tol=1e-12)
        assert math.isclose(float(row["pref_msa"]), original["pref_msa"], abs_tol=1e-12)
        msa_path = REPO_ROOT / "data" / "inputs" / "e2_targets" / "msas" / f"{original['id']}.a3m"
        expected_depth = sum(line.startswith(">") for line in msa_path.read_text().splitlines())
        assert int(row["msa_rows_including_query"]) == expected_depth
    assert len(xlsx_rows(SUBMISSION / "S2_Table.xlsx", "Template steering")) == 5
    gagb = xlsx_rows(SUBMISSION / "S2_Table.xlsx", "GA88-GB88 E49b")
    assert [row["target"] for row in gagb] == ["GA88", "GB88"]
    return "13 entries + 5 template targets + 2 designed-pair records verified"


@check("S3 all 30 cells and 480 structures")
def s3_values():
    summaries = xlsx_rows(SUBMISSION / "S3_Table.xlsx", "30-cell summaries")
    samples = xlsx_rows(SUBMISSION / "S3_Table.xlsx", "480 structure-level RMSD")
    archived = load_json("e428/results.json")
    assert len(summaries) == 30
    assert len(samples) == len(archived) == 480
    for row, original in zip(samples, archived):
        assert row["cif"] == original["cif"]
        expected = original["rmsdB"] - original["rmsdA"]
        assert math.isclose(float(row["delta_RMSD"]), expected, abs_tol=1e-9)
    counts = Counter((row["tid"], row["arm"]) for row in samples)
    assert len(counts) == 30 and set(counts.values()) == {16}
    return "30 cells x 16 structures; every sample delta RMSD recomputed"


@check("S4 noise, distillation, and interface rows")
def s4_values():
    noise = xlsx_rows(SUBMISSION / "S4_Table.xlsx", "E429 16 replicates")
    archived_noise = load_json("e429/results_A.json")["runs"] + load_json("e429/results_B.json")["runs"]
    assert len(noise) == len(archived_noise) == 16
    for row, original in zip(noise, archived_noise):
        assert math.isclose(float(row["val"]), original["val"], abs_tol=1e-12)
    assert len(xlsx_rows(SUBMISSION / "S4_Table.xlsx", "D1 functional recovery")) == 6
    interfaces = xlsx_rows(SUBMISSION / "S4_Table.xlsx", "CASP3 interface per mutant")
    assert len(interfaces) == 492
    counts = Counter(row["protocol"] for row in interfaces)
    assert sorted(counts.values()) == [100, 392], counts
    return "16 E429 + 6 functional-recovery + 492 interface records verified"


@check("LaTeX logs free of layout and reference failures")
def latex_logs():
    for stem in ("rf3_msa_template_audit_manuscript", "S1_text", "cover_letter"):
        log = (REPO_ROOT / "manuscript" / f"{stem}.log").read_text(errors="replace")
        for forbidden in ("Overfull", "Undefined control sequence", "Label(s) may have changed", "Reference undefined", "Citation undefined"):
            assert forbidden not in log, (stem, forbidden)
    return "No Overfull, undefined-control-sequence, or unresolved-reference diagnostics"


@check("SHA256SUMS complete and correct")
def checksums():
    lines = (SUBMISSION / "SHA256SUMS.txt").read_text().splitlines()
    expected_files = {
        path.relative_to(SUBMISSION).as_posix(): path
        for path in SUBMISSION.rglob("*")
        if path.is_file() and path.name not in {"SHA256SUMS.txt", "validation.json"}
    }
    recorded = {}
    for line in lines:
        digest, name = line.split("  ", 1)
        recorded[name] = digest
    assert set(recorded) == set(expected_files), (set(expected_files) - set(recorded), set(recorded) - set(expected_files))
    for name, digest in recorded.items():
        assert sha256(expected_files[name]) == digest, name
    return f"{len(recorded)} hashes match"


def main():
    report = {"status": "PASS", "checks": []}
    failures = 0
    for name, function in CHECKS:
        try:
            detail = function()
            report["checks"].append({"name": name, "status": "PASS", "detail": detail})
            print(f"PASS  {name:<48} {detail}")
        except Exception as error:
            failures += 1
            report["checks"].append({"name": name, "status": "FAIL", "error": str(error)})
            print(f"FAIL  {name:<48} {error}")
    report["status"] = "PASS" if failures == 0 else "FAIL"
    report["passed"] = len(CHECKS) - failures
    report["total"] = len(CHECKS)
    (SUBMISSION / "validation.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"\n{report['status']} — {report['passed']}/{report['total']} checks")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
