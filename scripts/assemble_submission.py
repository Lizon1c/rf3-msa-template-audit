#!/usr/bin/env python3
"""Assemble the PLOS Computational Biology initial-submission directory."""
import hashlib
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = REPO_ROOT / "submission"

UPLOAD_FILES = {
    "Manuscript.pdf": "manuscript/rf3_msa_template_audit_manuscript.pdf",
    "Cover_Letter.pdf": "manuscript/cover_letter.pdf",
    "Fig1.tif": "figures/Fig1_conformational_selection.tif",
    "Fig2.tif": "figures/Fig2_template_conflict.tif",
    "S1_Text.pdf": "manuscript/S1_text.pdf",
    "S1_Fig.tif": "figures/S1_Fig.tif",
    "S2_Fig.tif": "figures/S2_Fig.tif",
    "S1_Table.xlsx": "manuscript/supporting_information/S1_Table.xlsx",
    "S2_Table.xlsx": "manuscript/supporting_information/S2_Table.xlsx",
    "S3_Table.xlsx": "manuscript/supporting_information/S3_Table.xlsx",
    "S4_Table.xlsx": "manuscript/supporting_information/S4_Table.xlsx",
}

SOURCE_FILES = {
    "source/rf3_msa_template_audit_manuscript.tex": "manuscript/rf3_msa_template_audit_manuscript.tex",
    "source/S1_text.tex": "manuscript/S1_text.tex",
    "source/cover_letter.tex": "manuscript/cover_letter.tex",
    "source/supporting_information_captions.tex": "manuscript/supporting_information_captions.tex",
    "source/figures/Fig1_conformational_selection.pdf": "figures/Fig1_conformational_selection.pdf",
    "source/figures/Fig2_template_conflict.pdf": "figures/Fig2_template_conflict.pdf",
    "source/figures/S1_Fig.pdf": "figures/S1_Fig.pdf",
    "source/figures/S2_Fig.pdf": "figures/S2_Fig.pdf",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    SUBMISSION.mkdir(parents=True, exist_ok=True)
    for destination, source in {**UPLOAD_FILES, **SOURCE_FILES}.items():
        source_path = REPO_ROOT / source
        if not source_path.is_file() or source_path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty submission input: {source_path}")
        destination_path = SUBMISSION / destination
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination_path)
        print(f"{destination}\t{destination_path.stat().st_size}")
    for required in ("README.md", "Submission_Metadata.md"):
        if not (SUBMISSION / required).is_file():
            raise FileNotFoundError(f"Missing submission guide: {required}")

    checksum_path = SUBMISSION / "SHA256SUMS.txt"
    files = sorted(
        path
        for path in SUBMISSION.rglob("*")
        if path.is_file() and path.name not in {"SHA256SUMS.txt", "validation.json"}
    )
    checksum_path.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(SUBMISSION).as_posix()}\n" for path in files)
    )
    print(f"SHA256SUMS.txt\t{len(files)} files")


if __name__ == "__main__":
    main()
