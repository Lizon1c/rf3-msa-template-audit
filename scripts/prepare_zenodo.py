#!/usr/bin/env python3
import argparse
import hashlib
import json
import zipfile
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_archive(output, name, files, source_root, arc_prefix):
    path = output / name
    if path.exists():
        raise FileExistsError(f"{path} already exists; choose a new output directory")
    total = 0
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=True) as archive:
        for source in files:
            relative = source.relative_to(source_root)
            archive.write(source, f"{arc_prefix}/{relative}")
            total += source.stat().st_size
    return {
        "archive": name,
        "source_root": str(source_root),
        "files": len(files),
        "uncompressed_bytes": total,
        "archive_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def main():
    parser = argparse.ArgumentParser(description="Prepare analysis-ready Zenodo archives without uploading them")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--feature-root", type=Path, required=True, help="Root containing e2/ and e427_blocksplit_lightprobe/")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    features = args.feature_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    jobs = [
        (
            "e428_sample_structures.zip",
            sorted(p for p in (repo / "data/structures/e428").rglob("*") if p.is_file()),
            repo / "data/structures/e428",
            "e428",
        ),
        (
            "e427_feature_cache.zip",
            sorted((features / "e427_blocksplit_lightprobe/cache").glob("*.npy")),
            features / "e427_blocksplit_lightprobe/cache",
            "e427_feature_cache",
        ),
        (
            "e2_distogram_outputs.zip",
            sorted((features / "e2").glob("*/*/disto.pt")),
            features / "e2",
            "e2",
        ),
        (
            "canonical_results.zip",
            sorted(p for p in (repo / "data/results").rglob("*.json") if p.is_file()),
            repo / "data/results",
            "results",
        ),
    ]

    manifest = []
    for name, files, root, prefix in jobs:
        if not files:
            raise FileNotFoundError(f"No input files for {name} under {root}")
        record = add_archive(output, name, files, root, prefix)
        manifest.append(record)
        print(f"{name}: {record['files']} files, {record['archive_bytes'] / 1e6:.1f} MB")

    (output / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    checksums = "\n".join(f"{row['sha256']}  {row['archive']}" for row in manifest) + "\n"
    (output / "SHA256SUMS.txt").write_text(checksums)
    print(f"Total: {sum(row['archive_bytes'] for row in manifest) / 1e9:.3f} GB in {len(manifest)} archives")


if __name__ == "__main__":
    main()
