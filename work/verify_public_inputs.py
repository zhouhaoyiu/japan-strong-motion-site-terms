#!/usr/bin/env python3
"""Verify the exact public inputs used by the published analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "public_inputs_manifest.tsv"
COEFFICIENT_COLUMNS = {
    "Period",
    "a",
    "b1",
    "b2",
    "b3",
    "c1",
    "c2",
    "c3",
    "d",
    "pd",
    "Dlmin",
    "ps",
    "Vsmax",
    "gNE",
    "gSW",
    "PH",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty input manifest: {path}")
    return rows


def verify_archive(path: Path, members: str) -> None:
    if not members or members == "-":
        return
    required = {name for name in members.split(";") if name}
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
    missing = sorted(required - names)
    if missing:
        raise ValueError(f"{path}: missing archive members: {', '.join(missing)}")


def verify_coefficients(path: Path) -> None:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        header = next(csv.reader(handle))
    columns = {column.strip() for column in header}
    missing = sorted(COEFFICIENT_COLUMNS - columns)
    if missing:
        raise ValueError(f"{path}: missing coefficient columns: {', '.join(missing)}")


def verify(manifest: Path, skip_hash: bool) -> None:
    for row in read_manifest(manifest):
        path = ROOT / row["path"]
        if not path.is_file():
            raise FileNotFoundError(f"missing {row['role']}: {path}")
        expected_size = int(row["bytes"])
        if path.stat().st_size != expected_size:
            raise ValueError(f"{path}: expected {expected_size} bytes, found {path.stat().st_size}")
        if not skip_hash:
            observed_hash = sha256(path)
            if observed_hash != row["sha256"]:
                raise ValueError(f"{path}: SHA-256 mismatch")
        verify_archive(path, row["required_members"])
        if row["role"] == "mf2013_coefficients":
            verify_coefficients(path)
        print(f"verified {row['role']}: {row['path']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--skip-hash", action="store_true", help="check size and schema without hashing large files")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    verify(args.manifest, args.skip_hash)
