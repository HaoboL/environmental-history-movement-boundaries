#!/usr/bin/env python3
"""Regenerate the repository file manifest after intentional release changes."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata" / "FILE_MANIFEST.csv"
EXCLUDED = {".git", ".venv", "__pycache__", ".pytest_cache", ".mypy_cache", "output", "external_inputs"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    files = [
        path for path in ROOT.rglob("*")
        if path.is_file()
        and path != MANIFEST
        and not any(part in EXCLUDED for part in path.relative_to(ROOT).parts)
    ]
    temporary = MANIFEST.with_suffix(".csv.part")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["path", "bytes", "sha256"],
            lineterminator="\n",
        )
        writer.writeheader()
        for path in sorted(files):
            writer.writerow({
                "path": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            })
    temporary.replace(MANIFEST)
    print(f"MANIFEST_UPDATED files={len(files)} path={MANIFEST}")


if __name__ == "__main__":
    main()
