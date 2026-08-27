#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for the public release."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metadata/FILE_MANIFEST.csv"
EXCLUDED_PARTS = {".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".venv", "output"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def included(path: Path) -> bool:
    relative = path.relative_to(ROOT)
    return path != OUTPUT and not any(part in EXCLUDED_PARTS for part in relative.parts)


def main() -> int:
    paths = sorted(path for path in ROOT.rglob("*") if path.is_file() and included(path))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["path", "bytes", "sha256"])
        for path in paths:
            writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, digest(path)])
    print(f"wrote {OUTPUT.relative_to(ROOT)} with {len(paths)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
