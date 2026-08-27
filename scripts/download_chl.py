#!/usr/bin/env python3
"""Download frozen Copernicus Marine CHL subsets from public request manifests."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_ROOT = ROOT / "config" / "chl_requests"
DEFAULT_OUTPUT = ROOT / "external_inputs" / "environment"
DEFAULT_STATUS = ROOT / "external_inputs" / "download_status" / "chl"
REQUIRED_FIELDS = {
    "dataset_id", "variables", "request_lat_min", "request_lat_max",
    "request_lon_min", "request_lon_max", "start_date", "end_date", "output_filename",
}


def elapsed_text(seconds: float) -> str:
    hours, remainder = divmod(max(0, int(seconds)), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:d}h{minutes:02d}m{seconds:02d}s" if hours else f"{minutes:02d}m{seconds:02d}s"


def manifests_from_args(args: argparse.Namespace) -> list[Path]:
    if args.all:
        return sorted(MANIFEST_ROOT.glob("*.csv"))
    if not args.manifest:
        raise SystemExit("Use --all or provide one or more --manifest paths")
    return [path.resolve() for path in args.manifest]


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = REQUIRED_FIELDS - fields
        if missing:
            raise RuntimeError(f"{path} lacks required columns: {sorted(missing)}")
        return list(reader)


def command_for(row: dict[str, str], output_directory: Path) -> list[str]:
    command = ["copernicusmarine", "subset", "--dataset-id", row["dataset_id"]]
    for variable in row["variables"].split(","):
        if variable.strip():
            command.extend(["--variable", variable.strip()])
    command.extend([
        "--minimum-longitude", row["request_lon_min"],
        "--maximum-longitude", row["request_lon_max"],
        "--minimum-latitude", row["request_lat_min"],
        "--maximum-latitude", row["request_lat_max"],
        "--start-datetime", f"{row['start_date']}T00:00:00",
        "--end-datetime", f"{row['end_date']}T23:59:59",
        "--output-directory", str(output_directory),
        "--output-filename", row["output_filename"],
    ])
    return command


def write_status(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    fields = sorted({key for row in rows for key in row}) if rows else ["status"]
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="Run every manifest in config/chl_requests.")
    group.add_argument("--manifest", type=Path, action="append", help="Manifest path; repeat for several files.")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--status-root", type=Path, default=DEFAULT_STATUS)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--start-row", type=int, default=0)
    parser.add_argument("--max-downloads", type=int)
    parser.add_argument("--stop-on-error", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.dry_run and shutil.which("copernicusmarine") is None:
        raise SystemExit("copernicusmarine is not installed; run 'make setup' first")
    manifests = manifests_from_args(args)
    indexed: list[tuple[Path, int, dict[str, str]]] = []
    for manifest in manifests:
        for row_index, row in enumerate(read_manifest(manifest)):
            if row_index >= args.start_row:
                indexed.append((manifest, row_index, row))
    if args.max_downloads is not None:
        indexed = indexed[: args.max_downloads]
    if not indexed:
        raise SystemExit("No CHL requests selected")

    total = len(indexed)
    started = time.monotonic()
    statuses: dict[str, list[dict[str, Any]]] = {}
    errors = 0
    for sequence, (manifest, row_index, row) in enumerate(indexed, start=1):
        collection = manifest.stem.removesuffix("_requests")
        output_directory = args.output_root / collection
        output_directory.mkdir(parents=True, exist_ok=True)
        output = output_directory / row["output_filename"]
        status = "planned"
        command = command_for(row, output_directory)
        if output.is_file() and output.stat().st_size > 0 and not args.force:
            status = "existing"
        elif args.dry_run:
            status = "dry_run"
            print("[COMMAND] " + " ".join(command))
        else:
            if args.force and output.exists():
                output.unlink()
            process = subprocess.run(command, text=True, capture_output=True)
            status = "downloaded" if process.returncode == 0 and output.is_file() and output.stat().st_size > 0 else "error"
            if status == "error":
                errors += 1
        elapsed = time.monotonic() - started
        rate = sequence / elapsed if elapsed else 0.0
        eta = (total - sequence) / rate if rate else 0.0
        print(
            f"[CHL {sequence}/{total}] {100*sequence/total:.1f}% status={status} "
            f"file={output.name} elapsed={elapsed_text(elapsed)} eta={elapsed_text(eta)}",
            flush=True,
        )
        record = {
            "manifest": manifest.name,
            "manifest_row": row_index,
            "output_path": str(output),
            "status": status,
            "command": " ".join(command),
        }
        if not args.dry_run and status == "error":
            record["returncode"] = process.returncode
            record["stdout_tail"] = process.stdout[-2000:]
            record["stderr_tail"] = process.stderr[-4000:]
        statuses.setdefault(manifest.stem, []).append(record)
        write_status(args.status_root / f"{manifest.stem}_status.csv", statuses[manifest.stem])
        if status == "error" and args.stop_on_error:
            raise RuntimeError(f"Copernicus request failed: {output.name}")

    args.status_root.mkdir(parents=True, exist_ok=True)
    (args.status_root / "summary.json").write_text(
        json.dumps({
            "manifests": [str(path) for path in manifests],
            "requests": total,
            "errors": errors,
            "dry_run": args.dry_run,
            "elapsed_seconds": time.monotonic() - started,
        }, indent=2),
        encoding="utf-8",
    )
    if errors:
        raise SystemExit(f"CHL download completed with {errors} failed request(s)")
    print("CHL_DOWNLOAD_COMPLETE" if not args.dry_run else "CHL_DOWNLOAD_PLAN_COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
