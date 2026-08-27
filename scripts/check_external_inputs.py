#!/usr/bin/env python3
"""Check downloaded public sources and Copernicus CHL subsets."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "external_inputs" / "raw_sources"
ENVIRONMENT = ROOT / "external_inputs" / "environment"
MANIFESTS = ROOT / "config" / "chl_requests"


def pass_or_fail(label: str, ok: bool, detail: str) -> bool:
    print(f"[{'PASS' if ok else 'MISSING'}] {label}: {detail}")
    return ok


def main() -> int:
    checks: list[bool] = []
    goto = RAW / "goto_wandering_albatross"
    checks.append(pass_or_fail("Goto tracks", (goto / ".git").is_dir(), str(goto)))

    uesaka = list((RAW / "uesaka_crozet_multisensor" / "raw").glob("*_G.csv"))
    checks.append(pass_or_fail("Uesaka GPS", len(uesaka) > 0, f"{len(uesaka)} *_G.csv files"))

    shearwater = list((RAW / "short_tailed_shearwater" / "raw").glob("*"))
    checks.append(pass_or_fail("Short-tailed shearwater", any(path.is_file() for path in shearwater), f"{len(shearwater)} entries"))

    usgs_names = {
        "Deployments.csv", "LAAL_eObs.zip", "RFBO_iGotU.zip",
        "RFBO_tdr_FastLog.zip", "RFBO_tdr_WetDry.zip",
    }
    usgs_dir = RAW / "usgs_hawaiian_seabirds" / "raw"
    found_usgs = {path.name for path in usgs_dir.glob("*") if path.is_file()}
    checks.append(pass_or_fail("USGS seabirds", usgs_names <= found_usgs, f"{len(usgs_names & found_usgs)}/{len(usgs_names)} required files"))

    for manifest in sorted(MANIFESTS.glob("*.csv")):
        collection = manifest.stem.removesuffix("_requests")
        with manifest.open(newline="", encoding="utf-8-sig") as handle:
            expected = [row["output_filename"] for row in csv.DictReader(handle)]
        folder = ENVIRONMENT / collection
        present = sum((folder / name).is_file() and (folder / name).stat().st_size > 0 for name in expected)
        checks.append(pass_or_fail(f"CHL {collection}", present == len(expected), f"{present}/{len(expected)} subsets"))

    if all(checks):
        print("EXTERNAL_INPUT_CHECK_PASS")
        return 0
    print("EXTERNAL_INPUT_CHECK_INCOMPLETE")
    return 1


if __name__ == "__main__":
    sys.exit(main())
