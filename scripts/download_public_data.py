#!/usr/bin/env python3
"""Download the public animal-tracking sources used by the article.

The downloader uses the official Dryad and USGS ScienceBase APIs and the
authors' public Goto repository. Downloads are resumable and verified against
repository-reported byte counts and SHA-256 digests where available.
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote

import requests
from tqdm import tqdm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DESTINATION = ROOT / "external_inputs" / "raw_sources"
DRYAD_API = "https://datadryad.org/api/v2"
DRYAD_HOST = "https://datadryad.org"
SCIENCEBASE_API = "https://www.sciencebase.gov/catalog"
USER_AGENT = "environmental-history-movement-boundaries/1.0 (public-data reproduction)"

GOTO_REPOSITORY = "https://github.com/YusukeGoto510/Data_Wandering_Albatross.git"
UESAKA_DOI = "10.5061/dryad.tx95x6b2j"
SHEARWATER_DOI = "10.5061/dryad.j9k60"
USGS_PARENT_ID = "5b6cded4e4b0f5d578752dcf"
USGS_REQUIRED = {
    "Deployments.csv",
    "LAAL_eObs.zip",
    "RFBO_iGotU.zip",
    "RFBO_tdr_FastLog.zip",
    "RFBO_tdr_WetDry.zip",
}


def absolute_dryad_url(value: str) -> str:
    if value.startswith("http://") or value.startswith("https://"):
        return value
    return f"{DRYAD_HOST}/{value.lstrip('/')}"


def get_json(session: requests.Session, url: str, attempts: int = 5) -> dict[str, Any]:
    error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                url,
                headers={"Accept": "application/json", "User-Agent": USER_AGENT},
                timeout=(20, 120),
            )
            response.raise_for_status()
            return response.json()
        except Exception as exc:  # network retry boundary
            error = exc
            if attempt < attempts:
                time.sleep(min(30.0, 1.8 ** (attempt - 1)))
    raise RuntimeError(f"metadata request failed after {attempts} attempts: {url}") from error


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def human_bytes(value: int) -> str:
    size = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} TiB"


def download(
    session: requests.Session,
    url: str,
    destination: Path,
    expected_size: int | None,
    expected_sha256: str | None,
    dry_run: bool,
) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file():
        size_ok = expected_size is None or destination.stat().st_size == expected_size
        digest_ok = expected_sha256 is None or sha256(destination).lower() == expected_sha256.lower()
        if size_ok and digest_ok:
            print(f"[SKIP] {destination.name} ({human_bytes(destination.stat().st_size)})")
            return "existing_valid"
    if dry_run:
        print(f"[PLAN] {destination} ({human_bytes(expected_size or 0)}) <- {url}")
        return "planned"

    part = destination.with_suffix(destination.suffix + ".part")
    resume_at = part.stat().st_size if part.exists() else 0
    headers = {"User-Agent": USER_AGENT}
    if resume_at:
        headers["Range"] = f"bytes={resume_at}-"
    response = session.get(url, headers=headers, stream=True, timeout=(20, 300))
    response.raise_for_status()
    if resume_at and response.status_code != 206:
        resume_at = 0
        part.unlink(missing_ok=True)
    mode = "ab" if resume_at else "wb"
    total = expected_size or (resume_at + int(response.headers.get("Content-Length", 0))) or None
    with part.open(mode) as handle, tqdm(
        total=total,
        initial=resume_at,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        desc=destination.name[:48],
    ) as progress:
        for block in response.iter_content(chunk_size=8 << 20):
            if block:
                handle.write(block)
                progress.update(len(block))
    if expected_size is not None and part.stat().st_size != expected_size:
        raise RuntimeError(
            f"size mismatch for {destination.name}: {part.stat().st_size} != {expected_size}"
        )
    if expected_sha256 is not None and sha256(part).lower() != expected_sha256.lower():
        raise RuntimeError(f"SHA-256 mismatch for {destination.name}")
    os.replace(part, destination)
    print(f"[DONE] {destination} ({human_bytes(destination.stat().st_size)})")
    return "downloaded"


def dryad_files(session: requests.Session, doi: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    encoded = quote(f"doi:{doi}", safe="")
    dataset_url = f"{DRYAD_API}/datasets/{encoded}"
    dataset = get_json(session, dataset_url)
    version_link = dataset.get("_links", {}).get("stash:version", {}).get("href")
    if not version_link:
        raise RuntimeError(f"Dryad dataset has no current version link: {doi}")
    version = get_json(session, absolute_dryad_url(version_link))
    next_link = version.get("_links", {}).get("stash:files", {}).get("href")
    if not next_link:
        raise RuntimeError(f"Dryad version has no file listing: {doi}")
    files: list[dict[str, Any]] = []
    visited: set[str] = set()
    while next_link:
        url = absolute_dryad_url(next_link)
        if url in visited:
            raise RuntimeError(f"Dryad pagination loop: {url}")
        visited.add(url)
        page = get_json(session, url)
        files.extend(page.get("_embedded", {}).get("stash:files", []))
        next_link = page.get("_links", {}).get("next", {}).get("href")
    return {"dataset": dataset, "version": version}, files


def save_metadata(destination: Path, name: str, payload: Any) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    (destination / name).write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def download_dryad(
    session: requests.Session,
    doi: str,
    destination: Path,
    patterns: Iterable[str],
    dry_run: bool,
) -> tuple[int, int]:
    metadata, files = dryad_files(session, doi)
    selected = [
        item for item in files
        if any(fnmatch.fnmatch(str(item.get("path", "")), pattern) for pattern in patterns)
    ]
    if not selected:
        raise RuntimeError(f"no Dryad files matched {list(patterns)} for {doi}")
    save_metadata(destination / "metadata", "dryad_dataset_and_version.json", metadata)
    save_metadata(destination / "metadata", "dryad_selected_files.json", selected)
    total = sum(int(item.get("size", 0)) for item in selected)
    free = shutil.disk_usage(destination.parent).free
    print(f"[DRYAD] {doi}: {len(selected)} files, {human_bytes(total)}; free {human_bytes(free)}")
    if not dry_run and free < max(total + (2 << 30), int(total * 1.10)):
        raise RuntimeError(f"insufficient free space for {doi}")
    for index, item in enumerate(selected, start=1):
        path = str(item["path"])
        href = item.get("_links", {}).get("stash:download", {}).get("href")
        if not href:
            raise RuntimeError(f"Dryad file has no download link: {path}")
        print(f"[DRYAD {index}/{len(selected)}] {path}")
        digest = item.get("digest") if str(item.get("digestType", "")).lower() == "sha-256" else None
        download(
            session,
            absolute_dryad_url(href),
            destination / "raw" / path,
            int(item["size"]) if item.get("size") is not None else None,
            str(digest) if digest else None,
            dry_run,
        )
    return len(selected), total


def download_goto(destination: Path, dry_run: bool) -> None:
    if (destination / ".git").is_dir():
        print(f"[SKIP] Goto repository already present: {destination}")
        return
    if destination.exists() and any(destination.iterdir()):
        raise RuntimeError(f"Goto destination exists but is not a git clone: {destination}")
    command = ["git", "clone", "--depth", "1", GOTO_REPOSITORY, str(destination)]
    if dry_run:
        print("[PLAN] " + " ".join(command))
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)


def download_usgs(session: requests.Session, destination: Path, dry_run: bool) -> tuple[int, int]:
    children_url = f"{SCIENCEBASE_API}/items?parentId={USGS_PARENT_ID}&format=json&max=100"
    children = get_json(session, children_url)
    selected: list[dict[str, Any]] = []
    child_payloads: list[dict[str, Any]] = []
    for child in children.get("items", []):
        item = get_json(session, f"{SCIENCEBASE_API}/item/{child['id']}?format=json")
        child_payloads.append(item)
        for record in item.get("files", []):
            if record.get("name") in USGS_REQUIRED:
                selected.append(record)
    found = {str(item.get("name")) for item in selected}
    if found != USGS_REQUIRED:
        raise RuntimeError(f"USGS file inventory changed; missing={sorted(USGS_REQUIRED-found)}")
    save_metadata(destination / "metadata", "sciencebase_children.json", child_payloads)
    total = sum(int(item.get("size", 0)) for item in selected)
    print(f"[USGS] {len(selected)} files, {human_bytes(total)}")
    for index, item in enumerate(sorted(selected, key=lambda row: str(row["name"])), start=1):
        print(f"[USGS {index}/{len(selected)}] {item['name']}")
        download(
            session,
            str(item.get("downloadUri") or item["url"]),
            destination / "raw" / str(item["name"]),
            int(item["size"]) if item.get("size") is not None else None,
            None,
            dry_run,
        )
    return len(selected), total


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        choices=("all", "goto", "uesaka", "usgs", "shearwater"),
        default="all",
        help="Public source to download; 'all' downloads the article's core animal datasets.",
    )
    parser.add_argument(
        "--destination",
        type=Path,
        default=DEFAULT_DESTINATION,
        help="Download root (default: external_inputs/raw_sources).",
    )
    parser.add_argument(
        "--uesaka-profile",
        choices=("gps", "full"),
        default="gps",
        help="GPS downloads 55 *_G.csv files (~2.4 GB); full downloads all sensors (~63.5 GB).",
    )
    parser.add_argument(
        "--accept-large-download",
        action="store_true",
        help="Required with --uesaka-profile full to acknowledge the ~63.5-GB transfer.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Resolve metadata and print the plan without data transfer.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.uesaka_profile == "full" and not args.accept_large_download:
        raise SystemExit("Refusing the ~63.5-GB Uesaka full download without --accept-large-download")
    selected = {"goto", "uesaka", "usgs", "shearwater"} if args.dataset == "all" else {args.dataset}
    args.destination.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    if "goto" in selected:
        download_goto(args.destination / "goto_wandering_albatross", args.dry_run)
    if "uesaka" in selected:
        patterns = ("*",) if args.uesaka_profile == "full" else ("*_G.csv",)
        download_dryad(
            session,
            UESAKA_DOI,
            args.destination / "uesaka_crozet_multisensor",
            patterns,
            args.dry_run,
        )
    if "usgs" in selected:
        download_usgs(session, args.destination / "usgs_hawaiian_seabirds", args.dry_run)
    if "shearwater" in selected:
        download_dryad(
            session,
            SHEARWATER_DOI,
            args.destination / "short_tailed_shearwater",
            ("*",),
            args.dry_run,
        )
    print("PUBLIC_DATA_DOWNLOAD_COMPLETE" if not args.dry_run else "PUBLIC_DATA_DOWNLOAD_PLAN_COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
