#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "metadata/landing_timing_specification_cn.md"
INPUTS = ROOT / "external_inputs/landing_timing"
U10 = INPUTS / "uesaka_10_to_100m_events_with_landing.csv.gz"
UALL_BLIND = INPUTS / "uesaka_cross_scale_events.csv.gz"
UALL_OUTCOME = INPUTS / "uesaka_cross_scale_events_with_landing.csv.gz"
ULAND = INPUTS / "uesaka_gps_state_transitions.csv.gz"
UPHASE = INPUTS / "uesaka_phase_runs.csv.gz"
D09_LAND = INPUTS / "stomach_temperature_landing_sample.parquet"
D09_ING = INPUTS / "stomach_temperature_ingestion_events.parquet"
D02_EP = INPUTS / "behavior_episodes.parquet"
D02_BOUND = INPUTS / "behavior_episode_boundaries.parquet"

EXPECTED_HASHES = {
    U10: "86fcaa77f31e2ebc1489cc78e011a42eca7b86704bd2a7725140a1c26cbb02c9",
    UALL_BLIND: "848c479931cb198b03998d02681d4382c576e3853bc91187a1665b2946108ac1",
    UALL_OUTCOME: "6b6b66fe931c11f093233fb522afeb47140b3d8c13ce3038cdb090285a852ee0",
    ULAND: "1e3efb631cc4dd3393217f246e3ad85ebfa183b1cf819700950e8ef43acf7cdb",
    UPHASE: "46fce2ab1791b533f33b4d6e1d2787f023f60da25471a15f05ca07a287e1b90c",
    D09_LAND: "1167d104e43e28dd8905f1cd73754d67778e0f17720d295d1608a24747897340",
    D09_ING: "9560190b566fcbb80513393f88b50bc07b1c4c1f92e6ed4f98b3f70e708b8a14",
    D02_EP: "ddc0ddf160bbbc12ae9220a67241f0c14890406a4024b2d71a04ab94d20e59f9",
    D02_BOUND: "0c114b2849689f7381a85ed62a87034fa332436e896cfde4b0a8105355de8e73",
}

SEED = 2026082245
U_SCALES = [10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def affinity() -> list[int]:
    try:
        return sorted(os.sched_getaffinity(0))
    except Exception:
        return []


def write_progress(out: Path, stage: str, completed: int, total: int, started: float, note: str = "") -> None:
    elapsed = time.monotonic() - started
    rec = {
        "stage": stage,
        "completed": int(completed),
        "total": int(total),
        "percent": round(100.0 * completed / total, 2) if total else None,
        "elapsed_s": round(elapsed, 3),
        "throughput_per_s": round(completed / elapsed, 5) if elapsed > 0 else None,
        "eta_s": round((total - completed) * elapsed / completed, 3) if completed > 0 else None,
        "note": note,
    }
    with (out / "progress.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    print(json.dumps(rec, ensure_ascii=False), flush=True)


def truthy(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def to_epoch_s(s: pd.Series) -> np.ndarray:
    dt = pd.to_datetime(s, utc=True, errors="coerce")
    # Parquet columns may retain microsecond resolution while CSV-derived columns
    # are nanosecond-resolution. Normalize explicitly before integer conversion.
    ns = dt.astype("datetime64[ns, UTC]").astype("int64").to_numpy(dtype=np.float64)
    ns[pd.isna(dt).to_numpy()] = np.nan
    return ns / 1e9


def holm(pvals: list[float]) -> list[float]:
    a = np.asarray(pvals, dtype=float)
    out = np.full(len(a), np.nan)
    finite = np.isfinite(a)
    idx = np.where(finite)[0]
    if not len(idx):
        return out.tolist()
    order = idx[np.argsort(a[idx])]
    m = len(order)
    running = 0.0
    for rank, j in enumerate(order):
        running = max(running, (m - rank) * a[j])
        out[j] = min(1.0, running)
    return out.tolist()


def bootstrap_ci(values: np.ndarray, reps: int, seed: int) -> tuple[float, float, int]:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan, np.nan, 0
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=float)
    chunk = 1000
    for i in range(0, reps, chunk):
        n = min(chunk, reps - i)
        ix = rng.integers(0, len(v), size=(n, len(v)))
        draws[i : i + n] = np.mean(v[ix], axis=1)
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975)), int(len(v))


def nearest_flags(anchor_t: np.ndarray, event_t: np.ndarray, window_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    anchor_t = np.asarray(anchor_t, dtype=float)
    event_t = np.sort(np.asarray(event_t, dtype=float))
    n = len(anchor_t)
    pre_d = np.full(n, np.nan)
    post_d = np.full(n, np.nan)
    if not len(event_t):
        return np.zeros(n, dtype=bool), np.zeros(n, dtype=bool), pre_d, post_d
    left = np.searchsorted(event_t, anchor_t, side="left")
    right = np.searchsorted(event_t, anchor_t, side="right")
    ok_pre = left > 0
    pre_d[ok_pre] = anchor_t[ok_pre] - event_t[left[ok_pre] - 1]
    ok_post = right < len(event_t)
    post_d[ok_post] = event_t[right[ok_post]] - anchor_t[ok_post]
    return (pre_d > 0) & (pre_d <= window_s), (post_d > 0) & (post_d <= window_s), pre_d, post_d


def attach_actual_flags(anchors: pd.DataFrame, events_by_unit: dict[str, np.ndarray], window_s: float) -> pd.DataFrame:
    ans = anchors.copy().reset_index(drop=True)
    ans["pre"] = False
    ans["post"] = False
    ans["nearest_pre_s"] = np.nan
    ans["nearest_post_s"] = np.nan
    for unit, ix in ans.groupby("unit", sort=False).groups.items():
        pos = np.asarray(list(ix), dtype=int)
        pre, post, pre_d, post_d = nearest_flags(ans.loc[pos, "t"].to_numpy(float), events_by_unit.get(str(unit), np.array([])), window_s)
        ans.loc[pos, "pre"] = pre
        ans.loc[pos, "post"] = post
        ans.loc[pos, "nearest_pre_s"] = pre_d
        ans.loc[pos, "nearest_post_s"] = post_d
    ans["pre"] = ans["pre"].astype(bool)
    ans["post"] = ans["post"].astype(bool)
    return ans


def per_bird_values(flagged: pd.DataFrame) -> pd.DataFrame:
    x = flagged.groupby("bird", sort=True).agg(pre=("pre", "mean"), post=("post", "mean"), events=("pre", "size")).reset_index()
    x["diff"] = x["post"] - x["pre"]
    return x


def phase_test(
    flagged: pd.DataFrame,
    events_by_unit: dict[str, np.ndarray],
    spans: dict[str, tuple[float, float]],
    window_s: float,
    reps: int,
    seed: int,
    progress_cb=None,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    base = flagged.reset_index(drop=True).copy()
    birds = sorted(base["bird"].astype(str).unique())
    bird_map = {b: i for i, b in enumerate(birds)}
    bird_code = base["bird"].astype(str).map(bird_map).to_numpy(int)
    bird_n = np.bincount(bird_code, minlength=len(birds)).astype(float)
    groups = {str(u): np.asarray(list(ix), dtype=int) for u, ix in base.groupby("unit", sort=False).groups.items()}
    actual_bird = per_bird_values(base)
    actual_pre = float(actual_bird["pre"].mean())
    actual_post = float(actual_bird["post"].mean())
    actual_diff = actual_post - actual_pre
    rng = np.random.default_rng(seed)
    phase_pre = np.empty(reps)
    phase_post = np.empty(reps)
    phase_diff = np.empty(reps)
    bird_pre_acc = np.zeros(len(birds))
    bird_post_acc = np.zeros(len(birds))
    for r in range(reps):
        pre_all = np.zeros(len(base), dtype=float)
        post_all = np.zeros(len(base), dtype=float)
        for unit, pos in groups.items():
            ev = events_by_unit.get(unit, np.array([], dtype=float))
            if not len(ev):
                continue
            lo, hi = spans[unit]
            length = hi - lo
            if not np.isfinite(length) or length <= 2 * window_s:
                continue
            offset = rng.uniform(0.0, length)
            shifted = np.sort(((ev - lo + offset) % length) + lo)
            pre, post, _, _ = nearest_flags(base.loc[pos, "t"].to_numpy(float), shifted, window_s)
            pre_all[pos] = pre
            post_all[pos] = post
        bpre = np.bincount(bird_code, weights=pre_all, minlength=len(birds)) / bird_n
        bpost = np.bincount(bird_code, weights=post_all, minlength=len(birds)) / bird_n
        phase_pre[r] = float(np.mean(bpre))
        phase_post[r] = float(np.mean(bpost))
        phase_diff[r] = phase_post[r] - phase_pre[r]
        bird_pre_acc += bpre
        bird_post_acc += bpost
        if progress_cb is not None and ((r + 1) % max(1, reps // 20) == 0 or r + 1 == reps):
            progress_cb(r + 1, reps)
    exp_pre_bird = bird_pre_acc / reps
    exp_post_bird = bird_post_acc / reps
    actual_by_bird = actual_bird.set_index("bird").loc[birds]
    residual = pd.DataFrame({
        "bird": birds,
        "actual_pre": actual_by_bird["pre"].to_numpy(float),
        "actual_post": actual_by_bird["post"].to_numpy(float),
        "null_pre": exp_pre_bird,
        "null_post": exp_post_bird,
    })
    residual["resid_pre"] = residual["actual_pre"] - residual["null_pre"]
    residual["resid_post"] = residual["actual_post"] - residual["null_post"]
    residual["resid_diff"] = (residual["actual_post"] - residual["actual_pre"]) - (residual["null_post"] - residual["null_pre"])
    phase = pd.DataFrame({"rep": np.arange(reps), "pre": phase_pre, "post": phase_post, "diff": phase_diff})
    summary = {
        "actual_pre": actual_pre,
        "actual_post": actual_post,
        "actual_diff": actual_diff,
        "null_pre_mean": float(np.mean(phase_pre)),
        "null_post_mean": float(np.mean(phase_post)),
        "null_diff_mean": float(np.mean(phase_diff)),
        "resid_pre": float(residual["resid_pre"].mean()),
        "resid_post": float(residual["resid_post"].mean()),
        "resid_diff": float(residual["resid_diff"].mean()),
        "p_pre_upper": float((1 + np.sum(phase_pre >= actual_pre)) / (reps + 1)),
        "p_post_upper": float((1 + np.sum(phase_post >= actual_post)) / (reps + 1)),
        "p_diff_upper": float((1 + np.sum(phase_diff >= actual_diff)) / (reps + 1)),
    }
    return summary, phase, residual


def bin_summary(flagged: pd.DataFrame, events_by_unit: dict[str, np.ndarray], bins: list[tuple[float, float]], reps: int, seed: int) -> pd.DataFrame:
    rows = []
    for bi, (lo, hi) in enumerate(bins):
        x = flagged[["unit", "bird", "t"]].copy().reset_index(drop=True)
        pre_all = np.zeros(len(x), dtype=bool)
        post_all = np.zeros(len(x), dtype=bool)
        for unit, ix in x.groupby("unit", sort=False).groups.items():
            pos = np.asarray(list(ix), dtype=int)
            t = x.loc[pos, "t"].to_numpy(float)
            ev = np.sort(events_by_unit.get(str(unit), np.array([], dtype=float)))
            if not len(ev):
                continue
            post_left = np.searchsorted(ev, t + lo, side="right")
            post_right = np.searchsorted(ev, t + hi, side="right")
            pre_left = np.searchsorted(ev, t - hi, side="left")
            pre_right = np.searchsorted(ev, t - lo, side="left")
            post_all[pos] = post_right > post_left
            pre_all[pos] = pre_right > pre_left
        x["pre"] = pre_all
        x["post"] = post_all
        b = per_bird_values(x)
        ci_lo, ci_hi, n_birds = bootstrap_ci(b["diff"].to_numpy(float), reps, seed + bi)
        rows.append({
            "bin_low_s": lo,
            "bin_high_s": hi,
            "events": len(x),
            "birds": n_birds,
            "pre_bird_equal": float(b["pre"].mean()),
            "post_bird_equal": float(b["post"].mean()),
            "post_minus_pre": float(b["diff"].mean()),
            "ci_low": ci_lo,
            "ci_high": ci_hi,
        })
    return pd.DataFrame(rows)


def load_uesaka() -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, tuple[float, float]], pd.DataFrame]:
    cols = ["event_id", "delta_m", "segment_id", "bird_id", "deployment_id", "endpoint_segment_index", "endpoint_time_utc", "phase", "primary_risk_eligible_gps"]
    a = pd.read_csv(U10, usecols=cols)
    a = a[(a["delta_m"] <= 100) & (a["phase"] == "middle")]
    b = pd.read_csv(UALL_OUTCOME, usecols=cols)
    b = b[(b["delta_m"] > 100) & (b["phase"] == "middle")]
    x = pd.concat([a, b], ignore_index=True)
    x = x[truthy(x["primary_risk_eligible_gps"])].copy()
    x["delta_m"] = x["delta_m"].astype(int)
    x["endpoint_time_utc"] = pd.to_datetime(x["endpoint_time_utc"], utc=True)
    x = x.drop_duplicates(["delta_m", "event_id"])
    blind = pd.read_csv(UALL_BLIND, usecols=["event_id", "delta_m", "endpoint_time_utc", "phase"])
    blind = blind[(blind["phase"] == "middle") & (blind["delta_m"] > 100)].copy()
    blind["delta_m"] = blind["delta_m"].astype(int)
    key_out = set(zip(b["delta_m"].astype(int), b["event_id"].astype(str)))
    key_blind = set(zip(blind["delta_m"], blind["event_id"].astype(str)))
    if not key_out.issubset(key_blind):
        raise RuntimeError("Uesaka outcome keys are not a subset of blind endpoint keys")

    runs = pd.read_csv(UPHASE)
    runs = runs[runs["phase"] == "middle"].copy()
    merged = x.merge(runs, on=["segment_id", "bird_id", "deployment_id", "phase"], how="left", validate="many_to_many")
    inside = (merged["endpoint_segment_index"] >= merged["start_index"]) & (merged["endpoint_segment_index"] < merged["end_index_exclusive"])
    merged = merged[inside].copy()
    counts = merged.groupby(["delta_m", "event_id"]).size()
    if (counts != 1).any():
        raise RuntimeError("Each middle endpoint must map to exactly one middle phase run")
    merged["bilateral_middle_10m"] = (
        (merged["endpoint_segment_index"] - merged["start_index"] >= 660)
        & (merged["end_index_exclusive"] - merged["endpoint_segment_index"] > 660)
    )
    merged = merged[merged["bilateral_middle_10m"]].copy()
    merged["run_id"] = (
        merged["deployment_id"].astype(str) + "|" + merged["segment_id"].astype(str) + "|"
        + merged["start_index"].astype(str) + "|" + merged["end_index_exclusive"].astype(str)
    )
    et = to_epoch_s(merged["endpoint_time_utc"])
    merged["pred_run_start_s"] = et - (merged["endpoint_segment_index"] - merged["start_index"]).to_numpy(float)
    merged["pred_run_end_s"] = et + (merged["end_index_exclusive"] - merged["endpoint_segment_index"]).to_numpy(float)
    run_diag = merged.groupby("run_id").agg(
        start_min=("pred_run_start_s", "min"), start_max=("pred_run_start_s", "max"),
        end_min=("pred_run_end_s", "min"), end_max=("pred_run_end_s", "max"),
        deployment_id=("deployment_id", "first"), bird_id=("bird_id", "first"),
    ).reset_index()
    run_diag["start_spread_s"] = run_diag["start_max"] - run_diag["start_min"]
    run_diag["end_spread_s"] = run_diag["end_max"] - run_diag["end_min"]
    spans = {r.run_id: (float((r.start_min + r.start_max) / 2), float((r.end_min + r.end_max) / 2)) for r in run_diag.itertuples()}

    land = pd.read_csv(ULAND, usecols=["deployment_id", "bird_id", "event_type", "event_time_utc"])
    land = land[land["event_type"] == "landing_gps_conservative"].copy()
    land["t"] = to_epoch_s(land["event_time_utc"])
    events_by_unit: dict[str, np.ndarray] = {}
    for r in run_diag.itertuples():
        lo, hi = spans[r.run_id]
        vals = land.loc[(land["deployment_id"] == r.deployment_id) & (land["t"] >= lo) & (land["t"] <= hi), "t"].to_numpy(float)
        events_by_unit[r.run_id] = np.sort(vals)
    anchors = merged[["delta_m", "event_id", "bird_id", "run_id", "endpoint_time_utc"]].copy()
    anchors = anchors.rename(columns={"bird_id": "bird", "run_id": "unit"})
    anchors["t"] = to_epoch_s(anchors["endpoint_time_utc"])
    return anchors, events_by_unit, spans, run_diag


def load_d09() -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, tuple[float, float]], pd.DataFrame, pd.DataFrame]:
    land = pd.read_parquet(D09_LAND)
    ing = pd.read_parquet(D09_ING)
    land = land[truthy(land["primary_landing_qc_eligible"])].copy()
    for c in ["landing_time_dt", "stomach_recording_start", "stomach_recording_end", "bout_start", "bout_end"]:
        land[c] = pd.to_datetime(land[c], utc=True)
    ing["event_time_dt"] = pd.to_datetime(ing["event_time_dt"], utc=True)
    land = land[(land["landing_time_dt"] - land["stomach_recording_start"] >= pd.Timedelta(hours=1)) & (land["stomach_recording_end"] - land["landing_time_dt"] >= pd.Timedelta(hours=1))].copy()
    land["unit"] = land["deployment_id"].astype(str)
    land["bird"] = land["bird_id"].astype(str)
    land["t"] = to_epoch_s(land["landing_time_dt"])
    events_by_unit = {str(k): np.sort(to_epoch_s(g["event_time_dt"])) for k, g in ing.groupby("deployment_id")}
    spans = {}
    for dep, g in land.groupby("deployment_id"):
        spans[str(dep)] = (float(g["stomach_recording_start"].min().value / 1e9), float(g["stomach_recording_end"].max().value / 1e9))
    anchors = land[["landing_id", "deployment_id", "bird", "unit", "t", "landing_time_dt", "bout_start", "bout_end"]].copy()
    bout_rows = []
    for dep, lg in land.groupby("deployment_id"):
        ev = ing.loc[ing["deployment_id"] == dep, ["source_event_number", "event_time_dt"]]
        for lr in lg.itertuples():
            hit = ev[(ev["event_time_dt"] >= lr.bout_start) & (ev["event_time_dt"] <= lr.bout_end)]
            for er in hit.itertuples():
                bout_rows.append({
                    "landing_id": lr.landing_id,
                    "deployment_id": dep,
                    "bird_id": lr.bird_id,
                    "source_event_number": er.source_event_number,
                    "landing_time_utc": lr.landing_time_dt,
                    "ingestion_time_utc": er.event_time_dt,
                    "landing_to_ingestion_s": (er.event_time_dt - lr.landing_time_dt).total_seconds(),
                    "bout_duration_s": (lr.bout_end - lr.bout_start).total_seconds(),
                })
    return anchors, events_by_unit, spans, pd.DataFrame(bout_rows), ing


def load_d02(method: str) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, tuple[float, float]]]:
    ep = pd.read_parquet(D02_EP)
    ep = ep[ep["event_type"] == "ingestion_candidate_episode"].copy()
    ep["event_time_utc"] = pd.to_datetime(ep["event_time_utc"], utc=True)
    b = pd.read_parquet(D02_BOUND, columns=["bird_id", "trip_id", "method", "rd_phase", "boundary_time_utc"])
    b = b[(b["method"] == method) & (b["rd_phase"] == "endpoint")].copy()
    b["boundary_time_utc"] = pd.to_datetime(b["boundary_time_utc"], utc=True)
    b["unit"] = b["bird_id"].astype(str) + "|" + b["trip_id"].astype(str)
    ep["unit"] = ep["bird_id"].astype(str) + "|" + ep["trip_id"].astype(str)
    spans = {str(k): (float(g["boundary_time_utc"].min().value / 1e9), float(g["boundary_time_utc"].max().value / 1e9)) for k, g in b.groupby("unit")}
    events_by_unit = {str(k): np.sort(to_epoch_s(g["boundary_time_utc"])) for k, g in b.groupby("unit")}
    ep["t"] = to_epoch_s(ep["event_time_utc"])
    ep = ep[ep["unit"].isin(spans)].copy()
    ep = ep[ep.apply(lambda r: (r["t"] - spans[str(r["unit"])][0] >= 900) and (spans[str(r["unit"])][1] - r["t"] >= 900), axis=1)].copy()
    anchors = ep[["event_id", "bird_id", "unit", "t", "event_time_utc"]].rename(columns={"bird_id": "bird"})
    # phase_test assumes anchors are reference and events are response. Here anchors are ingestion;
    # swap labels afterward: its pre endpoint means ingestion occurs after an endpoint.
    return anchors, events_by_unit, spans


def summarize_test(flagged: pd.DataFrame, phase_summary: dict, residual: pd.DataFrame, boot_reps: int, seed: int) -> dict:
    b = per_bird_values(flagged)
    raw_lo, raw_hi, n_birds = bootstrap_ci(b["diff"].to_numpy(float), boot_reps, seed)
    rp_lo, rp_hi, _ = bootstrap_ci(residual["resid_pre"].to_numpy(float), boot_reps, seed + 1)
    ro_lo, ro_hi, _ = bootstrap_ci(residual["resid_post"].to_numpy(float), boot_reps, seed + 2)
    rd_lo, rd_hi, _ = bootstrap_ci(residual["resid_diff"].to_numpy(float), boot_reps, seed + 3)
    out = dict(phase_summary)
    out.update({
        "events": int(len(flagged)),
        "birds": int(n_birds),
        "pre_events": int(flagged["pre"].sum()),
        "post_events": int(flagged["post"].sum()),
        "raw_diff_ci_low": raw_lo,
        "raw_diff_ci_high": raw_hi,
        "resid_pre_ci_low": rp_lo,
        "resid_pre_ci_high": rp_hi,
        "resid_post_ci_low": ro_lo,
        "resid_post_ci_high": ro_hi,
        "resid_diff_ci_low": rd_lo,
        "resid_diff_ci_high": rd_hi,
    })
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--bootstrap-reps", type=int, default=20000)
    ap.add_argument("--phase-reps", type=int, default=5000)
    ap.add_argument("--validation-only", action="store_true")
    args = ap.parse_args()
    out = (ROOT / args.output).resolve() if not Path(args.output).is_absolute() else Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty output: {out}")
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    aff = affinity()
    if 1 in aff:
        raise RuntimeError(f"CPU1 forbidden, affinity={aff}")
    if aff and aff != [7]:
        raise RuntimeError(f"Formal/validation process must be pinned to CPU7, affinity={aff}")
    for p, expected in EXPECTED_HASHES.items():
        actual = sha256(p)
        if actual != expected:
            raise RuntimeError(f"Input hash mismatch: {p} {actual} != {expected}")
    write_progress(out, "start", 0, 1, started, f"validation_only={args.validation_only} affinity={aff}")

    u_anchor, u_events, u_spans, u_run_diag = load_uesaka()
    d09_anchor, d09_events, d09_spans, d09_bout, d09_ing = load_d09()
    d02_loaded = {m: load_d02(m) for m in ["RD100", "RD200"]}
    coverage_rows = []
    for scale, g in u_anchor.groupby("delta_m"):
        coverage_rows.append({"panel": "uesaka", "setting": str(int(scale)), "events": len(g), "birds": g["bird"].nunique(), "units": g["unit"].nunique()})
    coverage_rows.append({"panel": "d09", "setting": "landing_ingestion", "events": len(d09_anchor), "birds": d09_anchor["bird"].nunique(), "units": d09_anchor["unit"].nunique()})
    for method, (a, _, _) in d02_loaded.items():
        coverage_rows.append({"panel": "d02", "setting": method, "events": len(a), "birds": a["bird"].nunique(), "units": a["unit"].nunique()})
    pd.DataFrame(coverage_rows).to_csv(out / "coverage_inventory.csv", index=False)
    u_run_diag.to_csv(out / "uesaka_run_time_reconstruction_audit.csv", index=False)
    synth_anchor = np.array([100.0, 200.0])
    synth_events = np.array([110.0, 210.0])
    sp, so, _, _ = nearest_flags(synth_anchor, synth_events, 20.0)
    synthetic_pass = bool((~sp).all() and so.all())
    validation = {
        "status": "PASS" if synthetic_pass else "FAIL",
        "synthetic_post_detection": synthetic_pass,
        "uesaka_scales": sorted(u_anchor["delta_m"].unique().astype(int).tolist()),
        "uesaka_events": int(len(u_anchor)),
        "uesaka_birds": int(u_anchor["bird"].nunique()),
        "uesaka_run_start_spread_max_s": float(u_run_diag["start_spread_s"].max()),
        "uesaka_run_end_spread_max_s": float(u_run_diag["end_spread_s"].max()),
        "d09_landings_bilateral_60m": int(len(d09_anchor)),
        "d09_ingestions": int(len(d09_ing)),
        "d02_coverage": {m: {"events": int(len(a)), "birds": int(a["bird"].nunique())} for m, (a, _, _) in d02_loaded.items()},
    }
    (out / "validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    if not synthetic_pass:
        raise RuntimeError("Synthetic direction test failed")
    if sorted(u_anchor["delta_m"].unique().astype(int).tolist()) != U_SCALES:
        raise RuntimeError("Uesaka scale family mismatch")
    if args.validation_only:
        config = {
            "status": "VALIDATION_COMPLETE",
            "cpu_affinity": aff,
            "cpu1_excluded": 1 not in aff,
            "validation_only": True,
            "rd_rerun": False,
            "network_used": False,
            "paper1_modified": False,
            "elapsed_s": time.monotonic() - started,
        }
        (out / "run_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
        write_progress(out, "complete", 1, 1, started, "validation only; no real signed-lag effects computed")
        return

    u_rows = []
    u_phase_all = []
    u_bin_all = []
    u_event_all = []
    for si, scale in enumerate(U_SCALES):
        anchors = u_anchor[u_anchor["delta_m"] == scale].copy()
        flagged = attach_actual_flags(anchors, u_events, 600.0)
        ps, phase, residual = phase_test(
            flagged, u_events, u_spans, 600.0, args.phase_reps, SEED + scale,
            progress_cb=lambda c, t, s=scale, i=si: write_progress(out, "uesaka_phase", i * args.phase_reps + c, len(U_SCALES) * args.phase_reps, started, f"delta={s}"),
        )
        row = summarize_test(flagged, ps, residual, args.bootstrap_reps, SEED + 10000 + scale)
        row["delta_m"] = scale
        u_rows.append(row)
        phase.insert(0, "delta_m", scale)
        u_phase_all.append(phase)
        bins = bin_summary(flagged, u_events, [(0, 30), (30, 120), (120, 300), (300, 600)], args.bootstrap_reps, SEED + 20000 + scale)
        bins.insert(0, "delta_m", scale)
        u_bin_all.append(bins)
        u_event_all.append(flagged)
    us = pd.DataFrame(u_rows).sort_values("delta_m").reset_index(drop=True)
    us["holm15_direction_p"] = holm(us["p_diff_upper"].tolist())
    side_p = us["p_pre_upper"].tolist() + us["p_post_upper"].tolist()
    side_h = holm(side_p)
    us["holm30_pre_p"] = side_h[: len(us)]
    us["holm30_post_p"] = side_h[len(us) :]
    us["direction_pass"] = (us["raw_diff_ci_low"] > 0) & (us["resid_diff_ci_low"] > 0) & (us["holm15_direction_p"] <= 0.05)
    us["pre_side_pass"] = (us["resid_pre_ci_low"] > 0) & (us["holm30_pre_p"] <= 0.05)
    us["post_side_pass"] = (us["resid_post_ci_low"] > 0) & (us["holm30_post_p"] <= 0.05)
    us.to_csv(out / "uesaka_scale_direction_summary.csv", index=False)
    pd.concat(u_phase_all, ignore_index=True).to_csv(out / "uesaka_phase_null.csv.gz", index=False, compression="gzip")
    pd.concat(u_bin_all, ignore_index=True).to_csv(out / "uesaka_lag_bin_summary.csv", index=False)
    pd.concat(u_event_all, ignore_index=True).to_csv(out / "uesaka_endpoint_signed_lags.csv.gz", index=False, compression="gzip")

    d09_flag = attach_actual_flags(d09_anchor, d09_events, 3600.0)
    d09_ps, d09_phase, d09_resid = phase_test(
        d09_flag, d09_events, d09_spans, 3600.0, args.phase_reps, SEED + 30000,
        progress_cb=lambda c, t: write_progress(out, "d09_phase", c, t, started, "landing-ingestion"),
    )
    d09_row = summarize_test(d09_flag, d09_ps, d09_resid, args.bootstrap_reps, SEED + 31000)
    d09_row["raw_post_dominant_pass"] = bool(d09_row["raw_diff_ci_low"] > 0)
    d09_row["null_centered_post_dominant_pass"] = bool(d09_row["resid_diff_ci_low"] > 0 and d09_row["p_diff_upper"] <= 0.05)
    pd.DataFrame([d09_row]).to_csv(out / "d09_landing_ingestion_direction_summary.csv", index=False)
    d09_phase.to_csv(out / "d09_phase_null.csv.gz", index=False, compression="gzip")
    bin_summary(d09_flag, d09_events, [(0, 300), (300, 900), (900, 3600)], args.bootstrap_reps, SEED + 32000).to_csv(out / "d09_lag_bin_summary.csv", index=False)
    d09_bout.to_csv(out / "d09_wet_bout_ingestion_lags.csv.gz", index=False, compression="gzip")
    d09_flag.to_csv(out / "d09_landing_signed_lags.csv.gz", index=False, compression="gzip")

    d02_rows = []
    d02_phase_all = []
    d02_bin_all = []
    d02_methods = ["RD100", "RD200"]
    d02_estimable = [method for method in d02_methods if len(d02_loaded[method][0]) > 0]
    d02_progress_index = {method: i for i, method in enumerate(d02_estimable)}
    for mi, method in enumerate(d02_methods):
        a, ev, spans = d02_loaded[method]
        if len(a) == 0:
            d02_rows.append({
                "method": method,
                "status": "not_estimable_bilateral_15m_coverage_zero",
                "events": 0,
                "birds": 0,
                "p_ingestion_after_endpoint_upper": np.nan,
            })
            continue
        # anchor=ingestion; pre endpoint means ingestion follows endpoint. Rename in output.
        flag = attach_actual_flags(a, ev, 900.0)
        ps, phase, resid = phase_test(
            flag, ev, spans, 900.0, args.phase_reps, SEED + 40000 + mi,
            progress_cb=lambda c, t, m=method, i=d02_progress_index[method]: write_progress(
                out,
                "d02_phase",
                i * args.phase_reps + c,
                len(d02_estimable) * args.phase_reps,
                started,
                m,
            ),
        )
        row = summarize_test(flag, ps, resid, args.bootstrap_reps, SEED + 41000 + mi)
        row["method"] = method
        row["status"] = "estimated"
        row["ingestion_after_endpoint"] = row.pop("actual_pre")
        row["ingestion_before_endpoint"] = row.pop("actual_post")
        row["after_minus_before"] = -row.pop("actual_diff")
        row["null_after_minus_before"] = -row.pop("null_diff_mean")
        row["resid_after_minus_before"] = -row.pop("resid_diff")
        # summarize_test uses native event-after-anchor minus event-before-anchor.
        # Here the anchor is ingestion, so negate and swap CI bounds to report the
        # biologically readable ingestion-after-endpoint minus ingestion-before-endpoint.
        native_raw_low = row.pop("raw_diff_ci_low")
        native_raw_high = row.pop("raw_diff_ci_high")
        native_resid_low = row.pop("resid_diff_ci_low")
        native_resid_high = row.pop("resid_diff_ci_high")
        row["after_minus_before_ci_low"] = -native_raw_high
        row["after_minus_before_ci_high"] = -native_raw_low
        row["resid_after_minus_before_ci_low"] = -native_resid_high
        row["resid_after_minus_before_ci_high"] = -native_resid_low
        row["p_ingestion_after_endpoint_upper"] = float((1 + np.sum((-phase["diff"].to_numpy()) >= row["after_minus_before"])) / (args.phase_reps + 1))
        d02_rows.append(row)
        phase.insert(0, "method", method)
        phase["ingestion_after_minus_before"] = -phase["diff"]
        d02_phase_all.append(phase)
        bins = bin_summary(flag, ev, [(0, 300), (300, 900)], args.bootstrap_reps, SEED + 42000 + mi)
        bins.insert(0, "method", method)
        bins["ingestion_after_minus_before"] = -bins["post_minus_pre"]
        d02_bin_all.append(bins)
    d02s = pd.DataFrame(d02_rows)
    d02s["holm2_after_p"] = holm(d02s["p_ingestion_after_endpoint_upper"].tolist())
    d02s["after_direction_pass"] = (
        (d02s.get("after_minus_before_ci_low", pd.Series(np.nan, index=d02s.index)) > 0)
        & (d02s.get("resid_after_minus_before_ci_low", pd.Series(np.nan, index=d02s.index)) > 0)
        & (d02s["holm2_after_p"] <= 0.05)
    )
    d02s.to_csv(out / "d02_rd_ingestion_direction_summary.csv", index=False)
    if d02_phase_all:
        pd.concat(d02_phase_all, ignore_index=True).to_csv(out / "d02_phase_null.csv.gz", index=False, compression="gzip")
    else:
        pd.DataFrame().to_csv(out / "d02_phase_null.csv.gz", index=False, compression="gzip")
    if d02_bin_all:
        pd.concat(d02_bin_all, ignore_index=True).to_csv(out / "d02_lag_bin_summary.csv", index=False)
    else:
        pd.DataFrame().to_csv(out / "d02_lag_bin_summary.csv", index=False)

    config = {
        "status": "COMPLETE",
        "cpu_affinity": aff,
        "cpu1_excluded": 1 not in aff,
        "bootstrap_reps": args.bootstrap_reps,
        "phase_reps": args.phase_reps,
        "uesaka_events": int(len(u_anchor)),
        "uesaka_scales": U_SCALES,
        "d09_landings": int(len(d09_anchor)),
        "d09_ingestions": int(len(d09_ing)),
        "d02_eligible": {m: int(len(a)) for m, (a, _, _) in d02_loaded.items()},
        "prereg_sha256": sha256(PREREG),
        "script_sha256": sha256(Path(__file__)),
        "input_sha256": {str(p.relative_to(ROOT)): h for p, h in EXPECTED_HASHES.items()},
        "elapsed_s": time.monotonic() - started,
        "rd_rerun": False,
        "network_used": False,
        "paper1_modified": False,
    }
    (out / "run_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    write_progress(out, "complete", 1, 1, started, f"uesaka={len(u_anchor)} d09={len(d09_anchor)}")


if __name__ == "__main__":
    main()
