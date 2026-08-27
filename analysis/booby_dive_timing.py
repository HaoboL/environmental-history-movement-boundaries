#!/usr/bin/env python3
"""D126: RFBO D104 last-passage component by pre/post dive direction."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
import last_record_decomposition as last_record  # noqa: E402
import booby_behavior_context as behavior_context  # noqa: E402


PREREG = ROOT / "metadata/booby_dive_timing_specification_cn.md"
D49 = ROOT / "data/audit_inputs/rfbo_classified_event_features.csv.gz"
DEFAULT_OUT = ROOT / "results/booby_dive_timing"
SCALES = (250, 500, 1000, 2000)
GROUPS = ("none_observed", "pre_only", "post_only", "both", "coincident")
PRE = GROUPS.index("pre_only")
POST = GROUPS.index("post_only")
NONE = GROUPS.index("none_observed")
METRICS = last_record.A_NAMES
SEED = 1260826


def finite(value: Any) -> Any:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite(item) for item in value]
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def holm(values: list[float]) -> list[float]:
    raw = np.asarray(values, float)
    order = np.argsort(np.where(np.isfinite(raw), raw, 1.0), kind="stable")
    out = np.ones(len(raw), float)
    running = 0.0
    for rank, index in enumerate(order):
        value = raw[index] if np.isfinite(raw[index]) else 1.0
        running = max(running, min(1.0, (len(raw) - rank) * value))
        out[index] = running
    return out.tolist()


class Progress:
    def __init__(self, output: Path) -> None:
        self.output = output
        self.path = output / "progress.jsonl"
        self.started = time.monotonic()

    def emit(self, stage: str, done: int, total: int, **extra: Any) -> None:
        if stage == "rfbo_reconstruct" and done != total and done % 50 != 0:
            return
        elapsed = max(time.monotonic() - self.started, 1e-9)
        row = {
            "stage": stage, "completed": int(done), "total": int(total),
            "percent": round(100 * done / max(total, 1), 3), "elapsed_s": round(elapsed, 3),
            "eta_s": round((total - done) * elapsed / done, 3) if done and done < total else 0.0,
            **extra,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(finite(row), ensure_ascii=False) + "\n")
        print(json.dumps(finite(row), ensure_ascii=False), flush=True)


def bootstrap(unit: np.ndarray, reps: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    draws = unit[rng.integers(0, len(unit), size=(reps, len(unit)))].mean(axis=1)
    return np.quantile(draws, 0.025, axis=0), np.quantile(draws, 0.975, axis=0)


def phase_p_two_sided(observed: float, null: np.ndarray) -> float:
    center = float(null.mean())
    return float((1 + np.sum(np.abs(null - center) >= abs(observed - center))) / (len(null) + 1))


def attach_groups(records: list[dict[str, Any]]) -> tuple[pd.DataFrame, dict[str, Any]]:
    labels = pd.read_csv(D49, usecols=["delta_m", "event_id", "group"], low_memory=False)
    labels = labels.loc[labels.delta_m.isin(SCALES)].drop_duplicates()
    if labels.duplicated(["delta_m", "event_id"]).any():
        raise RuntimeError("D49 event key is not unique")
    group_map = labels.set_index(["delta_m", "event_id"]).group.to_dict()
    rows = []
    for record in records:
        key = (int(record["scale_m"]), str(record["event_id"]))
        if key not in group_map:
            raise RuntimeError(f"missing D49 group: {key}")
        group = str(group_map[key])
        if group not in GROUPS:
            raise RuntimeError(f"unexpected D49 group: {group}")
        record["direction_group"] = group
        rows.append({
            "scale_m": key[0], "event_id": key[1], "cluster": str(record["cluster"]),
            "segment_key": str(record["segment_key"]), "direction_group": group,
        })
    inventory = pd.DataFrame(rows)
    audit = {
        "records": len(records), "labels": len(labels), "exact_join": True,
        "groups": inventory.direction_group.value_counts().to_dict(),
    }
    return inventory, audit


def shifted_segment_sums(
    records: list[dict[str, Any]], reps: int, rng: np.random.Generator,
) -> tuple[np.ndarray, int, int]:
    group_index = {name: index for index, name in enumerate(GROUPS)}
    master: list[float] = []
    descriptors: list[tuple[int, int, np.ndarray, int]] = []
    for record in records:
        values = np.asarray(record["values"], float)
        master.extend(values.tolist())
        descriptors.append((
            len(master) - 1, len(values), np.asarray(record["position"], int),
            group_index[str(record["direction_group"])],
        ))
    tokens = np.asarray(master, float)
    cache: dict[int, np.ndarray | None] = {}
    output = np.zeros((reps, len(GROUPS), len(METRICS)), float)
    attempts = 0
    for replicate in range(reps):
        while True:
            attempts += 1
            if attempts > reps * max(1000, 20 * len(tokens)):
                raise RuntimeError("phase sampler exhausted attempts")
            offset = int(rng.integers(1, len(tokens)))
            if offset not in cache:
                sums = np.zeros((len(GROUPS), len(METRICS)), float)
                valid = True
                for end, width, positions, category in descriptors:
                    pseudo_end = (end - offset) % len(tokens)
                    indices = (pseudo_end - np.arange(width - 1, -1, -1, dtype=int)) % len(tokens)
                    window = tokens[indices]
                    if np.max(window) <= np.min(window):
                        valid = False
                        break
                    sums[category] += last_record.metrics_a(window, positions)
                cache[offset] = sums if valid else None
            if cache[offset] is not None:
                output[replicate] = np.asarray(cache[offset], float)
                break
    return output, len(tokens), len(cache)


def contrast_rows(
    name: str, left: int, right: int, observed_sum: dict[str, np.ndarray],
    null_sum: dict[str, np.ndarray], counts: dict[str, np.ndarray],
    phase_reps: int, bootstrap_reps: int, rng: np.random.Generator,
) -> list[dict[str, Any]]:
    paired = sorted(cluster for cluster in counts if counts[cluster][left] and counts[cluster][right])
    unit = np.row_stack([
        observed_sum[cluster][left] / counts[cluster][left]
        - observed_sum[cluster][right] / counts[cluster][right]
        for cluster in paired
    ])
    null = np.mean(np.stack([
        null_sum[cluster][:, left, :] / counts[cluster][left]
        - null_sum[cluster][:, right, :] / counts[cluster][right]
        for cluster in paired
    ], axis=1), axis=1)
    observed = unit.mean(axis=0)
    low, high = bootstrap(unit, bootstrap_reps, rng)
    rows = []
    for index, metric in enumerate(METRICS):
        rows.append({
            "contrast": name, "left": GROUPS[left], "right": GROUPS[right],
            "estimand": metric, "component": metric.split("_")[0], "tail": metric.split("_")[1],
            "paired_birds": len(paired),
            "left_events": int(sum(counts[cluster][left] for cluster in paired)),
            "right_events": int(sum(counts[cluster][right] for cluster in paired)),
            "observed_unit_equal": float(observed[index]),
            "bootstrap_ci_low": float(low[index]), "bootstrap_ci_high": float(high[index]),
            "phase_null_mean": float(null[:, index].mean()),
            "phase_p_two_sided": phase_p_two_sided(float(observed[index]), null[:, index]),
            "phase_reps": phase_reps, "bootstrap_reps": bootstrap_reps,
        })
    return rows


def scale_analysis(
    records: list[dict[str, Any]], scale: int, phase_reps: int, bootstrap_reps: int,
    rng: np.random.Generator, progress: Progress,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    use = [record for record in records if int(record["scale_m"]) == scale]
    by_segment: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in use:
        by_segment[str(record["segment_key"])].append(record)
    observed_sum: dict[str, np.ndarray] = defaultdict(lambda: np.zeros((len(GROUPS), len(METRICS)), float))
    null_sum: dict[str, np.ndarray] = defaultdict(lambda: np.zeros((phase_reps, len(GROUPS), len(METRICS)), float))
    counts: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(len(GROUPS), int))
    phase_audit = []
    group_index = {name: index for index, name in enumerate(GROUPS)}
    segments = list(by_segment.items())
    for number, (segment, segment_records) in enumerate(segments, 1):
        clusters = {str(record["cluster"]) for record in segment_records}
        if len(clusters) != 1:
            raise RuntimeError("segment maps to multiple birds")
        cluster = next(iter(clusters))
        for record in segment_records:
            category = group_index[str(record["direction_group"])]
            observed_sum[cluster][category] += np.asarray(record["observed"], float)
            counts[cluster][category] += 1
        shifted, tokens, unique_offsets = shifted_segment_sums(segment_records, phase_reps, rng)
        null_sum[cluster] += shifted
        phase_audit.append({
            "scale_m": scale, "segment_key": segment, "cluster": cluster,
            "events": len(segment_records), "tokens": tokens, "unique_offsets_attempted": unique_offsets,
        })
        if number % max(1, len(segments) // 10) == 0 or number == len(segments):
            progress.emit(f"phase_{scale}", number, len(segments), events=len(use))
    rows = []
    rows.extend(contrast_rows("post_minus_pre", POST, PRE, observed_sum, null_sum, counts, phase_reps, bootstrap_reps, rng))
    rows.extend(contrast_rows("post_minus_none", POST, NONE, observed_sum, null_sum, counts, phase_reps, bootstrap_reps, rng))
    rows.extend(contrast_rows("pre_minus_none", PRE, NONE, observed_sum, null_sum, counts, phase_reps, bootstrap_reps, rng))
    paired = sorted(cluster for cluster in counts if counts[cluster][POST] and counts[cluster][PRE])
    loo_unit = np.asarray([
        (observed_sum[cluster][POST] / counts[cluster][POST] - observed_sum[cluster][PRE] / counts[cluster][PRE])[8]
        for cluster in paired
    ])
    loo = np.asarray([np.delete(loo_unit, index).mean() for index in range(len(loo_unit))])
    audit = {
        "scale_m": scale, "events": len(use), "segments": len(segments), "birds": len(counts),
        "paired_pre_post_birds": len(paired), "loo_L_union_min": float(loo.min()),
        "loo_L_union_max": float(loo.max()), "loo_fraction_positive": float(np.mean(loo > 0)),
    }
    return rows, phase_audit, audit


def adjudicate(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame["holm_family"] = "diagnostic"
    frame["holm_phase_p"] = np.nan
    primary = frame.contrast.eq("post_minus_pre") & frame.estimand.eq("L_union")
    frame.loc[primary, "holm_family"] = "primary_post_minus_pre_L_union_4_scales"
    frame.loc[primary, "holm_phase_p"] = holm(frame.loc[primary, "phase_p_two_sided"].tolist())
    frame["coverage_gate"] = frame.paired_birds.ge(15)
    frame["primary_pass"] = False
    frame.loc[primary, "primary_pass"] = (
        frame.loc[primary, "coverage_gate"]
        & ((frame.loc[primary, "bootstrap_ci_low"] > 0) | (frame.loc[primary, "bootstrap_ci_high"] < 0))
        & frame.loc[primary, "holm_phase_p"].le(0.05)
    )
    primary_rows = frame.loc[primary].sort_values("scale_m")
    passed = primary_rows.loc[primary_rows.primary_pass, "scale_m"].astype(int).tolist()
    adjacent_pairs = []
    for left, right in zip(passed, passed[1:]):
        if right == left * 2:
            left_sign = np.sign(float(primary_rows.loc[primary_rows.scale_m.eq(left), "observed_unit_equal"].iloc[0]))
            right_sign = np.sign(float(primary_rows.loc[primary_rows.scale_m.eq(right), "observed_unit_equal"].iloc[0]))
            if left_sign == right_sign and left_sign != 0:
                adjacent_pairs.append([left, right])
    if adjacent_pairs:
        signs = [np.sign(float(primary_rows.loc[primary_rows.scale_m.eq(pair[0]), "observed_unit_equal"].iloc[0])) for pair in adjacent_pairs]
        verdict = "PRE_DIVE_LAST_PASSAGE_ENHANCED" if any(sign > 0 for sign in signs) else "POST_DIVE_PATH_RESET_ENHANCED"
    else:
        verdict = "NO_TEMPORAL_DIRECTION_SPECIFICITY"
    return frame, {
        "primary_pass_scales_m": passed, "same_sign_adjacent_pairs_m": adjacent_pairs,
        "verdict": verdict,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--phase-reps", type=int, default=999)
    parser.add_argument("--bootstrap-reps", type=int, default=20000)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        args.phase_reps = min(args.phase_reps, 19)
        args.bootstrap_reps = min(args.bootstrap_reps, 200)
    args.output.mkdir(parents=True, exist_ok=True)
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 is forbidden; affinity={affinity}")
    if not PREREG.exists():
        raise RuntimeError("D126 preregistration missing")
    progress = Progress(args.output)
    progress.emit("start", 0, 1, affinity=affinity, smoke=args.smoke)
    records, _, anchor = behavior_context.load_and_reconstruct(progress, False)
    inventory, join_audit = attach_groups(records)
    inventory.to_csv(args.output / "event_inventory.csv", index=False)

    rng = np.random.default_rng(SEED)
    rows = []
    phase_rows = []
    audits = []
    for number, scale in enumerate(SCALES, 1):
        part, phase, audit = scale_analysis(records, scale, args.phase_reps, args.bootstrap_reps, rng, progress)
        for row in part:
            row["scale_m"] = scale
        rows.extend(part); phase_rows.extend(phase); audits.append(audit)
        progress.emit("scale_complete", number, len(SCALES), scale_m=scale)
    results, verdict = adjudicate(pd.DataFrame(rows))
    results.to_csv(args.output / "direction_contrasts.csv", index=False)
    pd.DataFrame(phase_rows).to_csv(args.output / "phase_null_audit.csv.gz", index=False, compression="gzip")
    pd.DataFrame(audits).to_csv(args.output / "scale_audit.csv", index=False)
    summary = {
        "status": "D126_SMOKE_COMPLETE" if args.smoke else "D126_FORMAL_COMPLETE",
        "seed": SEED, "phase_reps": args.phase_reps, "bootstrap_reps": args.bootstrap_reps,
        "affinity": affinity, "cpu1_excluded": 1 not in affinity,
        "anchor_audit": anchor, "join_audit": join_audit, "scale_audit": audits, "verdict": verdict,
        "input_sha256": {
            "prereg": sha256(PREREG), "d49": sha256(D49), "last_record_events": sha256(behavior_context.D104_EVENTS),
            "behavior_events": sha256(behavior_context.D28_EVENTS), "last_record_script": sha256(Path(last_record.__file__)),
            "this_script": sha256(Path(__file__)),
        },
        "interpretation_limit": "TDR dive attempt timing around a last-passage event; not capture, sensory cue, prey truth, or Levy generation.",
    }
    (args.output / "summary.json").write_text(json.dumps(finite(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    progress.emit("complete", 1, 1, verdict=verdict["verdict"])


if __name__ == "__main__":
    main()
