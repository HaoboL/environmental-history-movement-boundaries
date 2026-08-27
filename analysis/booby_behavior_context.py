#!/usr/bin/env python3
"""D106: direct RFBO last-passage component x feeding-intent test.

The analysis is locked in
NODE_STATE/D106_RFBO_LAST_PASSAGE_FEEDING_INTENT_DIRECT_PREREG_20260825_CN.md.
It reconstructs D104 Family A sequences but never reruns RD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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


PREREG = ROOT / "metadata/booby_behavior_context_specification_cn.md"
D104_EVENTS = ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz"
D28_EVENTS = ROOT / "data/audit_inputs/rfbo_event_two_tail_metrics.csv.gz"
DEFAULT_OUT = ROOT / "results/booby_behavior_context"
SCALES = (250, 500, 1000, 2000)
BEHAVIORS = ("dive_near", "wet_only_near", "no_tdr_near")
CONTRASTS = (
    ("dive_minus_wet", "dive_near", "wet_only_near"),
    ("dive_minus_no_tdr", "dive_near", "no_tdr_near"),
)
METRICS = last_record.A_NAMES
SEED = 1060825


def finite(value: Any) -> Any:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite(v) for v in value]
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
        if stage == "rfbo_reconstruct" and done != total and done % 25 != 0:
            return
        elapsed = max(time.monotonic() - self.started, 1e-9)
        rate = done / elapsed if done else 0.0
        row = {
            "stage": stage,
            "completed": int(done),
            "total": int(total),
            "percent": round(100 * done / max(total, 1), 3),
            "elapsed_s": round(elapsed, 3),
            "throughput_per_s": round(rate, 5),
            "eta_s": round((total - done) / rate, 3) if rate > 0 and done < total else 0.0,
            **extra,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(finite(row), ensure_ascii=False) + "\n")
        print(json.dumps(finite(row), ensure_ascii=False), flush=True)


def bootstrap_summary(unit: np.ndarray, reps: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if unit.ndim != 2 or len(unit) < 2:
        raise RuntimeError("bootstrap unit matrix requires at least two units")
    observed = unit.mean(axis=0)
    draws = unit[rng.integers(0, len(unit), size=(reps, len(unit)))].mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975], axis=0)
    return observed, low, high


def load_and_reconstruct(progress: Progress, validate_only: bool) -> tuple[list[dict[str, Any]], pd.DataFrame, dict[str, Any]]:
    records, _, reconstruction_audit = last_record.rfbo_records(progress, smoke=False)
    if any(record["scale_m"] not in SCALES for record in records):
        raise RuntimeError("unexpected RFBO scale")

    behavior = pd.read_csv(D28_EVENTS, usecols=["event_id", "delta_m", "behavior_class"], low_memory=False)
    behavior["event_id"] = behavior.event_id.astype(str)
    behavior["scale_m"] = pd.to_numeric(behavior.delta_m, errors="raise").astype(int)
    behavior = behavior.drop(columns="delta_m").drop_duplicates(["scale_m", "event_id"])
    if behavior.duplicated(["scale_m", "event_id"]).any():
        raise RuntimeError("D28 behavior key is not unique")
    behavior_map = behavior.set_index(["scale_m", "event_id"]).behavior_class.to_dict()

    rows: list[dict[str, Any]] = []
    for record in records:
        key = (int(record["scale_m"]), str(record["event_id"]))
        if key not in behavior_map:
            raise RuntimeError(f"D104 event missing D28 label: {key}")
        label = str(behavior_map[key])
        if label not in BEHAVIORS:
            raise RuntimeError(f"unknown behavior label: {label}")
        record["behavior_class"] = label
        row = {
            "dataset": "rfbo", "scale_m": key[0], "event_id": key[1],
            "segment_key": str(record["segment_key"]), "cluster": str(record["cluster"]),
            "behavior_class": label,
        }
        row.update({name: float(value) for name, value in zip(METRICS, record["observed"])})
        rows.append(row)
    reconstructed = pd.DataFrame(rows)

    frozen = pd.read_csv(D104_EVENTS, low_memory=False)
    frozen = frozen.loc[frozen.dataset.eq("rfbo") & frozen.family.eq("A")].copy()
    frozen["event_id"] = frozen.event_id.astype(str)
    anchor = reconstructed.merge(
        frozen[["scale_m", "event_id", "segment_key", "cluster", *METRICS]],
        on=["scale_m", "event_id"], how="outer", validate="one_to_one", indicator=True,
        suffixes=("_new", "_frozen"),
    )
    if not anchor._merge.eq("both").all():
        raise RuntimeError(f"D104 anchor mismatch: {anchor._merge.value_counts().to_dict()}")
    metadata_mismatch = int(
        anchor.segment_key_new.ne(anchor.segment_key_frozen).sum()
        + anchor.cluster_new.ne(anchor.cluster_frozen).sum()
    )
    maximum_metric_difference = max(
        float(np.max(np.abs(anchor[f"{name}_new"] - anchor[f"{name}_frozen"]))) for name in METRICS
    )
    identity_error = float(np.max(np.abs(
        reconstructed[["E_high", "E_low", "E_union"]].to_numpy(float)
        - reconstructed[["R_high", "R_low", "R_union"]].to_numpy(float)
        - reconstructed[["L_high", "L_low", "L_union"]].to_numpy(float)
    )))
    if metadata_mismatch or maximum_metric_difference > 1e-12 or identity_error > 1e-12:
        raise RuntimeError("reconstructed D104 anchor failed")

    inventory = reconstructed.groupby(["scale_m", "behavior_class"], as_index=False).agg(
        events=("event_id", "size"), birds=("cluster", "nunique"), segments=("segment_key", "nunique"),
    )
    paired_rows = []
    for scale, frame in reconstructed.groupby("scale_m", sort=True):
        tab = frame.groupby(["cluster", "behavior_class"]).size().unstack(fill_value=0)
        paired_rows.append({
            "scale_m": int(scale),
            "paired_dive_wet_birds": int(((tab.dive_near > 0) & (tab.wet_only_near > 0)).sum()),
            "paired_dive_no_tdr_birds": int(((tab.dive_near > 0) & (tab.no_tdr_near > 0)).sum()),
        })
    paired = pd.DataFrame(paired_rows)
    inventory = inventory.merge(paired, on="scale_m", how="left", validate="many_to_one")
    audit = {
        "reconstructed_events": int(len(reconstructed)),
        "frozen_events": int(len(frozen)),
        "exact_event_join": True,
        "metadata_mismatch": metadata_mismatch,
        "maximum_metric_difference": maximum_metric_difference,
        "maximum_E_minus_R_minus_L_error": identity_error,
        "reconstruction_audit_rows": int(len(reconstruction_audit)),
        "validate_only": bool(validate_only),
    }
    return records, inventory, audit


def segment_phase_metrics(
    records: list[dict[str, Any]], reps: int, rng: np.random.Generator,
) -> tuple[np.ndarray, int, int, int]:
    """Return reps x behavior x metric sums under one common phase per segment."""
    master: list[float] = []
    descriptors: list[tuple[int, int, np.ndarray, int]] = []
    behavior_index = {name: index for index, name in enumerate(BEHAVIORS)}
    for record in records:
        values = np.asarray(record["values"], float)
        master.extend(values.tolist())
        descriptors.append((
            len(master) - 1, len(values), np.asarray(record["position"], int),
            behavior_index[str(record["behavior_class"])],
        ))
    array = np.asarray(master, float)
    if len(array) < 2:
        raise RuntimeError("segment phase master has fewer than two tokens")
    cache: dict[int, np.ndarray | None] = {}
    output = np.zeros((reps, len(BEHAVIORS), len(METRICS)), float)
    attempts = 0
    for replicate in range(reps):
        while True:
            attempts += 1
            if attempts > reps * max(1000, 20 * len(array)):
                raise RuntimeError("could not sample valid common phase")
            offset = int(rng.integers(1, len(array)))
            if offset not in cache:
                sums = np.zeros((len(BEHAVIORS), len(METRICS)), float)
                valid = True
                for end, width, positions, category in descriptors:
                    pseudo_end = (end - offset) % len(array)
                    index = (pseudo_end - np.arange(width - 1, -1, -1, dtype=int)) % len(array)
                    window = array[index]
                    if np.max(window) <= np.min(window):
                        valid = False
                        break
                    sums[category] += last_record.metrics_a(window, positions)
                cache[offset] = sums if valid else None
            value = cache[offset]
            if value is not None:
                output[replicate] = value
                break
    valid_offsets = int(sum(value is not None for value in cache.values()))
    return output, len(array), len(cache), valid_offsets


def scale_effects(
    records: list[dict[str, Any]], scale: int, phase_reps: int, bootstrap_reps: int,
    rng: np.random.Generator, progress: Progress,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    use = [record for record in records if int(record["scale_m"]) == scale]
    by_segment: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in use:
        by_segment[str(record["segment_key"])].append(record)

    observed_sum: dict[str, np.ndarray] = defaultdict(lambda: np.zeros((len(BEHAVIORS), len(METRICS)), float))
    null_sum: dict[str, np.ndarray] = defaultdict(lambda: np.zeros((phase_reps, len(BEHAVIORS), len(METRICS)), float))
    counts: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(len(BEHAVIORS), int))
    null_audit: list[dict[str, Any]] = []
    behavior_index = {name: index for index, name in enumerate(BEHAVIORS)}
    segments = list(by_segment.items())
    for number, (segment, segment_records) in enumerate(segments, 1):
        clusters = {str(record["cluster"]) for record in segment_records}
        if len(clusters) != 1:
            raise RuntimeError(f"segment maps to multiple birds: {segment}")
        cluster = next(iter(clusters))
        for record in segment_records:
            category = behavior_index[str(record["behavior_class"])]
            observed_sum[cluster][category] += np.asarray(record["observed"], float)
            counts[cluster][category] += 1
        shifted, tokens, attempted_offsets, valid_offsets = segment_phase_metrics(segment_records, phase_reps, rng)
        null_sum[cluster] += shifted
        null_audit.append({
            "scale_m": scale, "segment_key": segment, "cluster": cluster,
            "events": len(segment_records), "master_tokens": tokens,
            "attempted_unique_offsets": attempted_offsets, "valid_unique_offsets": valid_offsets,
        })
        if number % max(1, len(segments) // 20) == 0 or number == len(segments):
            progress.emit(f"phase_{scale}", number, len(segments), events=len(use))

    clusters = sorted(counts)
    group_rows: list[dict[str, Any]] = []
    modifier_rows: list[dict[str, Any]] = []
    for behavior in BEHAVIORS:
        category = behavior_index[behavior]
        eligible = [cluster for cluster in clusters if counts[cluster][category] > 0]
        unit = np.row_stack([observed_sum[cluster][category] / counts[cluster][category] for cluster in eligible])
        null = np.mean(np.stack([
            null_sum[cluster][:, category, :] / counts[cluster][category] for cluster in eligible
        ], axis=1), axis=1)
        observed, low, high = bootstrap_summary(unit, bootstrap_reps, rng)
        for metric_index, metric in enumerate(METRICS):
            p = float((1 + np.sum(null[:, metric_index] >= observed[metric_index])) / (phase_reps + 1))
            group_rows.append({
                "scale_m": scale, "behavior_class": behavior, "estimand": metric,
                "component": metric.split("_")[0], "tail": metric.split("_")[1],
                "events": int(sum(counts[cluster][category] for cluster in eligible)),
                "birds": len(eligible), "observed_unit_equal": float(observed[metric_index]),
                "bootstrap_ci_low": float(low[metric_index]), "bootstrap_ci_high": float(high[metric_index]),
                "phase_null_mean": float(null[:, metric_index].mean()), "phase_p_one_sided": p,
                "phase_reps": phase_reps, "bootstrap_reps": bootstrap_reps,
            })

    for contrast, left, right in CONTRASTS:
        left_index = behavior_index[left]; right_index = behavior_index[right]
        paired = [cluster for cluster in clusters if counts[cluster][left_index] > 0 and counts[cluster][right_index] > 0]
        unit = np.row_stack([
            observed_sum[cluster][left_index] / counts[cluster][left_index]
            - observed_sum[cluster][right_index] / counts[cluster][right_index]
            for cluster in paired
        ])
        null = np.mean(np.stack([
            null_sum[cluster][:, left_index, :] / counts[cluster][left_index]
            - null_sum[cluster][:, right_index, :] / counts[cluster][right_index]
            for cluster in paired
        ], axis=1), axis=1)
        observed, low, high = bootstrap_summary(unit, bootstrap_reps, rng)
        for metric_index, metric in enumerate(METRICS):
            p = float((1 + np.sum(null[:, metric_index] >= observed[metric_index])) / (phase_reps + 1))
            modifier_rows.append({
                "scale_m": scale, "contrast": contrast, "left": left, "right": right,
                "estimand": metric, "component": metric.split("_")[0], "tail": metric.split("_")[1],
                "paired_birds": len(paired), "observed_unit_equal": float(observed[metric_index]),
                "bootstrap_ci_low": float(low[metric_index]), "bootstrap_ci_high": float(high[metric_index]),
                "phase_null_mean": float(null[:, metric_index].mean()), "phase_p_one_sided": p,
                "phase_reps": phase_reps, "bootstrap_reps": bootstrap_reps,
            })
    return group_rows, modifier_rows, null_audit


def add_multiplicity_and_verdict(groups: pd.DataFrame, modifiers: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    groups["holm_family"] = "descriptive"
    groups["holm_phase_p"] = np.nan
    modifiers["holm_family"] = "descriptive"
    modifiers["holm_phase_p"] = np.nan

    primary_group = groups.behavior_class.eq("dive_near") & groups.estimand.eq("L_union")
    groups.loc[primary_group, "holm_family"] = "primary_dive_L_union_4_scales"
    groups.loc[primary_group, "holm_phase_p"] = holm(groups.loc[primary_group, "phase_p_one_sided"].tolist())
    side_group = groups.behavior_class.eq("dive_near") & groups.estimand.isin(["L_high", "L_low"])
    groups.loc[side_group, "holm_family"] = "secondary_dive_L_sides_8"
    groups.loc[side_group, "holm_phase_p"] = holm(groups.loc[side_group, "phase_p_one_sided"].tolist())

    primary_modifier = modifiers.contrast.eq("dive_minus_wet") & modifiers.estimand.eq("L_union")
    modifiers.loc[primary_modifier, "holm_family"] = "primary_dive_minus_wet_L_union_4_scales"
    modifiers.loc[primary_modifier, "holm_phase_p"] = holm(modifiers.loc[primary_modifier, "phase_p_one_sided"].tolist())
    side_modifier = modifiers.contrast.eq("dive_minus_wet") & modifiers.estimand.isin(["L_high", "L_low"])
    modifiers.loc[side_modifier, "holm_family"] = "secondary_dive_minus_wet_L_sides_8"
    modifiers.loc[side_modifier, "holm_phase_p"] = holm(modifiers.loc[side_modifier, "phase_p_one_sided"].tolist())
    no_tdr = modifiers.contrast.eq("dive_minus_no_tdr") & modifiers.component.eq("L")
    modifiers.loc[no_tdr, "holm_family"] = "secondary_dive_minus_no_tdr_L_all_12"
    modifiers.loc[no_tdr, "holm_phase_p"] = holm(modifiers.loc[no_tdr, "phase_p_one_sided"].tolist())
    non_l = modifiers.contrast.eq("dive_minus_wet") & modifiers.component.isin(["E", "R"])
    modifiers.loc[non_l, "holm_family"] = "secondary_dive_minus_wet_E_R_all_24"
    modifiers.loc[non_l, "holm_phase_p"] = holm(modifiers.loc[non_l, "phase_p_one_sided"].tolist())

    groups["coverage_gate"] = groups.events.ge(100) & groups.birds.ge(30)
    modifiers["coverage_gate"] = modifiers.paired_birds.ge(30)
    groups["primary_pass"] = False
    groups.loc[primary_group, "primary_pass"] = (
        groups.loc[primary_group, "coverage_gate"]
        & groups.loc[primary_group, "observed_unit_equal"].gt(0)
        & groups.loc[primary_group, "bootstrap_ci_low"].gt(0)
        & groups.loc[primary_group, "holm_phase_p"].le(0.05)
    )
    modifiers["primary_pass"] = False
    modifiers.loc[primary_modifier, "primary_pass"] = (
        modifiers.loc[primary_modifier, "coverage_gate"]
        & modifiers.loc[primary_modifier, "observed_unit_equal"].gt(0)
        & modifiers.loc[primary_modifier, "bootstrap_ci_low"].gt(0)
        & modifiers.loc[primary_modifier, "holm_phase_p"].le(0.05)
    )

    group_pass = set(groups.loc[primary_group & groups.primary_pass, "scale_m"].astype(int))
    modifier_pass = set(modifiers.loc[primary_modifier & modifiers.primary_pass, "scale_m"].astype(int))
    joint = sorted(group_pass & modifier_pass)
    adjacent = any(left in joint and right in joint for left, right in zip(SCALES[:-1], SCALES[1:]))
    verdict = {
        "primary_dive_L_union_pass_scales": sorted(group_pass),
        "primary_dive_minus_wet_L_union_pass_scales": sorted(modifier_pass),
        "joint_pass_scales": joint,
        "adjacent_joint_pass": adjacent,
        "verdict": "DIRECT_FEEDING_INTENT_LINK_SUPPORTED" if adjacent else "NO_DIRECT_FEEDING_INTENT_SPECIFICITY",
    }
    return groups, modifiers, verdict


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--phase-reps", type=int, default=999)
    parser.add_argument("--bootstrap-reps", type=int, default=5000)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    progress = Progress(args.output)

    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 is forbidden; affinity={affinity}")
    if not PREREG.exists():
        raise RuntimeError("preregistration is missing")
    rng = np.random.default_rng(SEED)
    records, inventory, anchor_audit = load_and_reconstruct(progress, args.validate_only)
    inventory.to_csv(args.output / "event_inventory.csv", index=False)
    with (args.output / "anchor_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(finite(anchor_audit), handle, ensure_ascii=False, indent=2)
    if args.validate_only:
        progress.emit("validate_only_complete", 1, 1, events=len(records))
        return

    group_rows: list[dict[str, Any]] = []
    modifier_rows: list[dict[str, Any]] = []
    null_audit_rows: list[dict[str, Any]] = []
    for scale_number, scale in enumerate(SCALES, 1):
        g, m, n = scale_effects(records, scale, args.phase_reps, args.bootstrap_reps, rng, progress)
        group_rows.extend(g); modifier_rows.extend(m); null_audit_rows.extend(n)
        progress.emit("scale_complete", scale_number, len(SCALES), scale_m=scale)
    groups = pd.DataFrame(group_rows)
    modifiers = pd.DataFrame(modifier_rows)
    groups, modifiers, verdict = add_multiplicity_and_verdict(groups, modifiers)
    groups.to_csv(args.output / "behavior_group_effects.csv", index=False)
    modifiers.to_csv(args.output / "behavior_modifiers.csv", index=False)
    pd.DataFrame(null_audit_rows).to_csv(args.output / "phase_null_audit.csv.gz", index=False, compression="gzip")

    summary = {
        "status": "D106_FORMAL_COMPLETE",
        "seed": SEED,
        "phase_reps": args.phase_reps,
        "bootstrap_reps": args.bootstrap_reps,
        "affinity": affinity,
        "cpu1_excluded": 1 not in affinity,
        "scales_m": list(SCALES),
        "anchor_audit": anchor_audit,
        "verdict": verdict,
        "input_sha256": {
            "prereg": sha256(PREREG), "d104_events": sha256(D104_EVENTS), "d28_events": sha256(D28_EVENTS),
            "last_record_script": sha256(Path(last_record.__file__)), "this_script": sha256(Path(__file__)),
            "rfbo_gps": sha256(last_record.RFBO_GPS), "rfbo_rd": sha256(last_record.RFBO_RD),
            "deployments": sha256(last_record.DEPLOYMENTS),
        },
    }
    with (args.output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(finite(summary), handle, ensure_ascii=False, indent=2)
    progress.emit("formal_complete", 1, 1, verdict=verdict["verdict"])


if __name__ == "__main__":
    main()
