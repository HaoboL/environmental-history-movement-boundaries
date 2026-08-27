#!/usr/bin/env python3
"""Independent D107 primary-gate audit; no RD detection and no source-script import."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
FORMAL = ROOT / "results/shearwater_behavior"
AUDIT_OUT = ROOT / "output/audits/shearwater_behavior"
SCALES = (500, 1000, 2000, 5000)
REPS = 999
BOOTS = 50000
SEED = 1070826


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def metrics(values: np.ndarray, record_runs: np.ndarray) -> np.ndarray:
    rank = pd.Series(values).rank(method="average").to_numpy(float) / len(values)
    flags = np.column_stack([rank > .8, rank <= .2, (rank > .8) | (rank <= .2)]).astype(float)
    endpoint = flags[-1]; overall = flags.mean(axis=0); record = flags[record_runs].mean(axis=0)
    return np.r_[endpoint - overall, record - overall, endpoint - record]


def holm(values: list[float]) -> list[float]:
    raw = np.asarray(values, float); order = np.argsort(raw, kind="stable")
    output = np.ones(len(raw)); running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(raw) - rank) * raw[index])); output[index] = running
    return output.tolist()


def slope_units(frame: pd.DataFrame, matrix: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    units = []; selectors = []
    for bird, group in frame.groupby("individual", sort=True):
        idx = group.matrix_index.to_numpy(int); x = group.future_specificity.to_numpy(float)
        if len(idx) < 3 or np.var(x) <= 0:
            continue
        centered = x - x.mean(); denom = float(centered @ centered)
        units.append(float(centered @ matrix[idx] / denom)); selectors.append((idx, centered, denom))
    return np.asarray(units, float), selectors


def group_units(frame: pd.DataFrame, matrix: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    use = frame.loc[frame.behavior_class.eq("forage_dominant")]
    units = []; selectors = []
    for bird, group in use.groupby("individual", sort=True):
        idx = group.matrix_index.to_numpy(int)
        units.append(float(matrix[idx].mean())); selectors.append(idx)
    return np.asarray(units, float), selectors


def bootstrap_ci(unit: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    draws = unit[rng.integers(0, len(unit), size=(BOOTS, len(unit)))].mean(axis=1)
    return tuple(float(v) for v in np.quantile(draws, [.025, .975]))


def main() -> None:
    AUDIT_OUT.mkdir(parents=True, exist_ok=True)
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 forbidden: {affinity}")
    events = pd.read_csv(FORMAL / "event_metrics.csv.gz", low_memory=False)
    formal_summary = json.loads((FORMAL / "summary.json").read_text(encoding="utf-8"))
    sequences = []
    with gzip.open(FORMAL / "event_sequences.jsonl.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            sequences.append(json.loads(line))
    seq = pd.DataFrame(sequences)
    for frame in (events, seq):
        frame["scale_m"] = pd.to_numeric(frame.scale_m, errors="raise").astype(int)
        frame["event_id"] = frame.event_id.astype(str)
        frame["run_key"] = frame.run_key.astype(str)
        frame["individual"] = frame.individual.astype(str)
    joined = events.merge(seq, on=["scale_m", "event_id", "run_key", "individual"], validate="one_to_one")
    if len(joined) != len(events):
        raise RuntimeError("sequence/event join failed")
    reconstructed = np.row_stack([
        metrics(np.asarray(row.values, float), np.asarray(row.record_runs, int)) for row in joined.itertuples(index=False)
    ])
    frozen = events[["E_high", "E_low", "E_union", "R_high", "R_low", "R_union", "L_high", "L_low", "L_union"]].to_numpy(float)
    max_metric_difference = float(np.max(np.abs(reconstructed - frozen)))
    if max_metric_difference > 1e-12:
        raise RuntimeError(f"metric reconstruction failed: {max_metric_difference}")
    observed_l = reconstructed[:, 8]
    joined["matrix_index"] = np.arange(len(joined))

    rng = np.random.default_rng(SEED)
    null = np.zeros((REPS, len(joined)), float)
    for (scale, run_key), group in joined.groupby(["scale_m", "run_key"], sort=False):
        descriptors = []; master = []
        for row in group.itertuples(index=False):
            value = np.asarray(row.values, float); master.extend(value.tolist())
            descriptors.append((len(master) - 1, len(value), np.asarray(row.record_runs, int), int(row.matrix_index)))
        master = np.asarray(master, float); cache = {}
        for replicate in range(REPS):
            while True:
                offset = int(rng.integers(1, len(master)))
                if offset not in cache:
                    trial = []; valid = True
                    for end, width, record_runs, index in descriptors:
                        pseudo_end = (end - offset) % len(master)
                        ix = (pseudo_end - np.arange(width - 1, -1, -1)) % len(master)
                        window = master[ix]
                        if np.max(window) <= np.min(window):
                            valid = False; break
                        trial.append((index, metrics(window, record_runs)[8]))
                    cache[offset] = trial if valid else None
                if cache[offset] is not None:
                    for index, value in cache[offset]:
                        null[replicate, index] = value
                    break

    slope_rows = []; group_rows = []
    for scale in SCALES:
        frame = joined.loc[joined.scale_m.eq(scale)]
        slope, slope_selectors = slope_units(frame, observed_l)
        slope_null = np.column_stack([
            (null[:, idx] * centered[None, :]).sum(axis=1) / denom for idx, centered, denom in slope_selectors
        ]).mean(axis=1)
        slope_ci = bootstrap_ci(slope, rng)
        slope_rows.append({
            "scale_m": scale, "events": int(sum(len(x[0]) for x in slope_selectors)),
            "individuals": len(slope), "effect": float(slope.mean()),
            "ci_low": slope_ci[0], "ci_high": slope_ci[1],
            "phase_p_raw": float((1 + np.sum(slope_null >= slope.mean())) / (REPS + 1)),
        })
        group, group_selectors = group_units(frame, observed_l)
        group_null = np.column_stack([null[:, idx].mean(axis=1) for idx in group_selectors]).mean(axis=1)
        ci = bootstrap_ci(group, rng)
        group_rows.append({
            "scale_m": scale, "events": int(sum(len(idx) for idx in group_selectors)),
            "individuals": len(group), "effect": float(group.mean()),
            "ci_low": ci[0], "ci_high": ci[1],
            "phase_p_raw": float((1 + np.sum(group_null >= group.mean())) / (REPS + 1)),
        })
    slopes = pd.DataFrame(slope_rows); groups = pd.DataFrame(group_rows)
    slopes["phase_p_holm"] = holm(slopes.phase_p_raw.tolist())
    groups["phase_p_holm"] = holm(groups.phase_p_raw.tolist())
    slopes["pass"] = (slopes.events >= 30) & (slopes.individuals >= 6) & (slopes.effect > 0) & (slopes.ci_low > 0) & (slopes.phase_p_holm <= .05)
    groups["pass"] = (groups.events >= 15) & (groups.individuals >= 5) & (groups.effect > 0) & (groups.ci_low > 0) & (groups.phase_p_holm <= .05)
    joint = sorted(set(slopes.loc[slopes["pass"], "scale_m"]) & set(groups.loc[groups["pass"], "scale_m"]))
    adjacent = [(a, b) for a, b in zip(SCALES[:-1], SCALES[1:]) if a in joint and b in joint]
    verdict = "DIRECT_FORAGING_STATE_LINK_SUPPORTED" if adjacent else "NO_DIRECT_FORAGING_STATE_LINK"
    if verdict != formal_summary["status"]:
        raise RuntimeError(f"audit verdict mismatch: {verdict} != {formal_summary['status']}")
    result = {
        "status": "D107_FORMAL_V1_AUDIT_PASS", "audited_verdict": verdict,
        "seed": SEED, "phase_reps": REPS, "bootstrap_reps": BOOTS,
        "cpu_affinity": affinity, "events": len(events),
        "maximum_reconstructed_metric_difference": max_metric_difference,
        "slope_pass_scales": slopes.loc[slopes["pass"], "scale_m"].astype(int).tolist(),
        "group_pass_scales": groups.loc[groups["pass"], "scale_m"].astype(int).tolist(),
        "joint_adjacent_pairs": adjacent,
        "input_hashes": {
            "event_metrics": sha256(FORMAL / "event_metrics.csv.gz"),
            "event_sequences": sha256(FORMAL / "event_sequences.jsonl.gz"),
            "formal_summary": sha256(FORMAL / "summary.json"),
            "audit_script": sha256(Path(__file__)),
        },
    }
    slopes.to_csv(AUDIT_OUT / "audit_primary_slopes.csv", index=False)
    groups.to_csv(AUDIT_OUT / "audit_primary_groups.csv", index=False)
    (AUDIT_OUT / "audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(slopes.to_string(index=False)); print(groups.to_string(index=False)); print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
