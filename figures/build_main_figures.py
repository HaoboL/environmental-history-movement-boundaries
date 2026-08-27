#!/usr/bin/env python3
"""Build Paper 2 V3 main figures from frozen result tables only.

The script performs display transformations and deterministic example
selection. It does not rerun segmentation, environmental extraction,
bootstrap, randomisation or model fitting.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import cmocean

from figure_utils import (  # noqa: E402
    add_panel_labels,
    audit_layout,
    export_figure,
    print_report,
    render_preview,
    setup_style,
)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT
OUT = ROOT / "output/figures"
D104 = ROOT / "results/last_record_decomposition"
D107 = ROOT / "results/shearwater_behavior"
D127 = ROOT / "results/laysan_same_event"
D128 = ROOT / "results/boundary_counterfactual"
D45 = ROOT / "results/landing_timing"
DENSE = ROOT / "external_inputs/goto_canonical_windows_continuous_dense_chl.csv.gz"
CANON = ROOT / "external_inputs/goto_canonical_events.csv"

BLUE = "#0072B2"          # observed last-passage evidence / relative low
ORANGE = "#D55E00"        # relative high / confirmation contrast
GREEN = "#009E73"
PURPLE = "#6A3D9A"        # combined-tail summaries only
SKY = "#56B4E9"
YELLOW = "#E69F00"        # observed boundary rho
GREY = "#737373"          # record/background comparison
LIGHT = "#D9D9D9"         # null distributions and unselected path
BLACK = "#222222"
OBSERVED = BLUE
NULL = "#8C8C8C"
RECORD = "#6F6F6F"
CONFIRMATION = "#B24745"
ROSE = "#B85C7A"          # counterfactual/reference accent
PALE_ROSE = "#F2D6DF"     # low-salience counterfactual/reference fill
PALE_BLUE = "#D7EAF4"     # low-salience observed/component fill
PALE_VIOLET = "#E9E1F0"   # total assembled from component contributions


@lru_cache(maxsize=None)
def file_sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finish(fig: plt.Figure, name: str, size: tuple[float, float],
           sources: list[Path], boundary: str, axes: list[plt.Axes]) -> dict:
    fig.canvas.draw()
    add_panel_labels(fig, axes=axes, style="nature")
    fig.canvas.draw()
    fig.set_layout_engine("none")
    issues = audit_layout(fig)
    verdict = print_report(issues)
    preview = render_preview(fig, str(OUT / f"{name}_preview.png"), dpi=180)
    exports = export_figure(
        fig, str(OUT / name), formats=["pdf", "svg", "png"],
        size_inches=size, dpi=300, grayscale_preview=False, tight=False,
    )
    grey = OUT / f"{name}_grayscale.png"
    Image.open(OUT / f"{name}.png").convert("L").save(grey, dpi=(300, 300))
    exports.append(str(grey))
    # Deterministic red/green-channel stress previews complement the grayscale
    # check. They are diagnostics only and never enter the manuscript.
    rgb = np.asarray(Image.open(OUT / f"{name}.png").convert("RGB"), dtype=float) / 255.0
    transforms = {
        "protanopia": np.array([[0.567, 0.433, 0.000], [0.558, 0.442, 0.000], [0.000, 0.242, 0.758]]),
        "deuteranopia": np.array([[0.625, 0.375, 0.000], [0.700, 0.300, 0.000], [0.000, 0.300, 0.700]]),
    }
    cvd_previews: dict[str, str] = {}
    for label, matrix in transforms.items():
        converted = np.clip(rgb @ matrix.T, 0, 1)
        path = OUT / f"{name}_{label}.png"
        Image.fromarray(np.uint8(np.round(converted * 255))).save(path, dpi=(300, 300))
        exports.append(str(path))
        cvd_previews[label] = str(path)
    minimum_font_pt = min(
        [text.get_fontsize() for text in fig.findobj(match=lambda obj: hasattr(obj, "get_fontsize"))]
        or [float("nan")]
    )
    audit = {
        "visual_qa_verdict": verdict,
        "visual_qa_issues": issues,
        "preview": preview,
        "exports": exports,
        "colour_vision_previews": cvd_previews,
        "minimum_text_size_pt": minimum_font_pt,
        "minimum_text_size_gate_pt": 5.0,
        "minimum_text_size_pass": bool(minimum_font_pt >= 5.0),
        "sources": [str(path) for path in sources],
        "source_sha256": {str(path): file_sha(path) for path in sources},
        "interpretation_boundary": boundary,
        "scientific_recomputation": False,
    }
    (OUT / f"{name}_build_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    plt.close(fig)
    if verdict == "FAIL" or minimum_font_pt < 5.0:
        raise RuntimeError(f"visual QA failed for {name}")
    return audit


def haversine_km(lat1: float, lon1: float, lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    radius = 6371.0088
    p1 = math.radians(lat1)
    p2 = np.radians(lat)
    dp = p2 - p1
    dl = np.radians(lon - lon1)
    a = np.sin(dp / 2) ** 2 + math.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * radius * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def strict_record_mask(values: np.ndarray) -> np.ndarray:
    mask = np.zeros(len(values), dtype=bool)
    running = -np.inf
    for index, value in enumerate(values):
        if np.isfinite(value) and value > running + 1e-10:
            mask[index] = True
            running = value
    return mask


def load_example() -> tuple[pd.DataFrame, pd.Series]:
    a_path = D104 / "family_a_event_metrics.csv.gz"
    b_path = D104 / "family_b_event_metrics.csv.gz"
    a = pd.read_csv(a_path).query("dataset == 'goto' and scale_m == 100")
    b = pd.read_csv(b_path).query("dataset == 'goto' and scale_m == 100")
    candidate = (
        a.drop(columns=["rho_tau_distinct_cell"], errors="ignore")
        .merge(b[["event_id", "C_union", "rho_tau_distinct_cell"]], on="event_id")
    )
    candidate = candidate[
        candidate["n_runs"].between(7, 25)
        & (candidate["n_record_runs"] >= 4)
        & candidate["rho_tau_distinct_cell"].fillna(False)
        & (candidate["L_union"] > 0)
        & (candidate["C_union"] > 0)
    ].copy()
    # Reuse the already frozen display event when available.  This is a
    # presentation cache, not an analysis input; it avoids rescanning the
    # dense interpolation table on every manuscript rebuild while retaining
    # the deterministic selection rule recorded in the source-data rows.
    display_cache = ROOT / "data/source_data/fig1_source_data.csv"
    if display_cache.exists():
        cached = pd.read_csv(display_cache)
        required = {"panel", "event_key", "is_endpoint", "is_trigger", "Time", "Lat", "Lon"}
        if required.issubset(cached.columns):
            event = cached[cached["panel"].eq("a_b")].copy()
            keys_in_cache = event["event_key"].dropna().astype(str).unique()
            if len(keys_in_cache) == 1 and keys_in_cache[0] in set(a["event_id"].astype(str)):
                selected_key = keys_in_cache[0]
                event["Time"] = pd.to_datetime(event["Time"], errors="coerce")
                event = event.dropna(subset=["Time", "Lat", "Lon"]).sort_values(
                    ["Time", "orig_idx_num"]
                ).drop_duplicates(subset=["Time", "Lat", "Lon"], keep="last")
                for column in ["is_step_start", "is_endpoint", "is_trigger",
                               "use_for_movement_step", "is_post_endpoint_to_trigger"]:
                    event[column] = event[column].astype(int)
                event["strict_record"] = event["strict_record"].astype(bool)
                summary = a.set_index("event_id").loc[selected_key]
                print(f"[example cache] reused {selected_key} ({len(event)} fixes)", flush=True)
                return event, summary
    canon = pd.read_csv(
        CANON,
        usecols=["event_key", "track_id", "segment_id", "start_orig_idx",
                 "endpoint_orig_idx", "trigger_orig_idx"],
    )
    candidate = candidate.merge(canon, left_on="event_id", right_on="event_key",
                                how="inner", validate="one_to_one")
    candidate["display_score"] = candidate["L_union"] + candidate["C_union"]
    candidate = candidate.sort_values(
        ["display_score", "n_runs", "event_id"], ascending=[False, True, True]
    ).head(60)
    keys = candidate["event_id"].tolist()
    wanted_pairs = set(zip(candidate["track_id"].astype(str), candidate["segment_id"].astype(int)))
    pieces = []
    usecols = [
        "track_id", "segment_id", "Time", "Lat", "Lon", "logCHL", "orig_idx",
        "edge_start_orig_idx", "edge_end_orig_idx", "interpolation_fraction_on_edge",
    ]
    started = time.monotonic()
    scanned = 0
    for chunk_index, chunk in enumerate(pd.read_csv(DENSE, usecols=usecols, chunksize=250_000), start=1):
        scanned += len(chunk)
        keep = [pair in wanted_pairs for pair in zip(chunk["track_id"].astype(str),
                                                      chunk["segment_id"].astype(int))]
        found = chunk.loc[keep]
        if len(found):
            pieces.append(found)
        if chunk_index == 1 or chunk_index % 10 == 0:
            print(
                f"[example scan] chunks={chunk_index}, rows={scanned:,}, "
                f"candidate rows={sum(len(x) for x in pieces):,}, "
                f"elapsed={time.monotonic()-started:.1f}s",
                flush=True,
            )
    if not pieces:
        raise RuntimeError("no dense rows found for deterministic Figure 1 candidates")
    dense = pd.concat(pieces, ignore_index=True)
    dense["position"] = np.where(
        dense["orig_idx"].notna(), dense["orig_idx"],
        dense["edge_start_orig_idx"] + dense["interpolation_fraction_on_edge"]
        * (dense["edge_end_orig_idx"] - dense["edge_start_orig_idx"]),
    )
    selected_key = None
    for key in keys:
        meta = candidate[candidate["event_id"].eq(key)].iloc[0]
        group = dense[
            dense["track_id"].astype(str).eq(str(meta["track_id"]))
            & dense["segment_id"].astype(int).eq(int(meta["segment_id"]))
            & dense["position"].between(float(meta["start_orig_idx"]),
                                         float(meta["trigger_orig_idx"]), inclusive="both")
        ]
        originals = group[group["orig_idx"].notna()]
        marks = set(originals["orig_idx"].astype(int))
        if (len(originals) >= 7 and int(meta["endpoint_orig_idx"]) in marks
                and int(meta["trigger_orig_idx"]) in marks):
            selected_key = key
            break
    if selected_key is None:
        raise RuntimeError("no Figure 1 candidate contained both rho and tau rows")
    meta = candidate[candidate["event_id"].eq(selected_key)].iloc[0]
    event = dense[
        dense["track_id"].astype(str).eq(str(meta["track_id"]))
        & dense["segment_id"].astype(int).eq(int(meta["segment_id"]))
        & dense["position"].between(float(meta["start_orig_idx"]),
                                     float(meta["trigger_orig_idx"]), inclusive="both")
        & dense["orig_idx"].notna()
    ].copy()
    event["event_key"] = selected_key
    event["orig_idx_num"] = event["orig_idx"]
    event["is_step_start"] = event["orig_idx"].eq(float(meta["start_orig_idx"])).astype(int)
    event["is_endpoint"] = event["orig_idx"].eq(float(meta["endpoint_orig_idx"])).astype(int)
    event["is_trigger"] = event["orig_idx"].eq(float(meta["trigger_orig_idx"])).astype(int)
    event["use_for_movement_step"] = event["orig_idx"].le(float(meta["endpoint_orig_idx"])).astype(int)
    event["is_post_endpoint_to_trigger"] = (
        event["orig_idx"].gt(float(meta["endpoint_orig_idx"]))
        & event["orig_idx"].le(float(meta["trigger_orig_idx"]))
    ).astype(int)
    event["Time"] = pd.to_datetime(event["Time"], errors="coerce")
    event = event.dropna(subset=["Time", "Lat", "Lon"]).sort_values(
        ["Time", "orig_idx_num"]
    ).drop_duplicates(subset=["Time", "Lat", "Lon"], keep="last")
    event["minutes"] = (event["Time"] - event["Time"].iloc[0]).dt.total_seconds() / 60
    event["radial_km"] = haversine_km(
        float(event["Lat"].iloc[0]), float(event["Lon"].iloc[0]),
        event["Lat"].to_numpy(float), event["Lon"].to_numpy(float),
    )
    event["strict_record"] = strict_record_mask(event["radial_km"].to_numpy())
    move = event["use_for_movement_step"].eq(1) & event["logCHL"].notna()
    event["path_chl_rank"] = np.nan
    event.loc[move, "path_chl_rank"] = event.loc[move, "logCHL"].rank(
        method="average", pct=True
    )
    summary = candidate.set_index("event_id").loc[selected_key]
    event["selection_rule"] = (
        "highest L_union+C_union among eligible 7-25-run events; first candidate "
        "with dense rho and tau support"
    )
    return event, summary


def errorbar_x(ax: plt.Axes, estimate: float, low: float, high: float,
               y: float, color: str, marker: str = "o", label: str | None = None) -> None:
    ax.errorbar(
        estimate, y, xerr=np.array([[estimate - low], [high - estimate]]),
        fmt=marker, color=color, ecolor=BLACK, elinewidth=0.75, capsize=2,
        markersize=4.2, markeredgecolor=BLACK, markeredgewidth=0.3, label=label,
    )


def build_figure1(event: pd.DataFrame, example: pd.Series) -> dict:
    summary_path = D104 / "formal_summary.csv"
    summary = pd.read_csv(summary_path)
    source = pd.concat([
        event.assign(panel="a_b", record_type="deterministic_observed_event"),
        summary.query("dataset == 'goto' and scale_m == 100").assign(
            panel="d_e", record_type="formal_effect"
        ),
    ], ignore_index=True, sort=False)
    source.to_csv(OUT / "fig1_source_data.csv", index=False)

    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    fig = plt.figure(figsize=(7.2, 6.7), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 0.72, 0.92])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, :])
    axd = fig.add_subplot(gs[2, 0])
    axe = fig.add_subplot(gs[2, 1])

    axa.plot(event["minutes"], event["radial_km"], color=OBSERVED, lw=1.2)
    records = event[event["strict_record"]]
    axa.scatter(records["minutes"], records["radial_km"], facecolor="white",
                edgecolor=RECORD, s=18, lw=0.65, zorder=4, label="Candidate record")
    rho = event[event["is_endpoint"].eq(1)].iloc[-1]
    tau = event[event["is_trigger"].eq(1)].iloc[-1]
    axa.scatter(rho["minutes"], rho["radial_km"], marker="*", s=68,
                color=YELLOW, edgecolor=BLACK, lw=0.45, zorder=5, label=r"Boundary $\rho$")
    axa.scatter(tau["minutes"], tau["radial_km"], marker="X", s=35,
                color=CONFIRMATION, edgecolor=BLACK, lw=0.35, zorder=5, label=r"Confirmation $\tau$")
    axa.set_xlabel("Time from event origin (min)")
    axa.set_ylabel("Radial distance (km)")
    axa.set_title("Candidate records", loc="left")
    axa.legend(frameon=False, fontsize=5.3, loc="best")

    lat0, lon0 = float(event["Lat"].iloc[0]), float(event["Lon"].iloc[0])
    event["x_km"] = (event["Lon"] - lon0) * 111.32 * math.cos(math.radians(lat0))
    event["y_km"] = (event["Lat"] - lat0) * 111.32
    rho = event[event["is_endpoint"].eq(1)].iloc[-1]
    tau = event[event["is_trigger"].eq(1)].iloc[-1]
    axb.plot(event["x_km"], event["y_km"], color=LIGHT, lw=0.85, zorder=1)
    visited = event[event["path_chl_rank"].notna()]
    scatter = axb.scatter(
        visited["x_km"], visited["y_km"], c=visited["path_chl_rank"],
        cmap=cmocean.cm.algae, vmin=0, vmax=1, s=18, edgecolor="none", zorder=3,
    )
    axb.scatter(rho["x_km"], rho["y_km"], marker="*", s=70,
                color=YELLOW, edgecolor=BLACK, lw=0.45, zorder=5)
    axb.scatter(tau["x_km"], tau["y_km"], marker="X", s=36,
                color=CONFIRMATION, edgecolor=BLACK, lw=0.35, zorder=5)
    axb.set_aspect("equal", adjustable="datalim")
    axb.set_xlabel("Local east (km)")
    axb.set_ylabel("Local north (km)")
    axb.set_title("Path-relative CHL", loc="left")
    colorbar = fig.colorbar(scatter, ax=axb, fraction=0.05, pad=0.02)
    colorbar.set_label("Within-path CHL rank")

    e_value = float(example["E_union"])
    r_value = float(example["R_union"])
    l_value = float(example["L_union"])
    xpos = np.arange(3)
    # Waterfall geometry makes the exact identity visible: R starts at zero,
    # L starts at R, and E is the resulting total.
    axc.bar(0, r_value, bottom=0, color=PALE_ROSE, width=0.56,
            edgecolor=ROSE, linewidth=0.8)
    axc.bar(1, l_value, bottom=r_value, color=PALE_BLUE, width=0.56,
            edgecolor=OBSERVED, linewidth=0.8)
    axc.bar(2, e_value, bottom=0, facecolor=PALE_VIOLET, width=0.56,
            edgecolor=PURPLE, linewidth=0.9)
    axc.plot([0.28, 0.72], [r_value, r_value], color=BLACK, lw=0.55)
    axc.plot([1.28, 1.72], [e_value, e_value], color=BLACK, lw=0.55)
    axc.axhline(0, color=BLACK, lw=0.65)
    axc.set_xticks(xpos, [r"Record $R$", r"Last record $L$", r"Endpoint $E$"])
    axc.set_ylabel("Combined-tail contrast")
    axc.set_title(r"Exact decomposition: $E=R+L$", loc="left")

    a_rows = summary.query(
        "dataset == 'goto' and scale_m == 100 and family == 'A' and "
        "estimand in ['R_high','R_low','R_union','L_high','L_low','L_union']"
    ).copy()
    tails = ["high", "low", "union"]
    y = np.arange(3)[::-1]
    for offset, component, color, marker in [(-0.11, "R", RECORD, "o"), (0.11, "L", OBSERVED, "D")]:
        for pos, tail in zip(y, tails):
            row = a_rows[a_rows["estimand"].eq(f"{component}_{tail}")].iloc[0]
            errorbar_x(axd, row["observed_unit_equal"], row["bootstrap_ci_low"],
                       row["bootstrap_ci_high"], pos + offset, color, marker,
                       component if tail == "high" else None)
            axd.scatter(row["phase_null_mean"], pos + offset, marker="x", s=18,
                        color=NULL, lw=0.8, zorder=4)
    axd.axvline(0, color=BLACK, lw=0.65, ls=(0, (3, 2)))
    axd.set_yticks(y, ["Upper tail", "Lower tail", "Either tail"])
    axd.set_xlabel("Observed contrast (95% CI)")
    axd.set_title("Selection component", loc="left")
    axd.scatter([], [], marker="x", s=18, color=NULL, lw=0.8, label="Phase-null mean")
    axd.legend(frameon=False, fontsize=5.3, loc="center right", ncol=1,
               handletextpad=0.5, labelspacing=0.3)

    b_rows = summary.query(
        "dataset == 'goto' and scale_m == 100 and family == 'B' and "
        "estimand in ['C_high','C_low','C_union']"
    ).copy()
    for pos, tail in zip(y, tails):
        row = b_rows[b_rows["estimand"].eq(f"C_{tail}")].iloc[0]
        errorbar_x(axe, row["observed_unit_equal"], row["bootstrap_ci_low"],
                   row["bootstrap_ci_high"], pos, OBSERVED, "o")
        axe.scatter(row["phase_null_mean"], pos, marker="x", s=18,
                    color=NULL, lw=0.8, zorder=4)
    axe.axvline(0, color=BLACK, lw=0.65, ls=(0, (3, 2)))
    axe.set_yticks(y, ["Upper tail", "Lower tail", "Either tail"])
    axe.set_xlabel("Boundary − confirmation contrast (95% CI)")
    axe.set_title("Boundary localisation", loc="left")

    for ax in [axa, axb, axc, axd, axe]:
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    return finish(
        fig, "fig1_last_passage_decomposition", (7.2, 6.7),
        [D104 / "family_a_event_metrics.csv.gz", D104 / "family_b_event_metrics.csv.gz",
         D104 / "formal_summary.csv", ROOT / "data/source_data/fig1_source_data.csv"],
        "rho is retrospective; path ranks are displayed only at visited locations; phase-null markers preserve path geometry and environmental autocorrelation.",
        [axa, axb, axc, axd, axe],
    )


def build_figure2() -> dict:
    d104_path = D104 / "formal_summary.csv"
    laysan_path = D127 / "laysan_phase_summary.csv"
    landing_path = D45 / "uesaka_lag_bin_summary.csv"
    shear_path = D107 / "behavior_group_effects.csv"
    d104 = pd.read_csv(d104_path)
    laysan = pd.read_csv(laysan_path)
    landing = pd.read_csv(landing_path).query("delta_m == 100").copy()
    shear = pd.read_csv(shear_path).query(
        "behavior_group == 'forage_dominant' and estimand == 'L_union'"
    ).copy()
    source = pd.concat([
        d104.query("family == 'A' and estimand == 'L_union'").assign(panel="a"),
        laysan.assign(panel="b"), landing.assign(panel="c"), shear.assign(panel="d"),
    ], ignore_index=True, sort=False)
    source.to_csv(OUT / "fig2_source_data.csv", index=False)

    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.7), constrained_layout=True)
    axa, axb, axc, axd = axs.flat

    repeated = d104.query("family == 'A' and estimand == 'L_union'").copy()
    repeated["label"] = repeated.apply(
        lambda row: "Wandering albatross\n100 m" if row["dataset"] == "goto"
        else f"Red-footed booby\n{int(row['scale_m']):,} m", axis=1
    )
    yy = np.arange(len(repeated))[::-1]
    colors = [OBSERVED for _ in repeated["dataset"]]
    markers = ["o" if value == "goto" else "D" for value in repeated["dataset"]]
    for yval, row, color, marker in zip(yy, repeated.itertuples(), colors, markers):
        errorbar_x(axa, row.observed_unit_equal, row.bootstrap_ci_low,
                   row.bootstrap_ci_high, yval, color, marker)
        axa.scatter(row.phase_null_mean, yval, marker="x", s=18,
                    color=NULL, lw=0.8, zorder=4)
    axa.axvline(0, color=BLACK, lw=0.65, ls=(0, (3, 2)))
    axa.set_yticks(yy, repeated["label"])
    axa.set_xlabel(r"$L_{union}$ (95% CI)")
    axa.set_title("Primary and booby systems", loc="left")
    axa.text(0.03, 0.04, "grey × = phase-null mean", transform=axa.transAxes, fontsize=5.3)

    laysan = laysan.query("estimand == 'L_union'").sort_values("scale_m")
    laysan_y = np.arange(len(laysan))[::-1]
    for yval, row in zip(laysan_y, laysan.itertuples()):
        errorbar_x(axb, row.observed_unit_equal, row.bootstrap_ci_low,
                   row.bootstrap_ci_high, yval, OBSERVED, "o")
        axb.scatter(row.phase_null_mean, yval, marker="x",
                    color=NULL, s=20, lw=0.8)
    axb.set_yticks(laysan_y, [f"{int(v):,}" for v in laysan["scale_m"]])
    axb.set_xlabel(r"$L_{union}$ (95% CI)")
    axb.set_ylabel("Drawdown scale (m)")
    axb.set_title("Laysan albatross", loc="left")
    axb.text(0.03, 0.05, "grey × = phase-null mean", transform=axb.transAxes, fontsize=5.3)

    landing["window"] = landing.apply(
        lambda row: f"{int(row.bin_low_s)}–{int(row.bin_high_s)} s", axis=1
    )
    x = np.arange(len(landing))
    axc.errorbar(
        x, landing["post_minus_pre"],
        yerr=np.vstack([landing["post_minus_pre"] - landing["ci_low"],
                        landing["ci_high"] - landing["post_minus_pre"]]),
        fmt="o", color=OBSERVED, ecolor=BLACK, elinewidth=0.75, capsize=2,
        markersize=4.2, markeredgecolor=BLACK, markeredgewidth=0.3,
    )
    axc.axhline(0, color=BLACK, lw=0.65, ls=(0, (3, 2)))
    axc.set_xticks(x, landing["window"], rotation=25, ha="right")
    axc.set_ylabel("Landing probability\n(post − pre)")
    axc.set_xlabel("Matched interval from 100-m boundary")
    axc.set_title("Landing timing", loc="left")

    shear = shear.sort_values("scale_m")
    shear_x = np.arange(len(shear))
    axd.scatter(shear_x, shear["phase_null_mean"], color=NULL,
                marker="x", s=20, linewidths=0.8, label="Phase-null mean")
    axd.errorbar(
        shear_x, shear["individual_equal_mean"],
        yerr=np.vstack([shear["individual_equal_mean"] - shear["bootstrap_ci_low"],
                        shear["bootstrap_ci_high"] - shear["individual_equal_mean"]]),
        fmt="D", color=OBSERVED, ecolor=BLACK, elinewidth=0.75, capsize=2,
        markersize=4.0, markeredgecolor=BLACK, markeredgewidth=0.3,
        label="Foraging-dominant events",
    )
    axd.axhline(0, color=BLACK, lw=0.65)
    axd.set_xticks(shear_x, [f"{int(v):,}" for v in shear["scale_m"]])
    axd.set_xlabel("Drawdown scale (m)")
    axd.set_ylabel(r"$L_{union}$ (95% CI)")
    axd.set_title("Foraging context", loc="left")
    axd.legend(frameon=False, fontsize=5.4, loc="best")

    for ax in axs.flat:
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    return finish(
        fig, "fig2_cross_system_context", (7.2, 5.7),
        [d104_path, laysan_path, landing_path, shear_path],
        "Cross-system union replication and behavioural context are shown; panels c-d do not establish a CHL-specific action, prey capture or sensory stimulus.",
        [axa, axb, axc, axd],
    )


def build_figure3() -> dict:
    model_path = D127 / "joint_model_summary.csv"
    grid_path = D127 / "conditional_L_3x3_grid.csv"
    models = pd.read_csv(model_path)
    grids = pd.read_csv(grid_path)
    source_rows: list[dict[str, object]] = []
    source_specs = [
        ("absolute CHL", "beta_absolute", "beta_absolute_ci_low", "beta_absolute_ci_high"),
        ("L_low", "beta_L_low", "beta_L_low_ci_low", "beta_L_low_ci_high"),
        ("L_high", "beta_L_high", "beta_L_high_ci_low", "beta_L_high_ci_high"),
    ]
    model_lookup = models.set_index(["dataset", "scale_m"])
    for row in models.to_dict("records"):
        common = {
            "figure": "Main Figure 3", "dataset": row["dataset"],
            "scale_m": row["scale_m"], "n_events": row["events"],
            "n_units": row["animals"], "bootstrap_reps": row["bootstrap_reps"],
        }
        for series, value, low, high in source_specs:
            source_rows.append({
                **common, "panel": "a", "record_type": "joint model coefficient",
                "series": series, "x_category": "", "y_category": "",
                "estimate": row[value], "ci_low": row[low], "ci_high": row[high],
            })
        source_rows.append({
            **common, "panel": "b", "record_type": "cell-standardised mean",
            "series": "L_union", "x_category": "", "y_category": "",
            "estimate": row["standardized_L_union"],
            "ci_low": row["standardized_L_union_ci_low"],
            "ci_high": row["standardized_L_union_ci_high"],
            "positive_cells": row["positive_grids"],
            "estimable_cells": row["estimable_grids"],
        })
    panel_map = {
        ("goto", 100): "c", ("usgs_laysan_albatross", 500): "d",
        ("usgs_laysan_albatross", 1000): "e",
        ("usgs_laysan_albatross", 2000): "f",
    }
    for row in grids.to_dict("records"):
        model = model_lookup.loc[(row["dataset"], row["scale_m"])]
        source_rows.append({
            "figure": "Main Figure 3", "panel": panel_map[(row["dataset"], row["scale_m"])],
            "record_type": "conditional tertile cell", "dataset": row["dataset"],
            "scale_m": row["scale_m"], "series": "L_union",
            "x_category": ["low", "mid", "high"][int(row["abs_tertile"])],
            "y_category": ["short", "mid", "long"][int(row["length_tertile"])],
            "estimate": row["L_union"], "ci_low": np.nan, "ci_high": np.nan,
            "n_events": model["events"], "n_units": model["animals"],
            "bootstrap_reps": model["bootstrap_reps"],
            "positive_cells": model["positive_grids"],
            "estimable_cells": model["estimable_grids"],
        })
    pd.DataFrame(source_rows).to_csv(OUT / "fig3_source_data.csv", index=False)

    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 4, height_ratios=[1.18, 0.82])
    axa = fig.add_subplot(gs[0, :3])
    axb = fig.add_subplot(gs[0, 3])
    heat_axes = [fig.add_subplot(gs[1, col]) for col in range(4)]

    group_labels = []
    for row in models.itertuples():
        group_labels.append(
            "Wandering\n100 m" if row.dataset == "goto"
            else f"Laysan\n{int(row.scale_m):,} m"
        )
    ybase = np.arange(len(models))[::-1]
    specs = [
        ("beta_absolute", "beta_absolute_ci_low", "beta_absolute_ci_high", RECORD, "o", "Absolute CHL"),
        ("beta_L_low", "beta_L_low_ci_low", "beta_L_low_ci_high", BLUE, "v", r"$L_{low}$"),
        ("beta_L_high", "beta_L_high_ci_low", "beta_L_high_ci_high", ORANGE, "^", r"$L_{high}$"),
    ]
    for offset, (value, low, high, color, marker, label) in zip([-0.18, 0, 0.18], specs):
        for pos, row in zip(ybase, models.itertuples()):
            errorbar_x(axa, getattr(row, value), getattr(row, low), getattr(row, high),
                       pos + offset, color, marker, label if pos == ybase[0] else None)
    axa.axvline(0, color=BLACK, lw=0.65, ls=(0, (3, 2)))
    axa.set_yticks(ybase, group_labels)
    axa.set_xlabel("Standardised coefficient (95% CI)")
    axa.set_title("Bout-length model", loc="left")
    axa.set_ylim(-0.5, 4.0)
    axa.legend(frameon=False, fontsize=5.5, ncol=3, loc="upper center",
               bbox_to_anchor=(0.55, 0.995))

    for pos, row in zip(ybase, models.itertuples()):
        errorbar_x(axb, row.standardized_L_union, row.standardized_L_union_ci_low,
                   row.standardized_L_union_ci_high, pos, PURPLE, "s")
    axb.axvline(0, color=BLACK, lw=0.65, ls=(0, (3, 2)))
    axb.set_yticks(ybase, ["" for _ in ybase])
    axb.set_xlabel(r"Standardised $L_{union}$")
    axb.set_title("Conditional mean", loc="left")

    vmin, vmax = float(grids["L_union"].min()), float(grids["L_union"].max())
    last_image = None
    for ax, ((dataset, scale), frame) in zip(
        heat_axes, grids.groupby(["dataset", "scale_m"], sort=False)
    ):
        matrix = frame.pivot(index="length_tertile", columns="abs_tertile", values="L_union").sort_index(ascending=False)
        last_image = ax.imshow(matrix.to_numpy(), cmap=cmocean.cm.algae,
                               vmin=vmin, vmax=vmax, aspect="equal")
        title = "Wandering 100 m" if dataset == "goto" else f"Laysan {int(scale):,} m"
        ax.set_title(title, fontsize=6.2)
        ax.set_xticks([0, 1, 2], ["Low", "Mid", "High"])
        ax.set_yticks([0, 1, 2], ["Long", "Mid", "Short"])
        ax.set_xlabel("Absolute CHL tertile")
        if ax is heat_axes[0]:
            ax.set_ylabel("Bout-length tertile")
        else:
            ax.set_ylabel("")
    colorbar = fig.colorbar(last_image, ax=heat_axes, fraction=0.025, pad=0.015)
    colorbar.set_label(r"Mean $L_{union}$")

    all_axes = [axa, axb] + heat_axes
    for ax in all_axes:
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    return finish(
        fig, "fig3_dual_reference", (7.2, 5.4), [model_path, grid_path],
        "Absolute and path-relative quantities are same-event retrospective statistics, not demonstrated sensory or cognitive reference frames.",
        all_axes,
    )


def build_figure4(event: pd.DataFrame) -> dict:
    null_path = D128 / "eligible_record_counterfactual_distribution.csv.gz"
    final_path = D128 / "final_summary.json"
    null = pd.read_csv(null_path)
    final = json.loads(final_path.read_text(encoding="utf-8"))
    observed_support = final["rho"]["rho_chl_support"]
    observed_length = final["rho"]["rho_chl_median_log_length"]
    pd.concat([
        event.assign(panel="a", record_type="candidate_boundary_example"),
        null.assign(panel="b_d", record_type="eligible_record_counterfactual"),
        pd.DataFrame([{
            "panel": "b_d", "record_type": "observed_boundary",
            "rho_chl_support": observed_support,
            "rho_chl_median_log_length": observed_length,
        }]),
    ], ignore_index=True, sort=False).to_csv(OUT / "fig4_source_data.csv", index=False)

    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.5), constrained_layout=True)
    axa, axb, axc, axd = axs.flat

    records = event[event["strict_record"] & event["use_for_movement_step"].eq(1)].copy()
    rho = event[event["is_endpoint"].eq(1)].iloc[-1]
    eligible = records[records["minutes"] < rho["minutes"]]
    alternative = eligible.iloc[len(eligible) // 2] if len(eligible) else records.iloc[0]
    axa.plot(event["minutes"], event["radial_km"], color=LIGHT, lw=1.0)
    axa.scatter(records["minutes"], records["radial_km"], facecolor="white",
                edgecolor=RECORD, s=20, lw=0.65, label="Eligible record")
    axa.scatter(rho["minutes"], rho["radial_km"], marker="*", s=75,
                color=YELLOW, edgecolor=BLACK, lw=0.4, zorder=5, label="Observed boundary")
    axa.scatter(alternative["minutes"], alternative["radial_km"], marker="D", s=34,
                facecolor="white", edgecolor=RECORD, lw=0.8, zorder=5,
                label="Example replacement")
    axa.set_xlabel("Time from event origin (min)")
    axa.set_ylabel("Radial distance (km)")
    axa.set_title("Boundary replacement", loc="left")
    axa.legend(frameon=False, fontsize=5.2, loc="best")

    axb.hist(null["rho_chl_support"], bins=28, color=PALE_ROSE,
             edgecolor="white", linewidth=0.4)
    axb.axvline(observed_support, color=OBSERVED, lw=1.7, label=f"Observed = {observed_support:.3f}")
    axb.axvline(null["rho_chl_support"].median(), color=BLACK, lw=0.8,
                ls=(0, (3, 2)), label=f"Null median = {null['rho_chl_support'].median():.3f}")
    axb.set_xlabel("CHL–relative Lomax-support correlation")
    axb.set_ylabel("Counterfactual partitions")
    axb.set_title("Tail-shape support", loc="left")
    axb.legend(frameon=False, fontsize=5.3, loc="best")
    axb.text(0.97, 0.07, "one-sided P = 0.004", transform=axb.transAxes,
             ha="right", va="bottom", fontsize=5.7)

    axc.hist(null["rho_chl_median_log_length"], bins=28, color=PALE_ROSE,
             edgecolor="white", linewidth=0.4)
    axc.axvline(observed_length, color=OBSERVED, lw=1.7, label=f"Observed = {observed_length:.3f}")
    axc.axvline(null["rho_chl_median_log_length"].median(), color=BLACK,
                lw=0.8, ls=(0, (3, 2)), label=f"Null median = {null['rho_chl_median_log_length'].median():.3f}")
    axc.set_xlabel("CHL–median-log-length correlation")
    axc.set_ylabel("Counterfactual partitions")
    axc.set_title("Median length", loc="left")
    axc.legend(frameon=False, fontsize=5.3, loc="best")
    axc.text(0.97, 0.07, "one-sided P = 0.999", transform=axc.transAxes,
             ha="right", va="bottom", fontsize=5.7)

    axd.scatter(null["rho_chl_support"], null["rho_chl_median_log_length"],
                s=8, alpha=0.58, color=PALE_ROSE, edgecolor="none",
                label="999 replacements")
    axd.scatter(observed_support, observed_length, marker="D", s=46,
                color=OBSERVED, edgecolor=BLACK, lw=0.4, zorder=5, label="Observed boundaries")
    axd.axvline(null["rho_chl_support"].median(), color=BLACK, lw=0.55, ls=(0, (3, 2)))
    axd.axhline(null["rho_chl_median_log_length"].median(), color=BLACK, lw=0.55, ls=(0, (3, 2)))
    axd.set_xlabel("CHL–tail-support correlation")
    axd.set_ylabel("CHL–median-length correlation")
    axd.set_title("Joint contrast", loc="left")
    axd.legend(frameon=False, fontsize=5.3, loc="best")

    for ax in axs.flat:
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    return finish(
        fig, "fig4_boundary_counterfactual", (7.2, 5.5),
        [null_path, final_path, ROOT / "data/source_data/fig1_source_data.csv"],
        "The eligible-record distribution is a statistical boundary-placement counterfactual with fixed event count; it is not an animal intervention and contains no renewal-density contrast.",
        [axa, axb, axc, axd],
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    overall = time.monotonic()
    print("[0/5 0%] selecting a deterministic observed event", flush=True)
    event, example = load_example()
    print(f"[1/5 20%] event selected; elapsed={time.monotonic()-overall:.1f}s", flush=True)
    audits = {}
    builders = [
        ("figure1", lambda: build_figure1(event.copy(), example)),
        ("figure2", build_figure2),
        ("figure3", build_figure3),
        ("figure4", lambda: build_figure4(event.copy())),
    ]
    for index, (name, builder) in enumerate(builders, start=2):
        print(f"[{index-1}/5 {(index-1)*20}%] building {name}", flush=True)
        audits[name] = builder()
        print(
            f"[{index}/5 {index*20}%] {name} complete; "
            f"elapsed={time.monotonic()-overall:.1f}s",
            flush=True,
        )
    (OUT / "MAIN_FIGURE_PACKAGE_AUDIT.json").write_text(
        json.dumps(audits, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"[5/5 100%] main figure package complete; elapsed={time.monotonic()-overall:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
