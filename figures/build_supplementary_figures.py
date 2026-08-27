#!/usr/bin/env python3
"""Build Paper 2 V3 supplementary figures from frozen result tables only.

All calculations in this file are display transformations or algebraic audit
readbacks. Segmentation, environmental extraction, bootstrap, randomisation,
behaviour classification and model fitting are not rerun.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cmocean

from build_main_figures import (
    BASE,
    BLACK,
    BLUE,
    OBSERVED,
    NULL,
    RECORD,
    D104,
    D107,
    D127,
    D128,
    D45,
    GREEN,
    GREY,
    LIGHT,
    ORANGE,
    OUT,
    PURPLE,
    finish,
    setup_style,
)

D106 = BASE / "results/booby_behavior_context"
D126 = BASE / "results/booby_dive_timing"


def save_source(frame: pd.DataFrame, name: str) -> Path:
    path = OUT / f"{name}_source_data.csv"
    frame.to_csv(path, index=False)
    return path


def horizontal_ci(ax: plt.Axes, row: pd.Series, y: float, color: str,
                  marker: str, label: str | None = None) -> None:
    estimate = float(row["observed_unit_equal"])
    low = float(row["bootstrap_ci_low"])
    high = float(row["bootstrap_ci_high"])
    ax.errorbar(
        estimate, y,
        xerr=np.array([[estimate - low], [high - estimate]]),
        fmt=marker, ms=3.7, color=color, markeredgecolor=BLACK,
        markeredgewidth=0.25, ecolor=BLACK, elinewidth=0.65, capsize=1.8,
        label=label,
    )
    ax.scatter(float(row["phase_null_mean"]), y, marker="x", s=18,
               color=NULL, linewidths=0.8, zorder=4)


def build_s1() -> dict:
    setup_style("nature", "en", use_sciplots=False, constrained_layout=False)
    a_path = D104 / "family_a_event_metrics.csv.gz"
    a = pd.read_csv(a_path)
    rows: list[dict[str, object]] = []
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=False)
    axa, axb = axs.flat
    labels = {"goto": "Wandering albatross, 100 m", "rfbo": "Red-footed booby, four scales"}
    colors = {"goto": OBSERVED, "rfbo": RECORD}
    linestyles = {"goto": "-", "rfbo": "--"}
    for dataset, group in a.groupby("dataset"):
        for panel, column, ax in [("a", "n_runs", axa), ("b", "n_record_runs", axb)]:
            values = np.sort(group[column].dropna().to_numpy(int))
            xvals = np.arange(int(values.min()), int(values.max()) + 1)
            fractions = np.array([(values >= value).mean() for value in xvals])
            ax.step(xvals, fractions, where="post", linewidth=1.15,
                    color=colors[dataset], linestyle=linestyles[dataset],
                    label=f"{labels[dataset]} (n={len(values):,})")
            rows.extend({
                "figure": "Supplementary Figure 1", "panel": panel,
                "system": labels[dataset], "metric": column,
                "x_value": int(xval), "fraction_events_ge": float(frac),
                "event_scale_rows": int(len(values)),
            } for xval, frac in zip(xvals, fractions))
    source_path = save_source(pd.DataFrame(rows), "figS1_observation_support")
    axa.axvline(3, color=NULL, lw=0.7, ls=(0, (3, 2)))
    axb.axvline(2, color=NULL, lw=0.7, ls=(0, (3, 2)))
    axa.set(xlabel="CHL runs", ylabel="Fraction of events ≥ value",
            title="Environmental observations")
    axb.set(xlabel="Record runs", ylabel="Fraction of events ≥ value",
            title="Candidate records")
    axa.legend(frameon=False, fontsize=5.2)
    for ax in axs.flat:
        ax.set_xscale("log")
        ax.set_ylim(0, 1.02)
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    fig.set_layout_engine(None)
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.19, top=0.87, wspace=0.18)
    return finish(
        fig, "figS1_observation_support", (7.2, 2.8),
        [a_path, source_path],
        "Curves describe qualified event-scale rows; the four booby scales are pooled for observation-support display and are not treated as independent animals.",
        list(axs.flat),
    )


def build_s2() -> dict:
    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    summary_path = D104 / "formal_summary.csv"
    frame = pd.read_csv(summary_path)
    source_path = save_source(frame, "figS2_complete_last_passage")
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.8), sharex="row")
    tails = ["high", "low", "union"]
    titles = ["Upper tail", "Lower tail", "Either tail"]
    systems = [
        ("goto", 100, "Wandering 100 m"),
        ("rfbo", 250, "Booby 250 m"),
        ("rfbo", 500, "Booby 500 m"),
        ("rfbo", 1000, "Booby 1,000 m"),
        ("rfbo", 2000, "Booby 2,000 m"),
    ]
    y_base = np.arange(len(systems))[::-1]
    for col, (tail, title) in enumerate(zip(tails, titles)):
        top = axs[0, col]
        bottom = axs[1, col]
        for yi, (dataset, scale, _) in zip(y_base, systems):
            for offset, component, color, marker in [
                (0.13, "R", RECORD, "o"), (-0.13, "L", OBSERVED, "D")
            ]:
                row = frame.query(
                    "dataset == @dataset and scale_m == @scale and family == 'A' and estimand == @component + '_' + @tail"
                ).iloc[0]
                horizontal_ci(top, row, yi + offset, color, marker,
                              label=component if yi == y_base[0] else None)
            crow = frame.query(
                "dataset == @dataset and scale_m == @scale and family == 'B' and estimand == 'C_' + @tail"
            ).iloc[0]
            horizontal_ci(bottom, crow, yi, OBSERVED, "o",
                          label=r"$\rho-\tau$" if yi == y_base[0] else None)
        for ax in (top, bottom):
            ax.axvline(0, color=BLACK, lw=0.6, ls=(0, (3, 2)))
            ax.set_yticks(y_base)
            if col == 0:
                ax.set_yticklabels([s[2] for s in systems])
            else:
                ax.set_yticklabels([])
            ax.grid(False)
        top.set_title(title, loc="left")
        top.set_xlabel("Observed component (95% CI)")
        bottom.set_xlabel("Contrast (95% CI)")
    axs[0, 0].legend(frameon=False, fontsize=5.2, loc="lower right")
    axs[1, 0].legend(frameon=False, fontsize=5.2, loc="lower right")
    fig.suptitle("Last-record selection and boundary localisation", x=0.06,
                 ha="left", fontsize=8.0)
    return finish(
        fig, "figS2_complete_last_passage", (7.2, 5.8),
        [summary_path, source_path],
        "Grey crosses are phase-null means; confidence intervals describe biological-unit uncertainty around observed values and are not phase-net confidence intervals.",
        list(axs.flat),
    )


def build_s3() -> dict:
    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    phase_path = D127 / "laysan_phase_summary.csv"
    grid_path = D127 / "conditional_L_3x3_grid.csv"
    model_path = D127 / "joint_model_summary.csv"
    phase = pd.read_csv(phase_path).query("estimand == 'L_union'").sort_values("scale_m")
    grid = pd.read_csv(grid_path).query("dataset == 'usgs_laysan_albatross'")
    models = pd.read_csv(model_path).query("dataset == 'usgs_laysan_albatross'").sort_values("scale_m")
    source_path = save_source(pd.concat([
        phase.assign(panel="a"), grid.assign(panel="b-d"), models.assign(panel="summary")
    ], ignore_index=True, sort=False), "figS3_laysan_conditional")
    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.25, 1, 1, 1])
    axa = fig.add_subplot(gs[0, 0])
    axes = [axa] + [fig.add_subplot(gs[0, i]) for i in range(1, 4)]
    y = np.arange(len(phase))[::-1]
    for yi, (_, row) in zip(y, phase.iterrows()):
        horizontal_ci(axa, row, yi, OBSERVED, "o")
    axa.axvline(0, color=BLACK, lw=0.6, ls=(0, (3, 2)))
    axa.set_yticks(y, [f"{int(v):,} m" for v in phase["scale_m"]])
    axa.set_xlabel("Effect")
    axa.set_title(r"Laysan $L_{union}$ replication", loc="left")
    for ax, scale in zip(axes[1:], [500, 1000, 2000]):
        sub = grid.query("scale_m == @scale")
        matrix = sub.pivot(index="length_tertile", columns="abs_tertile", values="L_union").sort_index(ascending=False)
        image = ax.imshow(matrix, cmap=cmocean.cm.algae, vmin=0.02, vmax=0.29,
                          aspect="equal")
        ax.set_title(f"{scale:,} m · 9/9 positive", fontsize=6.4)
        ax.set_xticks(range(3), ["Low", "Mid", "High"])
        ax.set_yticks(range(3), ["Long", "Mid", "Short"] if ax is axes[1] else [])
        ax.set_xlabel("Absolute CHL tertile")
        if ax is axes[1]:
            ax.set_ylabel("Bout-length tertile")
    cbar = fig.colorbar(image, ax=axes[1:], fraction=0.028, pad=0.025)
    cbar.set_label(r"Mean $L_{union}$")
    for ax in axes:
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    return finish(
        fig, "figS3_laysan_conditional", (7.2, 3.6),
        [phase_path, grid_path, model_path, source_path],
        "The three scales reuse the same 34 Laysan albatrosses and therefore represent within-system robustness rather than independent biological replications.",
        axes,
    )


def build_s4() -> dict:
    setup_style("nature", "en", use_sciplots=False, constrained_layout=False)
    landing_path = D45 / "uesaka_scale_direction_summary.csv"
    slope_path = D107 / "primary_slope_gate.csv"
    forage_path = D107 / "behavior_group_effects.csv"
    dive_path = D106 / "behavior_group_effects.csv"
    modifier_path = D106 / "behavior_modifiers.csv"
    direction_path = D126 / "direction_contrasts.csv"
    landing = pd.read_csv(landing_path).sort_values("delta_m")
    slope = pd.read_csv(slope_path).sort_values("scale_m")
    forage = pd.read_csv(forage_path).query(
        "behavior_group == 'forage_dominant' and estimand == 'L_union'"
    ).sort_values("scale_m")
    dive = pd.read_csv(dive_path).query(
        "behavior_class == 'dive_near' and estimand == 'L_union'"
    ).sort_values("scale_m")
    modifier = pd.read_csv(modifier_path).query(
        "contrast == 'dive_minus_wet' and estimand == 'L_union'"
    ).sort_values("scale_m")
    direction = pd.read_csv(direction_path).query(
        "contrast == 'post_minus_pre' and estimand == 'L_union'"
    ).sort_values("scale_m")
    source_rows: list[dict[str, object]] = []

    def append_rows(frame: pd.DataFrame, panel: str, system: str, estimand: str,
                    scale_col: str, estimate_col: str, low_col: str, high_col: str,
                    null_col: str | None, events_col: str | None,
                    units_col: str | None, p_col: str | None) -> None:
        for row in frame.to_dict("records"):
            source_rows.append({
                "figure": "Supplementary Figure 4", "panel": panel,
                "system": system, "scale_m": row[scale_col], "estimand": estimand,
                "estimate": row[estimate_col], "ci_low": row[low_col],
                "ci_high": row[high_col],
                "phase_null_mean": row.get(null_col) if null_col else np.nan,
                "n_events": row.get(events_col) if events_col else np.nan,
                "n_units": row.get(units_col) if units_col else np.nan,
                "holm_adjusted_p": row.get(p_col) if p_col else np.nan,
            })

    append_rows(landing, "a", "Wandering albatross", "landing post minus pre",
                "delta_m", "actual_diff", "raw_diff_ci_low", "raw_diff_ci_high",
                "null_diff_mean", "events", "birds", "holm15_direction_p")
    append_rows(forage, "b", "Short-tailed shearwater", "foraging-dominant footprint",
                "scale_m", "individual_equal_mean", "bootstrap_ci_low", "bootstrap_ci_high",
                "phase_null_mean", "events", "individuals", "phase_p_holm")
    append_rows(slope, "b", "Short-tailed shearwater", "foraging-minus-rest slope",
                "scale_m", "individual_equal_slope", "bootstrap_ci_low", "bootstrap_ci_high",
                "phase_null_mean", "events", "individuals", "phase_p_holm")
    append_rows(dive, "c", "Red-footed booby", "dive-near footprint",
                "scale_m", "observed_unit_equal", "bootstrap_ci_low", "bootstrap_ci_high",
                "phase_null_mean", "events", "birds", "holm_phase_p")
    append_rows(modifier, "c", "Red-footed booby", "dive minus wet-only",
                "scale_m", "observed_unit_equal", "bootstrap_ci_low", "bootstrap_ci_high",
                "phase_null_mean", None, "paired_birds", "holm_phase_p")
    append_rows(direction, "d", "Red-footed booby", "post-dive minus pre-dive",
                "scale_m", "observed_unit_equal", "bootstrap_ci_low", "bootstrap_ci_high",
                "phase_null_mean", None, "paired_birds", "holm_phase_p")
    source_path = save_source(pd.DataFrame(source_rows), "figS4_behaviour_context")

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4), constrained_layout=False)
    axa, axb, axc, axd = axs.flat

    def categorical_ci(ax: plt.Axes, frame: pd.DataFrame, estimate: str,
                       low: str, high: str, positions: np.ndarray, color: str,
                       marker: str, label: str | None = None) -> None:
        values = frame[estimate].to_numpy(float)
        ax.errorbar(positions, values,
                    yerr=np.vstack([values - frame[low].to_numpy(float),
                                    frame[high].to_numpy(float) - values]),
                    fmt=marker, color=color, ecolor=BLACK, elinewidth=0.55,
                    capsize=1.5, ms=3.5, markeredgecolor=BLACK,
                    markeredgewidth=0.25, label=label)

    pos = np.arange(len(landing))
    categorical_ci(axa, landing, "actual_diff", "raw_diff_ci_low", "raw_diff_ci_high",
                   pos, OBSERVED, "o")
    show = {0, 2, 6, 9, 12, 14}
    axa.set_xticks(pos, [f"{int(v):,}" if i in show else "" for i, v in enumerate(landing["delta_m"])])
    axa.set(ylabel=r"$\Delta$ landing" + "\nprobability", title="Landing timing")

    pos = np.arange(len(forage))
    categorical_ci(axb, forage, "individual_equal_mean", "bootstrap_ci_low", "bootstrap_ci_high",
                   pos - 0.10, OBSERVED, "D", "Foraging-dominant footprint")
    categorical_ci(axb, slope, "individual_equal_slope", "bootstrap_ci_low", "bootstrap_ci_high",
                   pos + 0.10, RECORD, "o", "Foraging-minus-rest slope")
    axb.set_xticks(pos, [f"{int(v):,}" for v in forage["scale_m"]])
    axb.set(ylabel="Effect", title="Foraging context")
    axb.legend(frameon=False, fontsize=5.0)

    pos = np.arange(len(dive))
    categorical_ci(axc, dive, "observed_unit_equal", "bootstrap_ci_low", "bootstrap_ci_high",
                   pos - 0.10, OBSERVED, "D", "Dive-near footprint")
    categorical_ci(axc, modifier, "observed_unit_equal", "bootstrap_ci_low", "bootstrap_ci_high",
                   pos + 0.10, RECORD, "o", "Dive minus wet-only")
    axc.set_xticks(pos, [f"{int(v):,}" for v in dive["scale_m"]])
    axc.set(ylabel="Effect", title="Dive context")
    axc.legend(frameon=False, fontsize=5.0)

    pos = np.arange(len(direction))
    categorical_ci(axd, direction, "observed_unit_equal", "bootstrap_ci_low", "bootstrap_ci_high",
                   pos, OBSERVED, "s")
    axd.set_xticks(pos, [f"{int(v):,}" for v in direction["scale_m"]])
    axd.set(ylabel="Post-dive − pre-dive", title="Dive timing")

    for ax in axs.flat:
        ax.axhline(0, color=BLACK, lw=0.6, ls=(0, (3, 2)))
        ax.set_xlabel("Drawdown scale (m)")
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    fig.set_layout_engine(None)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.09, top=0.94,
                        wspace=0.22, hspace=0.28)
    return finish(
        fig, "figS4_behaviour_context", (7.2, 5.4),
        [landing_path, slope_path, forage_path, dive_path, modifier_path,
         direction_path, source_path],
        "Positive within-context footprints are not equivalent to significant between-context specificity; dives indicate attempts or underwater activity, not confirmed capture.",
        list(axs.flat),
    )


def build_s5() -> dict:
    setup_style("nature", "en", use_sciplots=False, constrained_layout=True)
    audit_path = D128 / "segment_boundary_audit.csv"
    null_path = D128 / "eligible_record_counterfactual_distribution.csv.gz"
    final_path = D128 / "final_summary.json"
    audit = pd.read_csv(audit_path)
    null = pd.read_csv(null_path)
    final = json.loads(final_path.read_text(encoding="utf-8"))
    observed_support = float(final["rho"]["rho_chl_support"])
    observed_length = float(final["rho"]["rho_chl_median_log_length"])
    diagnostic = pd.DataFrame({
        "property": ["Valid observed partition", "Non-monotone confirmation", "Repeated confirmation"],
        "segments": [int(audit["event_count_match"].sum()),
                     int((audit["tau_nonincreasing_adjacent_pairs"] > 0).sum()),
                     int((audit["tau_duplicate_fixes"] > 0).sum())],
    })
    diagnostic["total_segments"] = len(audit)
    diagnostic["percent_segments"] = 100 * diagnostic["segments"] / len(audit)
    source_rows = [
        *diagnostic.assign(figure="Supplementary Figure 5", panel="a",
                           record_type="partition diagnostic").to_dict("records"),
        *null.assign(figure="Supplementary Figure 5", panel="b-d",
                     record_type="eligible-record replacement").to_dict("records"),
        {"figure": "Supplementary Figure 5", "panel": "b-d",
         "record_type": "observed boundaries", "replicate": "observed",
         "rho_chl_support": observed_support,
         "rho_chl_median_log_length": observed_length},
    ]
    source_path = save_source(pd.DataFrame(source_rows), "figS5_boundary_counterfactual")

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4))
    axa, axb, axc, axd = axs.flat
    ypos = np.arange(len(diagnostic))[::-1]
    axa.hlines(ypos, 0, diagnostic["percent_segments"], color=LIGHT, lw=2.2)
    axa.scatter(diagnostic["percent_segments"], ypos, color=[OBSERVED, RECORD, RECORD],
                marker="o", s=26, edgecolor=BLACK, linewidth=0.3, zorder=3)
    axa.set_yticks(ypos, diagnostic["property"])
    axa.set_xlim(0, 105)
    axa.set_xlabel("Segments (%)")
    axa.set_title("Partition diagnostic", loc="left")
    for yval, row in zip(ypos, diagnostic.itertuples()):
        axa.text(row.percent_segments + 2, yval, f"{row.segments}/{row.total_segments}",
                 va="center", ha="left", fontsize=5.4)

    axb.hist(null["rho_chl_support"], bins=28, color=LIGHT,
             edgecolor="white", linewidth=0.4)
    axb.axvline(observed_support, color=OBSERVED, lw=1.7, label="Observed")
    axb.axvline(null["rho_chl_support"].median(), color=BLACK, lw=0.8,
                ls=(0, (3, 2)), label="Replacement median")
    axb.set(xlabel="CHL–tail-support correlation", ylabel="Replacements",
            title="Tail-shape support")
    axb.legend(frameon=False, fontsize=5.2)

    axc.hist(null["rho_chl_median_log_length"], bins=28, color=LIGHT,
             edgecolor="white", linewidth=0.4)
    axc.axvline(observed_length, color=OBSERVED, lw=1.7, label="Observed")
    axc.axvline(null["rho_chl_median_log_length"].median(), color=BLACK, lw=0.8,
                ls=(0, (3, 2)), label="Replacement median")
    axc.set(xlabel="CHL–median-length correlation", ylabel="Replacements",
            title="Median length")
    axc.legend(frameon=False, fontsize=5.2)

    axd.scatter(null["rho_chl_support"], null["rho_chl_median_log_length"],
                s=7, alpha=0.22, color=GREY, edgecolor="none")
    axd.scatter(observed_support, observed_length, marker="D", s=42,
                color=OBSERVED, edgecolor=BLACK, lw=0.35, zorder=4)
    axd.axvline(null["rho_chl_support"].median(), color=BLACK, lw=0.55, ls=(0, (3, 2)))
    axd.axhline(null["rho_chl_median_log_length"].median(), color=BLACK, lw=0.55, ls=(0, (3, 2)))
    axd.set(xlabel="CHL–tail-support correlation",
            ylabel="CHL–median-length correlation",
            title="Joint contrast")
    for ax in axs.flat:
        ax.grid(False)
        ax.tick_params(direction="out", length=2.5, pad=2)
    return finish(
        fig, "figS5_boundary_counterfactual", (7.2, 5.4),
        [audit_path, null_path, final_path, source_path],
        "The reconstruction panels are validation readbacks; the counterfactual preserves event count and therefore cannot test renewal density.",
        list(axs.flat),
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    builders = [build_s1, build_s2, build_s3, build_s4, build_s5]
    audits: dict[str, dict] = {}
    print("[0/5 0%] supplementary figure build started", flush=True)
    for index, builder in enumerate(builders, start=1):
        name = f"figureS{index}"
        print(f"[{index-1}/5 {(index-1)*20}%] building {name}", flush=True)
        audits[name] = builder()
        print(f"[{index}/5 {index*20}%] {name} complete; elapsed={time.monotonic()-started:.1f}s", flush=True)
    (OUT / "SUPPLEMENTARY_FIGURE_PACKAGE_AUDIT.json").write_text(
        json.dumps(audits, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"[5/5 100%] supplementary package complete; elapsed={time.monotonic()-started:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
