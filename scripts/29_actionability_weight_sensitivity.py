#!/usr/bin/env python
"""
Assess sensitivity of the oligonucleotide perturbation-triage score to weight choices.

The perturbation-triage score is intended for validation planning. This script checks
whether the nominated knockdown-screening candidates remain stable under equal
weights, leave-one-layer-out scoring, and seeded +/-20% weight
perturbations.
"""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "oligonucleotide_actionability"
PLOT_DIR = OUT_DIR / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
MPLCONFIG_DIR = PROJECT_DIR / "results" / "_mpl_config"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


COMPONENTS = [
    "disease_direction_component",
    "cross_cohort_component",
    "cell_localization_component",
    "perturbation_component",
    "ml_component",
    "mirna_axis_component",
]

BASE_WEIGHTS = {
    "disease_direction_component": 1.30,
    "cross_cohort_component": 1.20,
    "cell_localization_component": 1.00,
    "perturbation_component": 1.00,
    "ml_component": 0.90,
    "mirna_axis_component": 0.80,
}

EXPECTED_KNOCKDOWN_SCREENING = {"CD24", "COL14A1", "PTGFRN", "ASPN", "CDH3"}


def score_table(df: pd.DataFrame, weights: dict[str, float], label: str) -> pd.DataFrame:
    out = df.copy()
    score = np.zeros(len(out), dtype=float)
    for comp, weight in weights.items():
        score += pd.to_numeric(out[comp], errors="coerce").fillna(0.0).to_numpy() * float(weight)
    score -= pd.to_numeric(out["context_penalty"], errors="coerce").fillna(0.0).to_numpy()
    out["sensitivity_scenario"] = label
    out["sensitivity_score"] = score
    out["sensitivity_rank"] = out["sensitivity_score"].rank(method="min", ascending=False)
    return out


def build_sensitivity() -> None:
    df = pd.read_csv(OUT_DIR / "oligonucleotide_actionability_index.csv")
    df = df[df["actionability_class"].ne("External TNIK bridge")].copy()
    for comp in COMPONENTS + ["context_penalty"]:
        df[comp] = pd.to_numeric(df[comp], errors="coerce").fillna(0.0)

    scenarios = []
    scenarios.append(score_table(df, BASE_WEIGHTS, "base_weights"))
    equal_weights = {comp: 1.0 for comp in COMPONENTS}
    scenarios.append(score_table(df, equal_weights, "equal_weights"))
    for comp in COMPONENTS:
        weights = BASE_WEIGHTS.copy()
        weights[comp] = 0.0
        scenarios.append(score_table(df, weights, f"leave_one_out__{comp}"))

    rng = np.random.default_rng(20260526)
    perturb_rows = []
    perturb_scores = []
    for i in range(1000):
        weights = {
            comp: BASE_WEIGHTS[comp] * rng.uniform(0.8, 1.2)
            for comp in COMPONENTS
        }
        label = f"weight_perturbation_{i + 1:04d}"
        perturb_rows.append({"sensitivity_scenario": label, **weights})
        perturb_scores.append(score_table(df, weights, label))

    all_scores = pd.concat(scenarios + perturb_scores, ignore_index=True)
    scenario_weights = pd.DataFrame(
        [
            {"sensitivity_scenario": "base_weights", **BASE_WEIGHTS},
            {"sensitivity_scenario": "equal_weights", **equal_weights},
            *[
                {
                    "sensitivity_scenario": f"leave_one_out__{comp}",
                    **{k: (0.0 if k == comp else v) for k, v in BASE_WEIGHTS.items()},
                }
                for comp in COMPONENTS
            ],
            *perturb_rows,
        ]
    )

    grouped = (
        all_scores.groupby(["gene_symbol", "actionability_class"], as_index=False)
        .agg(
            median_rank=("sensitivity_rank", "median"),
            mean_rank=("sensitivity_rank", "mean"),
            min_rank=("sensitivity_rank", "min"),
            max_rank=("sensitivity_rank", "max"),
            median_score=("sensitivity_score", "median"),
            top6_frequency=("sensitivity_rank", lambda x: float((x <= 6).mean())),
            top10_frequency=("sensitivity_rank", lambda x: float((x <= 10).mean())),
        )
    )
    base = all_scores[all_scores["sensitivity_scenario"].eq("base_weights")][["gene_symbol", "sensitivity_rank"]].rename(
        columns={"sensitivity_rank": "baseline_rank"}
    )
    equal = all_scores[all_scores["sensitivity_scenario"].eq("equal_weights")][["gene_symbol", "sensitivity_rank"]].rename(
        columns={"sensitivity_rank": "equal_weight_rank"}
    )
    leave_one = all_scores[all_scores["sensitivity_scenario"].str.startswith("leave_one_out__")]
    leave_one_summary = (
        leave_one.groupby("gene_symbol", as_index=False)
        .agg(
            leave_one_layer_out_min_rank=("sensitivity_rank", "min"),
            leave_one_layer_out_max_rank=("sensitivity_rank", "max"),
        )
    )
    perturb = all_scores[all_scores["sensitivity_scenario"].str.startswith("weight_perturbation_")]
    perturb_summary = (
        perturb.groupby("gene_symbol", as_index=False)
        .agg(
            perturbation_median_rank=("sensitivity_rank", "median"),
            perturbation_rank_5th_percentile=("sensitivity_rank", lambda x: float(np.percentile(x, 5))),
            perturbation_rank_95th_percentile=("sensitivity_rank", lambda x: float(np.percentile(x, 95))),
        )
    )
    rank_summary = (
        grouped.merge(base, on="gene_symbol", how="left")
        .merge(equal, on="gene_symbol", how="left")
        .merge(leave_one_summary, on="gene_symbol", how="left")
        .merge(perturb_summary, on="gene_symbol", how="left")
        .sort_values(["median_rank", "mean_rank", "gene_symbol"])
    )

    knockdown = rank_summary[rank_summary["actionability_class"].eq("Knockdown-screening candidate")].copy()
    expected_present = EXPECTED_KNOCKDOWN_SCREENING.issubset(set(knockdown["gene_symbol"]))
    stable_top10 = bool((knockdown["top10_frequency"] >= 0.95).all()) if not knockdown.empty else False
    stable_top6_any = bool((knockdown["top6_frequency"] > 0).all()) if not knockdown.empty else False
    qc = pd.DataFrame(
        [
            {"qc_item": "sensitivity_scenarios", "value": all_scores["sensitivity_scenario"].nunique(), "pass": True},
            {"qc_item": "weight_perturbation_iterations", "value": 1000, "pass": True},
            {"qc_item": "expected_knockdown_screening_present", "value": expected_present, "pass": expected_present},
            {"qc_item": "knockdown_screening_top10_frequency_ge_0_95", "value": stable_top10, "pass": stable_top10},
            {"qc_item": "knockdown_screening_has_nonzero_top6_frequency", "value": stable_top6_any, "pass": stable_top6_any},
        ]
    )

    all_scores.to_csv(OUT_DIR / "actionability_weight_sensitivity_all_scores.csv", index=False)
    scenario_weights.to_csv(OUT_DIR / "actionability_weight_sensitivity_scenario_weights.csv", index=False)
    rank_summary.to_csv(OUT_DIR / "actionability_weight_sensitivity_rank_summary.csv", index=False)
    qc.to_csv(OUT_DIR / "actionability_weight_sensitivity_qc.csv", index=False)

    plot_df = all_scores[
        all_scores["actionability_class"].isin(
            ["Knockdown-screening candidate", "Context-dependent candidate", "Restoration/pathway marker", "miRNA-axis hypothesis"]
        )
    ].copy()
    order = (
        rank_summary[rank_summary["gene_symbol"].isin(plot_df["gene_symbol"])]
        .sort_values("median_rank")["gene_symbol"]
        .tolist()
    )
    values = [plot_df.loc[plot_df["gene_symbol"].eq(g), "sensitivity_rank"].to_numpy() for g in order]
    colors = {
        "Knockdown-screening candidate": (181, 93, 96),
        "Context-dependent candidate": (217, 156, 68),
        "Restoration/pathway marker": (76, 120, 168),
        "miRNA-axis hypothesis": (78, 159, 135),
    }
    class_map = plot_df.drop_duplicates("gene_symbol").set_index("gene_symbol")["actionability_class"].to_dict()
    width, height = 2400, 1450
    left, right, top, bottom = 560, 120, 150, 180
    row_h = (height - top - bottom) / max(len(order), 1)
    x_min, x_max = 1.0, max(15.0, float(np.nanmax([np.nanmax(v) for v in values])))

    def x_pos(rank: float) -> int:
        return int(left + (rank - x_min) / (x_max - x_min) * (width - left - right))

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    try:
        font_title = ImageFont.truetype("arial.ttf", 50)
        font = ImageFont.truetype("arial.ttf", 34)
        font_small = ImageFont.truetype("arial.ttf", 28)
    except OSError:
        font_title = ImageFont.load_default()
        font = ImageFont.load_default()
        font_small = ImageFont.load_default()

    draw.text((left, 35), "Perturbation-triage score weight-sensitivity analysis", fill=(0, 0, 0), font=font_title)
    plot_top = top
    plot_bottom = height - bottom
    axis_y = plot_bottom + 20
    for rank in range(1, int(x_max) + 1):
        x = x_pos(rank)
        fill = (226, 232, 240) if rank % 2 == 0 else (241, 245, 249)
        draw.line((x, plot_top, x, plot_bottom), fill=fill, width=2)
        if rank in {1, 3, 6, 10, 15}:
            draw.text((x - 14, axis_y), str(rank), fill=(40, 40, 40), font=font_small)
    top6_x = x_pos(6)
    draw.line((top6_x, plot_top, top6_x, plot_bottom), fill=(90, 90, 90), width=4)
    draw.text((top6_x + 8, plot_top - 45), "Top 6", fill=(70, 70, 70), font=font_small)

    for idx, gene in enumerate(order):
        arr = np.sort(values[idx])
        y = int(plot_top + idx * row_h + row_h * 0.5)
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        lo, hi = np.percentile(arr, [5, 95])
        color = colors.get(class_map.get(gene, ""), (150, 150, 150))
        draw.text((40, y - 18), gene, fill=(0, 0, 0), font=font)
        draw.line((x_pos(lo), y, x_pos(hi), y), fill=(90, 90, 90), width=5)
        draw.rectangle((x_pos(q1), y - 18, x_pos(q3), y + 18), fill=color, outline=(70, 70, 70), width=2)
        draw.line((x_pos(med), y - 24, x_pos(med), y + 24), fill=(0, 0, 0), width=4)

    draw.text(
        (left, height - 95),
        "Rank across equal-weight, leave-one-layer-out, and +/-20% weight perturbation scenarios",
        fill=(0, 0, 0),
        font=font_small,
    )
    image.save(PLOT_DIR / "actionability_weight_sensitivity_rank_distribution.png")


if __name__ == "__main__":
    build_sensitivity()
    print(OUT_DIR)
