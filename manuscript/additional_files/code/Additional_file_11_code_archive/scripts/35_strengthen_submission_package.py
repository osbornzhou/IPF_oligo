#!/usr/bin/env python
"""
Add final submission-support details for reproducibility and statistical transparency.

Outputs:
- per-cohort bootstrap confidence intervals for Elastic Net external ROC/PR AUC;
- metric-definition notes clarifying cohort-level means versus pooled metrics;
- donor-level pseudobulk dot/box-style supplementary figure;
- pseudobulk display and low-cell-count policy notes.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


PROJECT_DIR = Path(__file__).resolve().parents[1]
ML_DIR = PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna"
PSEUDO_DIR = PROJECT_DIR / "results" / "single_cell_pseudobulk"
FIG_DIR = PROJECT_DIR / "manuscript" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_font(size: int, bold: bool = False):
    candidates = ["arialbd.ttf", "arial.ttf"] if bold else ["arial.ttf", "DejaVuSans.ttf"]
    for name in candidates:
        try:
            return ImageFont.truetype(name, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and sorted_values[j] == sorted_values[i]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def roc_auc_score_manual(y_true: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores).astype(float)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = average_ranks(s)
    pos_rank_sum = ranks[y == 1].sum()
    return float((pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def average_precision_manual(y_true: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores).astype(float)
    n_pos = int((y == 1).sum())
    if n_pos == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    y_sorted = y[order]
    cum_tp = np.cumsum(y_sorted == 1)
    ranks = np.arange(1, len(y_sorted) + 1)
    precisions = cum_tp / ranks
    return float(precisions[y_sorted == 1].sum() / n_pos)


def bootstrap_metric_ci(df: pd.DataFrame, metric_fn, rng: np.random.Generator, iterations: int = 2000) -> tuple[float, float, float]:
    y = df["true_label"].to_numpy(dtype=int)
    score = df["score"].to_numpy(dtype=float)
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    observed = metric_fn(y, score)
    values = []
    for _ in range(iterations):
        boot_idx = np.concatenate(
            [
                rng.choice(pos_idx, size=len(pos_idx), replace=True),
                rng.choice(neg_idx, size=len(neg_idx), replace=True),
            ]
        )
        values.append(metric_fn(y[boot_idx], score[boot_idx]))
    arr = np.asarray(values, dtype=float)
    return observed, float(np.nanpercentile(arr, 2.5)), float(np.nanpercentile(arr, 97.5))


def build_ml_ci_tables() -> None:
    pred = pd.read_csv(ML_DIR / "ml_external_validation_predictions.csv")
    pred = pred[pred["model"].eq("elastic_net")].copy()
    rng = np.random.default_rng(20260526)
    rows = []
    for label, sub in [*pred.groupby("series_id"), ("pooled_external", pred)]:
        series_id = str(label)
        roc, roc_lo, roc_hi = bootstrap_metric_ci(sub, roc_auc_score_manual, rng)
        pr, pr_lo, pr_hi = bootstrap_metric_ci(sub, average_precision_manual, rng)
        rows.append(
            {
                "series_id": series_id,
                "samples": int(len(sub)),
                "ipf_samples": int((sub["true_label"] == 1).sum()),
                "control_samples": int((sub["true_label"] == 0).sum()),
                "roc_auc": roc,
                "roc_auc_95ci_low": roc_lo,
                "roc_auc_95ci_high": roc_hi,
                "pr_auc": pr,
                "pr_auc_95ci_low": pr_lo,
                "pr_auc_95ci_high": pr_hi,
                "bootstrap_iterations": 2000,
                "bootstrap_design": "stratified resampling within IPF and control labels",
            }
        )
    pd.DataFrame(rows).to_csv(ML_DIR / "ml_external_validation_bootstrap_ci.csv", index=False, encoding="utf-8-sig")

    metric_notes = pd.DataFrame(
        [
            {
                "metric_or_term": "mean external ROC AUC / PR AUC / balanced accuracy",
                "definition": "Arithmetic mean of the corresponding per-validation-cohort metric across GSE110147, GSE150910, GSE53845, and GSE92592.",
            },
            {
                "metric_or_term": "minimum external ROC AUC",
                "definition": "Lowest per-cohort ROC AUC among the four validation cohorts.",
            },
            {
                "metric_or_term": "pooled external ROC AUC",
                "definition": "Metric computed after pooling individual external-validation predictions across cohorts; used for interpretability/calibration checks, not as the headline external metric.",
            },
            {
                "metric_or_term": "label permutation control",
                "definition": "Negative control for non-random predictive signal under the same modeling workflow; it does not exclude cohort-composition or platform effects.",
            },
            {
                "metric_or_term": "decision-curve analysis",
                "definition": "Exploratory threshold-utility check; no clinical decision threshold was predefined or claimed.",
            },
        ]
    )
    metric_notes.to_csv(ML_DIR / "ml_metric_definitions_and_interpretation_notes.csv", index=False, encoding="utf-8-sig")


def draw_pseudobulk_donor_dotplots() -> None:
    donor = pd.read_csv(PSEUDO_DIR / "single_cell_pseudobulk_donor_celltype_expression.csv")
    core = pd.read_csv(PSEUDO_DIR / "single_cell_pseudobulk_core_candidate_summary.csv")
    core = core.sort_values("ipf_minus_control_log1p_norm", ascending=True).reset_index(drop=True)
    width, height = 3000, 1850
    left, right, top, bottom = 680, 620, 170, 150
    plot_w = width - left - right
    row_h = (height - top - bottom) / len(core)
    x_min, x_max = 0.0, max(3.2, float(donor["log1p_norm_count_per_10k"].quantile(0.995)) + 0.5)

    def x_pos(value: float) -> int:
        return int(left + (value - x_min) / (x_max - x_min) * plot_w)

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title = load_font(56, bold=True)
    font = load_font(38)
    small = load_font(30)
    tiny = load_font(26)
    draw.text((left, 35), "Donor-level pseudobulk distributions for core candidates", fill=(0, 0, 0), font=title)
    draw.text((left, 100), "Dots show donor-level log1p normalized expression; bars mark group medians.", fill=(70, 70, 70), font=small)

    axis_top, axis_bottom = top - 15, height - bottom + 8
    for tick in np.arange(0, math.ceil(x_max) + 0.1, 0.5):
        x = x_pos(float(tick))
        draw.line((x, axis_top, x, axis_bottom), fill=(232, 237, 244), width=2)
        if abs(tick - round(tick)) < 1e-6:
            draw.text((x - 12, height - bottom + 25), f"{tick:g}", fill=(50, 50, 50), font=tiny)

    colors = {"Control": (115, 123, 135), "IPF": (181, 93, 96)}
    offsets = {"Control": -14, "IPF": 14}
    for i, row in core.iterrows():
        y = int(top + i * row_h + row_h * 0.5)
        label = f"{row['gene_symbol']} ({row['broad_celltype']}, {row['series_id']})"
        draw.text((35, y - 22), label, fill=(0, 0, 0), font=font)
        sub = donor[
            donor["series_id"].eq(row["series_id"])
            & donor["broad_celltype"].eq(row["broad_celltype"])
            & donor["gene_symbol"].eq(row["gene_symbol"])
        ].copy()
        for group in ["Control", "IPF"]:
            vals = sub.loc[sub["group"].eq(group), "log1p_norm_count_per_10k"].dropna().to_numpy(dtype=float)
            if len(vals) == 0:
                continue
            yy = y + offsets[group]
            for j, value in enumerate(vals):
                jitter = ((j % 7) - 3) * 1.7
                x = x_pos(float(value))
                draw.ellipse((x - 5, yy + jitter - 5, x + 5, yy + jitter + 5), fill=colors[group], outline=(255, 255, 255))
            med = float(np.median(vals))
            mx = x_pos(med)
            draw.line((mx, yy - 18, mx, yy + 18), fill=colors[group], width=6)
        draw.text(
            (width - 590, y - 18),
            f"n={int(row['ipf_donors'])}/{int(row['control_donors'])}; FDR={float(row['fdr_bh']):.1e}" if float(row["fdr_bh"]) < 1e-3 else f"n={int(row['ipf_donors'])}/{int(row['control_donors'])}; FDR={float(row['fdr_bh']):.3f}",
            fill=(45, 45, 45),
            font=small,
        )

    legend_x = left
    draw.ellipse((legend_x, height - 82, legend_x + 18, height - 64), fill=colors["Control"])
    draw.text((legend_x + 32, height - 88), "Control", fill=(60, 60, 60), font=small)
    draw.ellipse((legend_x + 190, height - 82, legend_x + 208, height - 64), fill=colors["IPF"])
    draw.text((legend_x + 222, height - 88), "IPF", fill=(60, 60, 60), font=small)
    draw.text((left + 580, height - 88), "Donor-level log1p normalized expression", fill=(0, 0, 0), font=small)

    png = PSEUDO_DIR / "single_cell_pseudobulk_core_candidate_donor_dotplots.png"
    pdf = PSEUDO_DIR / "single_cell_pseudobulk_core_candidate_donor_dotplots.pdf"
    image.save(png)
    image.save(pdf, "PDF", resolution=300.0)
    image.save(FIG_DIR / "Additional_Figure_S12_pseudobulk_donor_dotplots.png")

    figure_note = pd.DataFrame(
        [
            {
                "figure": "Additional Figure S12",
                "file_png": str((FIG_DIR / "Additional_Figure_S12_pseudobulk_donor_dotplots.png").relative_to(PROJECT_DIR)).replace("\\", "/"),
                "file_pdf": str(pdf.relative_to(PROJECT_DIR)).replace("\\", "/"),
                "description": "Donor-level dot/median plots for the core pseudobulk comparisons shown in Figure 6C.",
            }
        ]
    )
    figure_note.to_csv(PSEUDO_DIR / "single_cell_pseudobulk_donor_dotplot_manifest.csv", index=False, encoding="utf-8-sig")


def build_pseudobulk_policy_notes() -> None:
    core = pd.read_csv(PSEUDO_DIR / "single_cell_pseudobulk_core_candidate_summary.csv")
    notes = pd.DataFrame(
        [
            {
                "item": "statistical_test",
                "value": "two-sided Welch test comparing donor-level log1p normalized pseudobulk expression between IPF and control donors",
            },
            {
                "item": "fdr_scope",
                "value": "Benjamini-Hochberg correction across all tested core candidate x dataset x broad-celltype pseudobulk comparisons",
            },
            {
                "item": "display_rule_for_figure_6c",
                "value": "Displayed comparisons are direction-consistent with bulk evidence, have at least 10 IPF and 10 control donors, and pass BH-adjusted FDR < 0.01.",
            },
            {
                "item": "donor_celltype_cell_count_policy",
                "value": "No additional minimum per-donor cell-count threshold was imposed after broad-celltype aggregation; donor counts and per-donor cell-count summaries are provided for audit.",
            },
            {
                "item": "minimum_displayed_ipf_donors",
                "value": int(core["ipf_donors"].min()),
            },
            {
                "item": "minimum_displayed_control_donors",
                "value": int(core["control_donors"].min()),
            },
        ]
    )
    notes.to_csv(PSEUDO_DIR / "single_cell_pseudobulk_method_notes.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    build_ml_ci_tables()
    draw_pseudobulk_donor_dotplots()
    build_pseudobulk_policy_notes()
    print(ML_DIR / "ml_external_validation_bootstrap_ci.csv")
    print(PSEUDO_DIR / "single_cell_pseudobulk_core_candidate_donor_dotplots.pdf")


if __name__ == "__main__":
    main()
