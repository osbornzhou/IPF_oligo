#!/usr/bin/env python
"""
Build an oligonucleotide-focused perturbation-triage score.

The score is an experimental-prioritization aid, not a therapeutic-readiness
or causal-effect score. It separates knockdown-screening upregulated candidates,
context-dependent candidates, restoration/pathway markers, miRNA-axis
hypotheses, and the externally motivated TNIK bridge.
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

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ModuleNotFoundError:
    HAS_MATPLOTLIB = False


if HAS_MATPLOTLIB:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

COLORS = {
    "Knockdown-screening candidate": "#B55D60",
    "Context-dependent candidate": "#D99C44",
    "Restoration/pathway marker": "#4C78A8",
    "miRNA-axis hypothesis": "#4E9F87",
    "External TNIK bridge": "#7B5EA7",
}

KNOCKDOWN_SCREENING = {"COL14A1", "CD24", "PTGFRN", "ASPN", "CDH3"}
CONTEXT_DEPENDENT = {"SPP1", "COL1A1", "COL3A1", "POSTN"}
RESTORATION_MARKERS = {"GPX3", "NECAB1"}


def read_csv(rel_path: str) -> pd.DataFrame:
    return pd.read_csv(PROJECT_DIR / rel_path)


def norm01(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    vmax = values.max()
    vmin = values.min()
    if vmax == vmin:
        return values * 0
    return (values - vmin) / (vmax - vmin)


def classify_gene(row: pd.Series) -> str:
    gene = str(row["gene_symbol"])
    direction = str(row.get("direction", ""))
    exact_axis = float(row.get("exact_axis_count", 0) or 0)
    if gene == "TNIK":
        return "External TNIK bridge"
    if gene in KNOCKDOWN_SCREENING:
        return "Knockdown-screening candidate"
    if gene in CONTEXT_DEPENDENT:
        return "Context-dependent candidate"
    if gene in RESTORATION_MARKERS:
        return "Restoration/pathway marker"
    if exact_axis > 0:
        return "miRNA-axis hypothesis"
    if direction == "downregulated":
        return "Restoration/pathway marker"
    if direction == "upregulated" and bool(row.get("in_ml_panel", False)):
        return "Knockdown-screening candidate"
    if direction == "upregulated":
        return "Context-dependent candidate"
    return "Context-dependent candidate"


def short_evidence(row: pd.Series) -> str:
    layers = []
    if float(row.get("same_direction_fdr_sig_count", 0) or 0) >= 2:
        layers.append("robust")
    if bool(row.get("in_ml_panel", False)):
        layers.append("ML")
    if float(row.get("exact_axis_count", 0) or 0) > 0:
        layers.append("exact miRNA axis")
    elif float(row.get("arm_agnostic_axis_count", 0) or 0) > 0:
        layers.append("exploratory miRNA axis")
    if float(row.get("max_abs_single_cell_delta", 0) or 0) > 0:
        layers.append("single-cell")
    if float(row.get("perturbation_priority_proxy_score", 0) or 0) > 0:
        layers.append("perturbation-priority")
    if str(row.get("gene_symbol")) == "TNIK":
        layers.append("external bridge")
    return "; ".join(layers)


def build_actionability() -> tuple[pd.DataFrame, pd.DataFrame]:
    priority = read_csv("results/submission_enhancements/final_target_priority_integrated.csv")
    perturb = read_csv("results/mechanistic_extension/perturbation_priority_proxy_summary.csv")
    axes = read_csv("results/submission_enhancements/mirna_mrna_axes_evidence_graded.csv")

    perturb_cols = [
        "target_gene",
        "perturbation_priority_proxy_score",
        "oligonucleotide_interpretation",
        "claim_strength",
    ]
    merged = priority.merge(
        perturb[perturb_cols],
        left_on="gene_symbol",
        right_on="target_gene",
        how="left",
    )
    merged["perturbation_priority_proxy_score"] = pd.to_numeric(
        merged["perturbation_priority_proxy_score"], errors="coerce"
    ).fillna(0.0)
    merged["actionability_class"] = merged.apply(classify_gene, axis=1)

    merged["disease_direction_component"] = np.select(
        [
            merged["direction"].eq("upregulated"),
            merged["direction"].eq("downregulated"),
        ],
        [1.0, 0.35],
        default=0.0,
    )
    merged["cross_cohort_component"] = np.minimum(
        pd.to_numeric(merged["same_direction_fdr_sig_count"], errors="coerce").fillna(0) / 4.0,
        1.0,
    )
    merged["cell_localization_component"] = np.minimum(
        pd.to_numeric(merged["single_cell_dataset_count"], errors="coerce").fillna(0) / 2.0,
        1.0,
    )
    merged["perturbation_component"] = norm01(merged["perturbation_priority_proxy_score"].clip(lower=0))
    merged["mirna_axis_component"] = np.minimum(
        pd.to_numeric(merged["exact_axis_count"], errors="coerce").fillna(0)
        + 0.25 * pd.to_numeric(merged["arm_agnostic_axis_count"], errors="coerce").fillna(0),
        1.0,
    )
    merged["ml_component"] = pd.to_numeric(merged["ml_selection_frequency"], errors="coerce").fillna(0.0)
    merged["context_penalty"] = np.where(
        merged["gene_symbol"].isin(["COL1A1", "COL3A1", "SPP1", "POSTN"]) | merged["direction"].eq("downregulated"),
        0.15,
        0.0,
    )
    merged["oligo_actionability_score"] = (
        1.30 * merged["disease_direction_component"]
        + 1.20 * merged["cross_cohort_component"]
        + 1.00 * merged["cell_localization_component"]
        + 1.00 * merged["perturbation_component"]
        + 0.90 * merged["ml_component"]
        + 0.80 * merged["mirna_axis_component"]
        - merged["context_penalty"]
    )
    merged.loc[merged["gene_symbol"].eq("TNIK"), "oligo_actionability_score"] = 1.5
    merged["actionability_evidence_layers"] = merged.apply(short_evidence, axis=1)
    merged["recommended_validation_route"] = np.select(
        [
            merged["actionability_class"].eq("Knockdown-screening candidate"),
            merged["actionability_class"].eq("Context-dependent candidate"),
            merged["actionability_class"].eq("Restoration/pathway marker"),
            merged["actionability_class"].eq("miRNA-axis hypothesis"),
            merged["actionability_class"].eq("External TNIK bridge"),
        ],
        [
            "siRNA/ASO knockdown screen in localized disease-relevant cells",
            "context-specific knockdown or pathway-modulation screen",
            "restoration/pathway readout; avoid simple knockdown framing",
            "miRNA mimic or axis-specific repression validation",
            "test TNIK reagent as external Wnt/TNIK bridge, not primary discovery candidate",
        ],
        default="exploratory follow-up",
    )

    show_genes = [
        "COL14A1",
        "POSTN",
        "CD24",
        "PTGFRN",
        "ASPN",
        "CDH3",
        "SPP1",
        "COL1A1",
        "COL3A1",
        "GPX3",
        "NECAB1",
        "CLDN1",
        "MNS1",
        "RPGRIP1L",
        "TNIK",
    ]
    table = merged[merged["gene_symbol"].isin(show_genes)].copy()
    if "TNIK" not in set(table["gene_symbol"]):
        tnik = perturb[perturb["target_gene"].eq("TNIK")].head(1)
        table = pd.concat(
            [
                table,
                pd.DataFrame(
                    [
                        {
                            "gene_symbol": "TNIK",
                            "direction": "context-dependent",
                            "actionability_class": "External TNIK bridge",
                            "oligo_actionability_score": 1.5,
                            "actionability_evidence_layers": "validation DE; Wnt/TNIK bridge; external reagent",
                            "recommended_validation_route": "test TNIK reagent as external Wnt/TNIK bridge, not primary discovery candidate",
                            "top_single_cell_context": "",
                            "perturbation_priority_proxy_score": float(tnik["perturbation_priority_proxy_score"].iloc[0]) if not tnik.empty else np.nan,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    present = set(table["gene_symbol"])
    for _, axis in axes[axes["recommended_manuscript_role"].eq("main_text_prioritized_axis")].iterrows():
        gene = axis["target_gene"]
        if gene in present:
            continue
        table = pd.concat(
            [
                table,
                pd.DataFrame(
                    [
                        {
                            "gene_symbol": gene,
                            "direction": "upregulated",
                            "actionability_class": "miRNA-axis hypothesis",
                            "oligo_actionability_score": 2.2,
                            "actionability_evidence_layers": "exact miRNA axis",
                            "recommended_validation_route": "miRNA mimic or axis-specific repression validation",
                            "top_single_cell_context": "",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

    table = table.sort_values(
        ["actionability_class", "oligo_actionability_score"],
        ascending=[True, False],
    )
    table.to_csv(OUT_DIR / "oligonucleotide_actionability_index.csv", index=False, encoding="utf-8-sig")

    qc = pd.DataFrame(
        [
            {
                "candidates_in_actionability_map": int(len(table)),
                "actionability_classes": int(table["actionability_class"].nunique()),
            "knockdown_screening_candidates": int((table["actionability_class"] == "Knockdown-screening candidate").sum()),
                "restoration_or_pathway_markers": int((table["actionability_class"] == "Restoration/pathway marker").sum()),
                "mirna_axis_hypotheses": int((table["actionability_class"] == "miRNA-axis hypothesis").sum()),
                "external_tnik_bridge_present": bool((table["actionability_class"] == "External TNIK bridge").any()),
                "qc1_input_priority_table_pass": len(priority) > 100,
                "qc2_class_assignment_pass": table["actionability_class"].notna().all(),
                "qc3_score_nonmissing_pass": table["oligo_actionability_score"].notna().all(),
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_input_priority_table_pass"] & qc["qc2_class_assignment_pass"] & qc["qc3_score_nonmissing_pass"]
    qc.to_csv(OUT_DIR / "oligonucleotide_actionability_qc.csv", index=False, encoding="utf-8-sig")
    return table, qc


def plot_actionability(table: pd.DataFrame) -> None:
    if not HAS_MATPLOTLIB:
        plot_actionability_pil(table)
        return

    plot = table.copy()
    order = {
        "Knockdown-screening candidate": 0,
        "Context-dependent candidate": 1,
        "Restoration/pathway marker": 2,
        "miRNA-axis hypothesis": 3,
        "External TNIK bridge": 4,
    }
    plot["class_order"] = plot["actionability_class"].map(order)
    plot = plot.sort_values(["class_order", "oligo_actionability_score"], ascending=[False, True])

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    colors = plot["actionability_class"].map(COLORS).fillna("#7A8793")
    ax.barh(plot["gene_symbol"], plot["oligo_actionability_score"], color=colors)
    ax.set_xlabel("Oligonucleotide validation-planning score")
    ax.set_title("Oligonucleotide-focused perturbation-triage map")
    ax.grid(axis="x", color="#E0E0E0", lw=0.5)
    for y, (_, row) in enumerate(plot.iterrows()):
        ax.text(
            row["oligo_actionability_score"] + 0.05,
            y,
            row["actionability_class"].replace(" candidate", ""),
            va="center",
            fontsize=6.5,
            color="#333333",
        )
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=color, label=label)
        for label, color in COLORS.items()
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "figure_7_oligonucleotide_actionability_map.png", dpi=300)
    fig.savefig(PLOT_DIR / "figure_7_oligonucleotide_actionability_map.pdf")
    plt.close(fig)


def pil_font(size: int, bold: bool = False):
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def plot_actionability_pil(table: pd.DataFrame) -> None:
    plot = table.copy()
    order = {
        "Knockdown-screening candidate": 0,
        "Context-dependent candidate": 1,
        "Restoration/pathway marker": 2,
        "miRNA-axis hypothesis": 3,
        "External TNIK bridge": 4,
    }
    plot["class_order"] = plot["actionability_class"].map(order)
    plot = plot.sort_values(["class_order", "oligo_actionability_score"], ascending=[False, True])

    width, height = 2460, 1980
    left, top, right, bottom = 360, 190, 1780, 1780
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = pil_font(42, True)
    label_font = pil_font(27)
    tick_font = pil_font(24)
    small_font = pil_font(22)

    draw.text((left, 70), "Oligonucleotide-focused perturbation-triage map", fill=(20, 20, 20), font=title_font)
    max_score = max(float(plot["oligo_actionability_score"].max()), 4.5)
    x_max = 4.8
    rows = list(plot.iterrows())
    row_gap = (bottom - top) / max(1, len(rows))
    bar_h = max(38, int(row_gap * 0.62))

    for x in np.arange(0, x_max + 0.01, 1.0):
        px = left + int((right - left) * x / x_max)
        draw.line((px, top - 12, px, bottom + 8), fill=(224, 224, 224), width=2)
        draw.text((px - 10, bottom + 26), f"{x:.0f}", fill=(60, 60, 60), font=tick_font)
    draw.line((left, bottom + 8, right, bottom + 8), fill=(80, 80, 80), width=2)
    draw.text((left + 360, bottom + 80), "Oligonucleotide validation-planning score", fill=(35, 35, 35), font=label_font)

    for i, (_, row) in enumerate(rows):
        y = int(top + i * row_gap + (row_gap - bar_h) / 2)
        gene = str(row["gene_symbol"])
        score = float(row["oligo_actionability_score"])
        cls = str(row["actionability_class"])
        x1 = left + int((right - left) * score / x_max)
        color = COLORS.get(cls, "#7A8793")
        draw.text((40, y + 5), gene, fill=(20, 20, 20), font=label_font)
        draw.rectangle((left, y, x1, y + bar_h), fill=color)
        draw.text((x1 + 18, y + 5), cls.replace(" candidate", ""), fill=(40, 40, 40), font=small_font)

    out_png = PLOT_DIR / "figure_7_oligonucleotide_actionability_map.png"
    image.save(out_png, dpi=(300, 300))


def main() -> None:
    table, qc = build_actionability()
    plot_actionability(table)
    print(OUT_DIR)
    print(table[["gene_symbol", "actionability_class", "oligo_actionability_score"]].to_string(index=False))
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
