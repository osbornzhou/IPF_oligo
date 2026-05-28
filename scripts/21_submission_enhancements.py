#!/usr/bin/env python
"""
Generate submission-strengthening analyses:

1. Evidence grading for miRNA-mRNA axes.
2. A transparent final target-priority table.
3. Machine-learning interpretability, calibration, and decision-curve outputs.
4. A compact cell-type disease-model summary.

All scores are evidence-prioritization scores, not causal effect estimates.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "submission_enhancements"
PLOT_DIR = OUT_DIR / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
MPLCONFIG_DIR = PROJECT_DIR / "results" / "_mpl_config"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore", message="X has feature names, but SimpleImputer was fitted without feature names")
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
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

COL_UP = "#B55D60"
COL_DOWN = "#4C78A8"
COL_ACCENT = "#3D7F6F"
COL_NEUTRAL = "#7A8793"


def read_csv(rel_path: str) -> pd.DataFrame:
    return pd.read_csv(PROJECT_DIR / rel_path)


def savefig(fig: plt.Figure, name: str) -> None:
    for ext in ["png", "pdf"]:
        fig.savefig(PLOT_DIR / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def norm01(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    vmax = values.max()
    if vmax <= 0:
        return values * 0
    return values / vmax


def grade_axes() -> tuple[pd.DataFrame, pd.DataFrame]:
    axes = read_csv("results/mirna_mrna_axes/top100_robust_mirna_mrna_axes.csv")
    axes = axes.copy()
    exact = axes["match_type"].eq("exact")
    functional = axes["support_types"].astype(str).str.contains("Functional", case=False, na=False)
    strong_assay = axes["experiments"].astype(str).str.contains("Luciferase|qRT-PCR|Western|Reporter", case=False, na=False)
    axes["axis_evidence_grade"] = np.select(
        [
            exact & functional,
            axes["match_type"].eq("arm_agnostic") & functional & strong_assay,
            axes["match_type"].eq("arm_agnostic") & functional,
        ],
        ["high_exact", "moderate_arm_agnostic_functional_assay", "exploratory_arm_agnostic"],
        default="exploratory",
    )
    axes["recommended_manuscript_role"] = np.where(
        axes["axis_evidence_grade"].eq("high_exact"),
        "main_text_prioritized_axis",
        "supplementary_exploratory_axis",
    )
    axes["validation_requirement"] = np.where(
        axes["match_type"].eq("exact"),
        "validate target repression in IPF-relevant cells",
        "confirm mature arm by small-RNA sequencing or arm-specific qPCR before mechanistic claims",
    )
    axes["axis_evidence_weight"] = np.select(
        [
            axes["axis_evidence_grade"].eq("high_exact"),
            axes["axis_evidence_grade"].eq("moderate_arm_agnostic_functional_assay"),
        ],
        [1.0, 0.5],
        default=0.25,
    )
    qc = pd.DataFrame(
        [
            {
                "total_axes": len(axes),
                "exact_axes": int((axes["match_type"] == "exact").sum()),
                "arm_agnostic_axes": int((axes["match_type"] == "arm_agnostic").sum()),
                "main_text_prioritized_axes": int((axes["recommended_manuscript_role"] == "main_text_prioritized_axis").sum()),
                "supplementary_exploratory_axes": int((axes["recommended_manuscript_role"] == "supplementary_exploratory_axis").sum()),
                "qc1_axis_table_present": len(axes) > 0,
                "qc2_exact_not_mixed_with_arm_pass": set(axes["match_type"]).issubset({"exact", "arm_agnostic"}),
                "qc3_roles_assigned_pass": axes["recommended_manuscript_role"].notna().all(),
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_axis_table_present"] & qc["qc2_exact_not_mixed_with_arm_pass"] & qc["qc3_roles_assigned_pass"]
    axes.to_csv(OUT_DIR / "mirna_mrna_axes_evidence_graded.csv", index=False, encoding="utf-8-sig")
    qc.to_csv(OUT_DIR / "mirna_mrna_axes_evidence_grading_qc.csv", index=False, encoding="utf-8-sig")
    return axes, qc


def build_priority_table(graded_axes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    robust = read_csv("results/robust_candidates/robust_mrna_candidates_strict.csv")
    ml = read_csv("results/models/ml_outputs_discovery_only_mrna/ml_final_biomarker_panel.csv")
    hubs = read_csv("results/ppi_network/string_ppi_hub_genes_robust_mrna_strict_medium_confidence.csv")
    gene_sets = read_csv("results/enrichment/gene_sets_used.csv")
    sc_top = read_csv("results/single_cell_validation/single_cell_clean_top_celltype_gene_changes.csv")

    genes = set(robust["gene_symbol"].dropna().astype(str))
    genes |= set(ml["feature"].dropna().astype(str))
    genes |= set(hubs["gene_symbol"].dropna().astype(str))
    genes |= set(graded_axes["target_gene"].dropna().astype(str))
    genes |= set(sc_top["gene_symbol"].dropna().astype(str))

    robust_idx = robust.drop_duplicates("gene_symbol").set_index("gene_symbol")
    ml_idx = ml.drop_duplicates("feature").set_index("feature")
    hubs_idx = hubs.drop_duplicates("gene_symbol").set_index("gene_symbol")
    enrichment_counts = gene_sets.groupby("gene_symbol")["gene_set"].nunique()
    sc_summary = (
        sc_top.groupby("gene_symbol")
        .agg(
            max_abs_single_cell_delta=("abs_delta", "max"),
            top_single_cell_context=("fine_celltype", lambda x: ";".join(pd.Series(x).drop_duplicates().astype(str).head(3))),
            single_cell_dataset_count=("series_id", "nunique"),
        )
    )
    axis_summary = (
        graded_axes.groupby("target_gene")
        .agg(
            exact_axis_count=("match_type", lambda x: int((x == "exact").sum())),
            arm_agnostic_axis_count=("match_type", lambda x: int((x == "arm_agnostic").sum())),
            max_axis_score=("axis_score", "max"),
            max_axis_evidence_weight=("axis_evidence_weight", "max"),
            regulating_mirnas=("candidate_mirna", lambda x: ";".join(pd.Series(x).drop_duplicates().astype(str))),
        )
    )

    rows = []
    for gene in sorted(genes):
        robust_row = robust_idx.loc[gene] if gene in robust_idx.index else None
        ml_row = ml_idx.loc[gene] if gene in ml_idx.index else None
        hub_row = hubs_idx.loc[gene] if gene in hubs_idx.index else None
        sc_row = sc_summary.loc[gene] if gene in sc_summary.index else None
        axis_row = axis_summary.loc[gene] if gene in axis_summary.index else None
        rows.append(
            {
                "gene_symbol": gene,
                "discovery_logFC": robust_row["discovery_logFC"] if robust_row is not None else (ml_row["discovery_logFC"] if ml_row is not None and "discovery_logFC" in ml_row else np.nan),
                "discovery_adj_p": robust_row["discovery_adj_p"] if robust_row is not None else (ml_row["discovery_adj_p"] if ml_row is not None and "discovery_adj_p" in ml_row else np.nan),
                "same_direction_fdr_sig_count": robust_row["same_direction_fdr_sig_count"] if robust_row is not None else (ml_row["same_direction_fdr_sig_count"] if ml_row is not None and "same_direction_fdr_sig_count" in ml_row else 0),
                "robust_score": robust_row["robust_score"] if robust_row is not None else 0,
                "in_ml_panel": gene in ml_idx.index,
                "ml_selection_frequency": ml_row["overall_selection_frequency"] if ml_row is not None else 0,
                "ml_panel_score": ml_row["panel_score"] if ml_row is not None and "panel_score" in ml_row else 0,
                "ppi_hub_rank": hub_row["hub_rank"] if hub_row is not None else np.nan,
                "ppi_hub_score": hub_row["hub_score"] if hub_row is not None else 0,
                "enrichment_gene_set_count": enrichment_counts.get(gene, 0),
                "exact_axis_count": axis_row["exact_axis_count"] if axis_row is not None else 0,
                "arm_agnostic_axis_count": axis_row["arm_agnostic_axis_count"] if axis_row is not None else 0,
                "axis_evidence_weight": axis_row["max_axis_evidence_weight"] if axis_row is not None else 0,
                "regulating_mirnas": axis_row["regulating_mirnas"] if axis_row is not None else "",
                "max_abs_single_cell_delta": sc_row["max_abs_single_cell_delta"] if sc_row is not None else 0,
                "single_cell_dataset_count": sc_row["single_cell_dataset_count"] if sc_row is not None else 0,
                "top_single_cell_context": sc_row["top_single_cell_context"] if sc_row is not None else "",
            }
        )

    table = pd.DataFrame(rows)
    table["robust_component"] = 2.0 * norm01(table["robust_score"])
    table["ml_component"] = 2.0 * table["ml_selection_frequency"].astype(float)
    table["ppi_component"] = 1.0 * norm01(table["ppi_hub_score"])
    table["enrichment_component"] = 1.0 * np.minimum(pd.to_numeric(table["enrichment_gene_set_count"], errors="coerce").fillna(0) / 3.0, 1.0)
    table["axis_component"] = 1.0 * np.minimum(table["axis_evidence_weight"].astype(float), 1.0)
    table["single_cell_component"] = 1.5 * norm01(table["max_abs_single_cell_delta"])
    table["oligo_tractability_component"] = np.where(table["discovery_logFC"].fillna(0).abs() >= 1, 0.5, 0)
    table["final_priority_score"] = table[
        [
            "robust_component",
            "ml_component",
            "ppi_component",
            "enrichment_component",
            "axis_component",
            "single_cell_component",
            "oligo_tractability_component",
        ]
    ].sum(axis=1)
    table["direction"] = np.where(table["discovery_logFC"] > 0, "upregulated", np.where(table["discovery_logFC"] < 0, "downregulated", "not_available"))
    table["suggested_oligonucleotide_strategy"] = np.select(
        [
            table["direction"].eq("upregulated"),
            table["direction"].eq("downregulated"),
        ],
        [
            "siRNA/ASO knockdown candidate",
            "restoration pathway or avoid direct knockdown",
        ],
        default="insufficient direction evidence",
    )
    table["priority_tier"] = pd.cut(
        table["final_priority_score"].rank(method="first", ascending=False),
        bins=[0, 10, 30, 80, np.inf],
        labels=["Tier 1", "Tier 2", "Tier 3", "Exploratory"],
    )
    table = table.sort_values("final_priority_score", ascending=False)

    qc = pd.DataFrame(
        [
            {
                "candidate_genes_scored": len(table),
                "tier1_genes": int((table["priority_tier"] == "Tier 1").sum()),
                "genes_with_ml_panel_evidence": int(table["in_ml_panel"].sum()),
                "genes_with_exact_axis": int((table["exact_axis_count"] > 0).sum()),
                "genes_with_single_cell_evidence": int((table["max_abs_single_cell_delta"] > 0).sum()),
                "qc1_gene_universe_pass": len(table) >= len(robust),
                "qc2_score_nonnegative_pass": bool((table["final_priority_score"] >= 0).all()),
                "qc3_tiers_assigned_pass": table["priority_tier"].notna().all(),
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_gene_universe_pass"] & qc["qc2_score_nonnegative_pass"] & qc["qc3_tiers_assigned_pass"]
    table.to_csv(OUT_DIR / "final_target_priority_integrated.csv", index=False, encoding="utf-8-sig")
    qc.to_csv(OUT_DIR / "final_target_priority_qc.csv", index=False, encoding="utf-8-sig")

    top = table.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.barh(top["gene_symbol"], top["final_priority_score"], color=np.where(top["direction"].eq("upregulated"), COL_UP, COL_DOWN))
    ax.set_xlabel("Integrated priority score")
    ax.set_title("Top integrated oligonucleotide-candidate targets")
    savefig(fig, "final_target_priority_top20")
    return table, qc


def load_external_common_matrices() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    datasets = ["GSE110147", "GSE150910", "GSE53845", "GSE92592"]
    frames = []
    for sid in datasets:
        df = read_csv(f"results/models/feature_matrices/{sid}_common_discovery_only_mrna_matrix.csv").copy()
        df["series_id"] = sid
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True)
    feature_cols = [c for c in merged.columns if c not in {"sample_id", "group", "label", "series_id"}]
    y = merged["label"].astype(int).to_numpy()
    return merged, y, feature_cols


def ml_explainability_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model_path = PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna" / "artifacts" / "elastic_net_final_model.joblib"
    model = joblib.load(model_path)
    ext, y, feature_cols = load_external_common_matrices()
    x = ext[feature_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    base_scores = model.predict_proba(x)[:, 1]
    base_auc = roc_auc_score(y, base_scores)

    perm = permutation_importance(model, x, y, n_repeats=20, random_state=20260524, scoring="roc_auc", n_jobs=1)
    importance = pd.DataFrame(
        {
            "feature": feature_cols,
            "permutation_auc_drop_mean": perm.importances_mean,
            "permutation_auc_drop_sd": perm.importances_std,
        }
    ).sort_values("permutation_auc_drop_mean", ascending=False)
    importance.to_csv(OUT_DIR / "ml_external_permutation_importance.csv", index=False, encoding="utf-8-sig")

    frac_pos, mean_pred = calibration_curve(y, base_scores, n_bins=10, strategy="quantile")
    calibration = pd.DataFrame({"mean_predicted_probability": mean_pred, "observed_ipf_fraction": frac_pos})
    calibration.to_csv(OUT_DIR / "ml_external_calibration_curve.csv", index=False, encoding="utf-8-sig")

    thresholds = np.arange(0.05, 0.96, 0.05)
    n = len(y)
    prevalence = y.mean()
    dca_rows = []
    for pt in thresholds:
        pred = base_scores >= pt
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        nb_model = tp / n - fp / n * (pt / (1 - pt))
        nb_all = prevalence - (1 - prevalence) * (pt / (1 - pt))
        dca_rows.append({"threshold": pt, "model_net_benefit": nb_model, "treat_all_net_benefit": nb_all, "treat_none_net_benefit": 0.0})
    dca = pd.DataFrame(dca_rows)
    dca.to_csv(OUT_DIR / "ml_external_decision_curve.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    top = importance.head(20).iloc[::-1]
    axes[0].barh(top["feature"], top["permutation_auc_drop_mean"], xerr=top["permutation_auc_drop_sd"], color=COL_ACCENT)
    axes[0].set_xlabel("External AUC decrease")
    axes[0].set_title("Permutation importance")
    axes[1].plot([0, 1], [0, 1], ls="--", c=COL_NEUTRAL, lw=1)
    axes[1].plot(mean_pred, frac_pos, marker="o", c=COL_UP)
    axes[1].set_xlabel("Mean predicted probability")
    axes[1].set_ylabel("Observed IPF fraction")
    axes[1].set_title("External calibration")
    axes[2].plot(dca["threshold"], dca["model_net_benefit"], label="model", c=COL_UP)
    axes[2].plot(dca["threshold"], dca["treat_all_net_benefit"], label="treat all", c=COL_NEUTRAL)
    axes[2].plot(dca["threshold"], dca["treat_none_net_benefit"], label="treat none", c="#333333", ls="--")
    axes[2].set_xlabel("Threshold probability")
    axes[2].set_ylabel("Net benefit")
    axes[2].set_title("Decision curve")
    axes[2].legend(frameon=False)
    fig.suptitle(f"Machine-learning interpretability and utility checks (pooled external AUC={base_auc:.3f})", y=1.03)
    savefig(fig, "ml_interpretability_calibration_decision_curve")

    qc = pd.DataFrame(
        [
            {
                "model_artifact": str(model_path),
                "external_samples": n,
                "external_features": len(feature_cols),
                "pooled_external_auc": base_auc,
                "top_permutation_feature": importance.iloc[0]["feature"],
                "calibration_bins": len(calibration),
                "decision_curve_thresholds": len(dca),
                "qc1_model_artifact_pass": model_path.exists(),
                "qc2_importance_nonempty_pass": len(importance) == len(feature_cols),
                "qc3_probability_range_pass": bool(((base_scores >= 0) & (base_scores <= 1)).all()),
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_model_artifact_pass"] & qc["qc2_importance_nonempty_pass"] & qc["qc3_probability_range_pass"]
    qc.to_csv(OUT_DIR / "ml_interpretability_qc.csv", index=False, encoding="utf-8-sig")
    return importance, calibration, qc


def biological_model_summary(priority: pd.DataFrame, graded_axes: pd.DataFrame) -> pd.DataFrame:
    sc_top = read_csv("results/single_cell_validation/single_cell_clean_top_celltype_gene_changes.csv")
    rows = [
        {
            "module": "Macrophage/myeloid activation",
            "supporting_genes": "SPP1",
            "main_cellular_context": "GSE135893 macrophages; GSE136831 macrophage/myeloid",
            "evidence": "largest positive single-cell disease-control shifts; STRING hub evidence",
            "oligonucleotide_angle": "candidate knockdown or pathway-modulation target after functional validation",
        },
        {
            "module": "Stromal matrix remodeling",
            "supporting_genes": "COL1A1;COL1A2;COL3A1;POSTN;COL14A1;ASPN",
            "main_cellular_context": "myofibroblast/stromal and mesenchymal compartments",
            "evidence": "robust mRNA upregulation; PPI hub and enrichment support; single-cell localization",
            "oligonucleotide_angle": "upregulated mRNAs are conceptually compatible with siRNA/ASO knockdown screening",
        },
        {
            "module": "Stromal antioxidant/metabolic loss",
            "supporting_genes": "GPX3",
            "main_cellular_context": "fibroblast and myofibroblast populations",
            "evidence": "ML stable feature; robust downregulation; strong single-cell negative delta",
            "oligonucleotide_angle": "not a direct knockdown candidate; suggests restoration or pathway-protection hypothesis",
        },
        {
            "module": "Epithelial/ciliary remodeling",
            "supporting_genes": "TPPP3;DNAI1;TEKT1;RSPH4A;MNS1;RPGRIP1L",
            "main_cellular_context": "epithelial and ciliated/axonemal programs",
            "evidence": "cilium/axoneme enrichment and PPI hub evidence; selected epithelial single-cell localization",
            "oligonucleotide_angle": "candidate epithelial program for mechanistic validation, not yet therapeutically validated",
        },
        {
            "module": "Higher-confidence miRNA-target axis",
            "supporting_genes": "hsa-miR-375 -> CLDN1/MNS1/RPGRIP1L",
            "main_cellular_context": "bulk miRNA/mRNA inverse evidence; cell context requires targeted validation",
            "evidence": f"{int((graded_axes['axis_evidence_grade'] == 'high_exact').sum())} exact mature-miRNA axes; arm-agnostic axes downgraded to supplementary exploratory evidence",
            "oligonucleotide_angle": "miRNA replacement or target-axis perturbation hypothesis requiring isoform-specific validation",
        },
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "biological_model_summary.csv", index=False, encoding="utf-8-sig")
    qc = pd.DataFrame(
        [
            {
                "modules": len(summary),
                "priority_table_genes": len(priority),
                "single_cell_rows_used": len(sc_top),
                "exact_axes": int((graded_axes["match_type"] == "exact").sum()),
                "qc1_modules_present": len(summary) >= 4,
                "qc2_no_arm_claim_as_high_pass": bool((graded_axes.loc[graded_axes["match_type"].eq("arm_agnostic"), "recommended_manuscript_role"] == "supplementary_exploratory_axis").all()),
                "qc3_oligo_angle_present": summary["oligonucleotide_angle"].notna().all(),
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_modules_present"] & qc["qc2_no_arm_claim_as_high_pass"] & qc["qc3_oligo_angle_present"]
    qc.to_csv(OUT_DIR / "biological_model_summary_qc.csv", index=False, encoding="utf-8-sig")

    def wrap_text(value: str, width: int) -> str:
        import textwrap

        wrapped_lines = []
        for line in str(value).splitlines():
            wrapped_lines.extend(textwrap.wrap(line, width=width, break_long_words=False) or [""])
        return "\n".join(wrapped_lines)

    fig, ax = plt.subplots(figsize=(12, 7.4))
    ax.axis("off")
    ax.text(0.02, 0.95, "Biological interpretation model for prioritized IPF candidates", fontweight="bold", fontsize=12, va="top")
    ax.text(
        0.02,
        0.905,
        "Each module links robust bulk signals, machine-learning priority, single-cell context, and a cautious oligonucleotide development angle.",
        fontsize=8,
        color="#555555",
        va="top",
    )
    y_positions = np.linspace(0.78, 0.14, len(summary))
    for y, (_, row) in zip(y_positions, summary.iterrows()):
        ax.add_patch(
            plt.Rectangle(
                (0.02, y - 0.07),
                0.96,
                0.12,
                facecolor="#F7F9FA",
                edgecolor="#D5D8DC",
                linewidth=0.7,
            )
        )
        ax.text(0.04, y + 0.025, wrap_text(row["module"], 32), fontweight="bold", fontsize=8.6, va="top")
        details = (
            f"Genes/axis: {row['supporting_genes']}\n"
            f"Cell context: {row['main_cellular_context']}\n"
            f"Oligonucleotide angle: {row['oligonucleotide_angle']}"
        )
        ax.text(0.34, y + 0.032, wrap_text(details, 100), fontsize=7.0, va="top", linespacing=1.18)
    savefig(fig, "biological_model_summary")
    return summary


def export_excel() -> None:
    with pd.ExcelWriter(OUT_DIR / "submission_enhancements_summary.xlsx", engine="openpyxl") as writer:
        for name in [
            "mirna_mrna_axes_evidence_graded",
            "mirna_mrna_axes_evidence_grading_qc",
            "final_target_priority_integrated",
            "final_target_priority_qc",
            "ml_external_permutation_importance",
            "ml_external_calibration_curve",
            "ml_external_decision_curve",
            "ml_interpretability_qc",
            "biological_model_summary",
            "biological_model_summary_qc",
        ]:
            pd.read_csv(OUT_DIR / f"{name}.csv").to_excel(writer, sheet_name=name[:31], index=False)


def main() -> None:
    graded_axes, axis_qc = grade_axes()
    priority, priority_qc = build_priority_table(graded_axes)
    ml_importance, calibration, ml_qc = ml_explainability_outputs()
    model_summary = biological_model_summary(priority, graded_axes)
    export_excel()
    print("Submission enhancement outputs written to:")
    print(OUT_DIR)
    print("\nQC summary:")
    print(pd.concat([axis_qc, priority_qc, ml_qc], axis=0, ignore_index=True).to_string(index=False))
    print("\nTop integrated targets:")
    print(priority.head(15)[["gene_symbol", "final_priority_score", "priority_tier", "direction", "in_ml_panel", "exact_axis_count", "max_abs_single_cell_delta"]].to_string(index=False))


if __name__ == "__main__":
    main()
