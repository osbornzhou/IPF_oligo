#!/usr/bin/env python
"""
Validate published IPF biomarker signatures under the same train/external-validation framework.

This is not a direct reproduction of external model coefficients, because most
published studies do not provide deployable model objects. Instead, each
published gene signature is used as a fixed feature set, trained only on
GSE32537, and externally validated on GSE110147/GSE150910/GSE53845/GSE92592.
"""

from __future__ import annotations

import os
import math
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
MPLCONFIG_DIR = PROJECT_DIR / "results" / "_mpl_config"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


MATRIX_DIR = PROJECT_DIR / "results" / "models" / "feature_matrices"
OUT_DIR = PROJECT_DIR / "results" / "published_signature_validation"
PLOT_DIR = OUT_DIR / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

TRAINING_DATASET = "GSE32537"
VALIDATION_DATASETS = ["GSE110147", "GSE150910", "GSE53845", "GSE92592"]


PUBLISHED_SIGNATURES = {
    "metabolism_demrg_frontiers_2023": {
        "genes": ["ENPP3", "ENTPD1", "PDE7B", "GPX3", "PNMT", "POLR3H"],
        "reference": "Identification and validation of metabolism-related hub genes in idiopathic pulmonary fibrosis",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10010493/",
    },
    "ann_six_gene_scirep_2023": {
        "genes": ["CDH3", "DIO2", "ADAMTS14", "HS6ST2", "IL13RA2", "IGFL2"],
        "reference": "Artificial neural network identified the significant genes to distinguish idiopathic pulmonary fibrosis",
        "url": "https://www.nature.com/articles/s41598-023-28536-w",
    },
    "ppi_hub_frontiers_2021": {
        "genes": ["COL1A1", "COL1A2", "COL3A1", "COL14A1", "COL15A1", "POSTN", "SPP1", "MMP1", "MMP7", "ASPN", "CDH2", "CTSK"],
        "reference": "Identification of Hub Genes and Pathways Associated With Idiopathic Pulmonary Fibrosis via Bioinformatics Analysis",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8406749/",
    },
    "explainable_ml_reported_genes_2023": {
        "genes": ["MMP7", "COL15A1", "COMP", "IL13RA2", "PAPSS2", "COL1A1", "COL17A1", "COL5A2"],
        "reference": "An explainable machine learning-driven proposal of pulmonary fibrosis biomarkers",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10049879/",
    },
    "classic_mmp_biomarkers": {
        "genes": ["MMP1", "MMP7"],
        "reference": "MMP1 and MMP7 as potential peripheral blood biomarkers in idiopathic pulmonary fibrosis",
        "url": "https://pubmed.ncbi.nlm.nih.gov/18447576/",
    },
    "cxcl14_shap_reported_2025": {
        "genes": ["CXCL14"],
        "reference": "SHAP-based feature selection reports CXCL14 as an IPF predictor",
        "url": "https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2025.1608078/full",
    },
}


def normalize_gene(x: object) -> str:
    text = "" if pd.isna(x) else str(x).strip().upper()
    return text if text and text != "NA" else ""


def load_gene_matrix(series_id: str) -> pd.DataFrame:
    path = MATRIX_DIR / f"{series_id}_gene_level_matrix.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    rename = {col: normalize_gene(col) for col in df.columns if col not in {"sample_id", "group", "label"}}
    df = df.rename(columns=rename)
    return df


def score_metrics(y_true: np.ndarray, scores: np.ndarray, preds: np.ndarray) -> dict:
    tn, fp, fn, tp = confusion_matrix(y_true, preds, labels=[0, 1]).ravel()
    return {
        "roc_auc": roc_auc_score(y_true, scores) if len(np.unique(y_true)) == 2 else np.nan,
        "pr_auc": average_precision_score(y_true, scores) if len(np.unique(y_true)) == 2 else np.nan,
        "balanced_accuracy": balanced_accuracy_score(y_true, preds),
        "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
        "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


def fit_signature_model(train_df: pd.DataFrame, genes: list[str]) -> Pipeline:
    x = train_df[genes].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    y = train_df["label"].to_numpy(dtype=int)
    model = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=20260524,
                ),
            ),
        ]
    )
    model.fit(x, y)
    return model


def main() -> None:
    matrices = {sid: load_gene_matrix(sid) for sid in [TRAINING_DATASET] + VALIDATION_DATASETS}
    ml_dir = PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna"
    our_panel = pd.read_csv(ml_dir / "ml_final_biomarker_panel.csv")
    our_panel_genes = {normalize_gene(x) for x in our_panel["feature"]}

    rows = []
    qc_rows = []
    signature_rows = []

    for name, info in PUBLISHED_SIGNATURES.items():
        genes = [normalize_gene(g) for g in info["genes"] if normalize_gene(g)]
        common = sorted(set(genes).intersection(*(set(matrices[sid].columns) for sid in matrices)))
        overlap_our = sorted(set(genes).intersection(our_panel_genes))
        required_available = max(1, math.ceil(0.8 * len(genes)))
        for gene in genes:
            signature_rows.append(
                {
                    "signature_name": name,
                    "gene_symbol": gene,
                    "in_our_final_panel": gene in our_panel_genes,
                    "reference": info["reference"],
                    "url": info["url"],
                }
            )

        qc_rows.append(
            {
                "signature_name": name,
                "reported_genes": len(genes),
                "common_available_genes": len(common),
                "required_available_genes": required_available,
                "complete_signature_available": len(common) == len(genes),
                "available_genes": ";".join(common),
                "missing_genes": ";".join(sorted(set(genes) - set(common))),
                "overlap_with_our_panel": len(overlap_our),
                "overlap_genes": ";".join(overlap_our),
                "reference": info["reference"],
                "url": info["url"],
                "qc1_signature_defined_pass": len(genes) > 0,
                "qc2_feature_availability_pass": len(common) >= required_available,
                "qc3_model_validation_possible": len(common) >= 1,
            }
        )

        if len(common) < 1:
            continue

        model = fit_signature_model(matrices[TRAINING_DATASET], common)
        for sid in VALIDATION_DATASETS:
            df = matrices[sid]
            x = df[common].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
            y = df["label"].to_numpy(dtype=int)
            scores = model.predict_proba(x)[:, 1]
            preds = model.predict(x)
            row = score_metrics(y, scores, preds)
            row.update(
                {
                    "signature_name": name,
                    "series_id": sid,
                    "feature_count": len(common),
                    "features_used": ";".join(common),
                    "reference": info["reference"],
                    "url": info["url"],
                }
            )
            rows.append(row)

    perf = pd.DataFrame(rows)
    signature_table = pd.DataFrame(signature_rows)
    qc = pd.DataFrame(qc_rows)
    if not perf.empty:
        summary = (
            perf.groupby("signature_name", as_index=False)
            .agg(
                mean_external_roc_auc=("roc_auc", "mean"),
                min_external_roc_auc=("roc_auc", "min"),
                mean_external_pr_auc=("pr_auc", "mean"),
                mean_external_balanced_accuracy=("balanced_accuracy", "mean"),
                feature_count=("feature_count", "first"),
            )
            .sort_values("mean_external_roc_auc", ascending=False)
        )
    else:
        summary = pd.DataFrame()

    qc["triple_qc_pass"] = qc["qc1_signature_defined_pass"] & qc["qc2_feature_availability_pass"] & qc["qc3_model_validation_possible"]

    comparison = summary.copy() if not summary.empty else pd.DataFrame()
    if not comparison.empty:
        comparison.insert(0, "comparator_type", "published_signature")
        own_summary = pd.read_csv(ml_dir / "ml_model_performance_external_validation_summary.csv")
        best_own = own_summary.sort_values("mean_external_roc_auc", ascending=False).iloc[0]
        own_row = pd.DataFrame(
            [
                {
                    "comparator_type": "our_discovery_only_model",
                    "signature_name": f"our_{best_own['model']}",
                    "mean_external_roc_auc": best_own["mean_external_roc_auc"],
                    "min_external_roc_auc": best_own["min_external_roc_auc"],
                    "mean_external_pr_auc": best_own["mean_external_pr_auc"],
                    "mean_external_balanced_accuracy": best_own["mean_external_balanced_accuracy"],
                    "feature_count": len(our_panel_genes),
                }
            ]
        )
        comparison = pd.concat([own_row, comparison], ignore_index=True).sort_values("mean_external_roc_auc", ascending=False)

    perf.to_csv(OUT_DIR / "published_signature_external_validation_by_dataset.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT_DIR / "published_signature_external_validation_summary.csv", index=False, encoding="utf-8-sig")
    qc.to_csv(OUT_DIR / "published_signature_validation_qc.csv", index=False, encoding="utf-8-sig")
    signature_table.to_csv(OUT_DIR / "published_signature_gene_overlap.csv", index=False, encoding="utf-8-sig")
    comparison.to_csv(OUT_DIR / "ml_vs_published_signature_comparison.csv", index=False, encoding="utf-8-sig")

    if not summary.empty:
        fig, ax = plt.subplots(figsize=(8.4, 4.8))
        plot_df = summary.sort_values("mean_external_roc_auc")
        ax.barh(plot_df["signature_name"], plot_df["mean_external_roc_auc"], color="#5b84a4")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Mean external ROC AUC")
        ax.set_title("Published IPF signatures evaluated in this project")
        fig.tight_layout()
        fig.savefig(OUT_DIR / "plots" / "published_signature_external_auc_barplot.png", dpi=180)
        plt.close(fig)

    if not comparison.empty:
        fig, ax = plt.subplots(figsize=(8.8, 5.2))
        plot_df = comparison.sort_values("mean_external_roc_auc")
        colors = ["#b55d60" if x == "our_discovery_only_model" else "#5b84a4" for x in plot_df["comparator_type"]]
        ax.barh(plot_df["signature_name"], plot_df["mean_external_roc_auc"], color=colors)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Mean external ROC AUC")
        ax.set_title("Our model versus published IPF signatures")
        fig.tight_layout()
        fig.savefig(OUT_DIR / "plots" / "ml_vs_published_signature_external_auc.png", dpi=180)
        plt.close(fig)

    with pd.ExcelWriter(OUT_DIR / "published_signature_validation_summary.xlsx", engine="openpyxl") as writer:
        qc.to_excel(writer, sheet_name="qc", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        comparison.to_excel(writer, sheet_name="model_comparison", index=False)
        perf.to_excel(writer, sheet_name="by_dataset", index=False)
        signature_table.to_excel(writer, sheet_name="gene_overlap", index=False)

    print(qc.to_string(index=False))
    print("\nSummary:")
    print(summary.to_string(index=False))
    print("\nModel comparison:")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
