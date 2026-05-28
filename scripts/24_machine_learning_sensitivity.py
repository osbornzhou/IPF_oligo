#!/usr/bin/env python
"""
Machine-learning sensitivity analyses for the IPF manuscript.

Adds two model-sensitivity checks:
1. Leave-one-validation-cohort-out sensitivity using the already locked
   discovery-trained Elastic Net external predictions.
2. Random 25-gene robust-panel baselines trained in GSE32537 and externally
   validated on the four held-out cohorts.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "ml_sensitivity"
PLOT_DIR = OUT_DIR / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

os.environ["MPLCONFIGDIR"] = str(PROJECT_DIR / "results" / "_mpl_config")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import joblib
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


RANDOM_STATE = 20260525
N_RANDOM_PANELS = 500
TRAINING_DATASET = "GSE32537"
VALIDATION_DATASETS = ["GSE110147", "GSE150910", "GSE53845", "GSE92592"]
MATRIX_DIR = PROJECT_DIR / "results" / "models" / "feature_matrices"
ML_DIR = PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna"


def load_matrix(series_id: str) -> pd.DataFrame:
    path = MATRIX_DIR / f"{series_id}_common_robust_mrna_matrix.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def metric_summary(y_true: np.ndarray, score: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": roc_auc_score(y_true, score) if len(np.unique(y_true)) == 2 else np.nan,
        "pr_auc": average_precision_score(y_true, score) if len(np.unique(y_true)) == 2 else np.nan,
        "balanced_accuracy": balanced_accuracy_score(y_true, pred),
    }


def make_fixed_elastic_net(c_value: float, l1_ratio: float, random_state: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    class_weight="balanced",
                    C=c_value,
                    l1_ratio=l1_ratio,
                    max_iter=8000,
                    random_state=random_state,
                ),
            ),
        ]
    )


def extract_final_elastic_net_params() -> tuple[float, float]:
    artifact = ML_DIR / "artifacts" / "elastic_net_final_model.joblib"
    model = joblib.load(artifact)
    clf = model.named_steps["clf"]
    return float(clf.C), float(clf.l1_ratio)


def leave_one_validation_sensitivity() -> tuple[pd.DataFrame, pd.DataFrame]:
    pred = pd.read_csv(ML_DIR / "ml_external_validation_predictions.csv")
    pred = pred[pred["model"] == "elastic_net"].copy()

    per_cohort_rows = []
    for series_id, df in pred.groupby("series_id"):
        y = df["true_label"].to_numpy()
        score = df["score"].to_numpy()
        pred_label = df["prediction"].to_numpy()
        row = metric_summary(y, score, pred_label)
        row.update({"series_id": series_id, "samples": len(df), "ipf": int(y.sum()), "control": int((y == 0).sum())})
        per_cohort_rows.append(row)
    per_cohort = pd.DataFrame(per_cohort_rows).sort_values("series_id")

    rows = []
    for omitted in VALIDATION_DATASETS:
        keep = pred[pred["series_id"] != omitted]
        y = keep["true_label"].to_numpy()
        score = keep["score"].to_numpy()
        pred_label = keep["prediction"].to_numpy()
        pooled = metric_summary(y, score, pred_label)
        kept_cohorts = [x for x in VALIDATION_DATASETS if x != omitted]
        cohort_metrics = per_cohort[per_cohort["series_id"].isin(kept_cohorts)]
        pooled.update(
            {
                "omitted_validation_cohort": omitted,
                "kept_validation_cohorts": ";".join(kept_cohorts),
                "samples": len(keep),
                "mean_cohort_roc_auc": cohort_metrics["roc_auc"].mean(),
                "min_cohort_roc_auc": cohort_metrics["roc_auc"].min(),
                "mean_cohort_pr_auc": cohort_metrics["pr_auc"].mean(),
                "mean_cohort_balanced_accuracy": cohort_metrics["balanced_accuracy"].mean(),
            }
        )
        rows.append(pooled)

    loo = pd.DataFrame(rows)
    per_cohort.to_csv(OUT_DIR / "elastic_net_external_per_cohort_performance.csv", index=False, encoding="utf-8-sig")
    loo.to_csv(OUT_DIR / "leave_one_validation_cohort_sensitivity.csv", index=False, encoding="utf-8-sig")
    return per_cohort, loo


def random_robust_gene_panel_baseline() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RANDOM_STATE)
    c_value, l1_ratio = extract_final_elastic_net_params()

    matrices = {series_id: load_matrix(series_id) for series_id in [TRAINING_DATASET] + VALIDATION_DATASETS}
    feature_cols = [
        c
        for c in matrices[TRAINING_DATASET].columns
        if c not in {"sample_id", "group", "label"}
    ]
    common_features = set(feature_cols)
    for df in matrices.values():
        common_features &= set(c for c in df.columns if c not in {"sample_id", "group", "label"})
    common_features = sorted(common_features)

    panel = pd.read_csv(ML_DIR / "ml_final_biomarker_panel.csv")
    observed_panel = [g for g in panel["feature"].astype(str).head(25).tolist() if g in common_features]
    if len(observed_panel) < 10:
        raise ValueError("Observed panel has too few genes in common robust matrix.")

    train_df = matrices[TRAINING_DATASET]
    y_train = train_df["label"].to_numpy(dtype=int)

    rows = []
    pred_rows = []

    def fit_and_score(features: list[str], label: str, panel_id: int, random_state: int) -> dict[str, float]:
        model = make_fixed_elastic_net(c_value, l1_ratio, random_state)
        model.fit(train_df[features].to_numpy(dtype=float), y_train)
        cohort_rocs = []
        cohort_prs = []
        cohort_baccs = []
        pooled_y = []
        pooled_score = []
        pooled_pred = []
        for series_id in VALIDATION_DATASETS:
            df = matrices[series_id]
            y = df["label"].to_numpy(dtype=int)
            x = df[features].to_numpy(dtype=float)
            score = model.predict_proba(x)[:, 1]
            pred = model.predict(x)
            m = metric_summary(y, score, pred)
            cohort_rocs.append(m["roc_auc"])
            cohort_prs.append(m["pr_auc"])
            cohort_baccs.append(m["balanced_accuracy"])
            pooled_y.extend(y.tolist())
            pooled_score.extend(score.tolist())
            pooled_pred.extend(pred.tolist())
            pred_rows.append(
                {
                    "panel_type": label,
                    "panel_id": panel_id,
                    "series_id": series_id,
                    "roc_auc": m["roc_auc"],
                    "pr_auc": m["pr_auc"],
                    "balanced_accuracy": m["balanced_accuracy"],
                }
            )
        pooled = metric_summary(np.array(pooled_y), np.array(pooled_score), np.array(pooled_pred))
        return {
            "panel_type": label,
            "panel_id": panel_id,
            "feature_count": len(features),
            "features": ";".join(features),
            "mean_external_roc_auc": float(np.mean(cohort_rocs)),
            "min_external_roc_auc": float(np.min(cohort_rocs)),
            "mean_external_pr_auc": float(np.mean(cohort_prs)),
            "mean_external_balanced_accuracy": float(np.mean(cohort_baccs)),
            "pooled_external_roc_auc": pooled["roc_auc"],
            "pooled_external_pr_auc": pooled["pr_auc"],
            "pooled_external_balanced_accuracy": pooled["balanced_accuracy"],
        }

    rows.append(fit_and_score(observed_panel, "observed_25_gene_panel_refit", 0, RANDOM_STATE))
    for i in range(1, N_RANDOM_PANELS + 1):
        features = rng.choice(common_features, size=len(observed_panel), replace=False).tolist()
        rows.append(fit_and_score(features, "random_robust_25_gene_panel", i, RANDOM_STATE + i))

    summary = pd.DataFrame(rows)
    cohort_perf = pd.DataFrame(pred_rows)

    observed = summary[summary["panel_type"] == "observed_25_gene_panel_refit"].iloc[0]
    random_df = summary[summary["panel_type"] == "random_robust_25_gene_panel"]
    deployed_perf = pd.read_csv(ML_DIR / "ml_model_performance_external_validation_summary.csv")
    deployed_elastic = deployed_perf[deployed_perf["model"] == "elastic_net"].iloc[0]

    comparison = pd.DataFrame(
        [
            {
                "comparison_metric": "mean_external_roc_auc",
                "observed_25_gene_panel_refit": observed["mean_external_roc_auc"],
                "random_mean": random_df["mean_external_roc_auc"].mean(),
                "random_sd": random_df["mean_external_roc_auc"].std(),
                "random_95th_percentile": random_df["mean_external_roc_auc"].quantile(0.95),
                "random_max": random_df["mean_external_roc_auc"].max(),
                "empirical_p_random_ge_observed": (1 + (random_df["mean_external_roc_auc"] >= observed["mean_external_roc_auc"]).sum())
                / (len(random_df) + 1),
                "deployed_elastic_net_mean_external_roc_auc": deployed_elastic["mean_external_roc_auc"],
            },
            {
                "comparison_metric": "min_external_roc_auc",
                "observed_25_gene_panel_refit": observed["min_external_roc_auc"],
                "random_mean": random_df["min_external_roc_auc"].mean(),
                "random_sd": random_df["min_external_roc_auc"].std(),
                "random_95th_percentile": random_df["min_external_roc_auc"].quantile(0.95),
                "random_max": random_df["min_external_roc_auc"].max(),
                "empirical_p_random_ge_observed": (1 + (random_df["min_external_roc_auc"] >= observed["min_external_roc_auc"]).sum())
                / (len(random_df) + 1),
                "deployed_elastic_net_mean_external_roc_auc": deployed_elastic["min_external_roc_auc"],
            },
        ]
    )

    summary.to_csv(OUT_DIR / "random_robust_25_gene_panel_baseline.csv", index=False, encoding="utf-8-sig")
    cohort_perf.to_csv(OUT_DIR / "random_robust_25_gene_panel_cohort_performance.csv", index=False, encoding="utf-8-sig")
    comparison.to_csv(OUT_DIR / "random_robust_25_gene_panel_baseline_summary.csv", index=False, encoding="utf-8-sig")
    return summary, cohort_perf, comparison


def plot_outputs(per_cohort: pd.DataFrame, loo: pd.DataFrame, random_summary: pd.DataFrame, random_comparison: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.bar(per_cohort["series_id"], per_cohort["roc_auc"], color="#4f7ea8")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("External ROC AUC")
    ax.set_title("Elastic Net performance by validation cohort")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "elastic_net_per_cohort_auc.png", dpi=220)
    fig.savefig(PLOT_DIR / "elastic_net_per_cohort_auc.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    labels = [f"omit {x}" for x in loo["omitted_validation_cohort"]]
    ax.bar(labels, loo["mean_cohort_roc_auc"], color="#7a9a5b")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Mean cohort ROC AUC")
    ax.set_title("Leave-one-validation-cohort sensitivity")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "leave_one_validation_cohort_sensitivity.png", dpi=220)
    fig.savefig(PLOT_DIR / "leave_one_validation_cohort_sensitivity.pdf")
    plt.close(fig)

    random_df = random_summary[random_summary["panel_type"] == "random_robust_25_gene_panel"]
    observed = random_summary[random_summary["panel_type"] == "observed_25_gene_panel_refit"].iloc[0]
    deployed = random_comparison.loc[
        random_comparison["comparison_metric"] == "mean_external_roc_auc",
        "deployed_elastic_net_mean_external_roc_auc",
    ].iloc[0]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.hist(random_df["mean_external_roc_auc"], bins=30, color="#b8c7d9", edgecolor="white")
    ax.axvline(observed["mean_external_roc_auc"], color="#b42318", linewidth=2, label="Observed 25-gene refit")
    ax.axvline(deployed, color="#1f4e79", linewidth=2, linestyle="--", label="Deployed Elastic Net")
    ax.set_xlabel("Mean external ROC AUC")
    ax.set_ylabel("Random panel count")
    ax.set_title("Random robust 25-gene panel baseline")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "random_robust_25_gene_panel_auc_distribution.png", dpi=220)
    fig.savefig(PLOT_DIR / "random_robust_25_gene_panel_auc_distribution.pdf")
    plt.close(fig)


def write_excel() -> None:
    with pd.ExcelWriter(OUT_DIR / "machine_learning_sensitivity_summary.xlsx", engine="openpyxl") as writer:
        for csv_name in [
            "elastic_net_external_per_cohort_performance",
            "leave_one_validation_cohort_sensitivity",
            "random_robust_25_gene_panel_baseline_summary",
            "random_robust_25_gene_panel_baseline",
            "random_robust_25_gene_panel_cohort_performance",
        ]:
            pd.read_csv(OUT_DIR / f"{csv_name}.csv").to_excel(writer, sheet_name=csv_name[:31], index=False)


def qc_outputs(per_cohort: pd.DataFrame, loo: pd.DataFrame, random_summary: pd.DataFrame, random_comparison: pd.DataFrame) -> None:
    random_df = random_summary[random_summary["panel_type"] == "random_robust_25_gene_panel"]
    observed_rows = random_summary[random_summary["panel_type"] == "observed_25_gene_panel_refit"]
    qc = pd.DataFrame(
        [
            {
                "check": "validation_cohorts_present",
                "pass": set(per_cohort["series_id"]) == set(VALIDATION_DATASETS),
                "value": ";".join(per_cohort["series_id"].astype(str)),
            },
            {
                "check": "leave_one_rows",
                "pass": len(loo) == len(VALIDATION_DATASETS),
                "value": len(loo),
            },
            {
                "check": "random_panel_count",
                "pass": len(random_df) == N_RANDOM_PANELS,
                "value": len(random_df),
            },
            {
                "check": "observed_panel_present",
                "pass": len(observed_rows) == 1,
                "value": len(observed_rows),
            },
            {
                "check": "auc_values_in_range",
                "pass": random_summary["mean_external_roc_auc"].between(0, 1).all() and per_cohort["roc_auc"].between(0, 1).all(),
                "value": "0_to_1",
            },
            {
                "check": "empirical_p_computed",
                "pass": random_comparison["empirical_p_random_ge_observed"].between(0, 1).all(),
                "value": ";".join(random_comparison["empirical_p_random_ge_observed"].round(4).astype(str)),
            },
        ]
    )
    qc.to_csv(OUT_DIR / "machine_learning_sensitivity_qc.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    per_cohort, loo = leave_one_validation_sensitivity()
    random_summary, random_cohort, random_comparison = random_robust_gene_panel_baseline()
    plot_outputs(per_cohort, loo, random_summary, random_comparison)
    write_excel()
    qc_outputs(per_cohort, loo, random_summary, random_comparison)
    print(OUT_DIR)
    print(random_comparison.to_string(index=False))


if __name__ == "__main__":
    main()
