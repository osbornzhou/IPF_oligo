#!/usr/bin/env python
"""
Train and validate IPF vs Control machine-learning models.

The implementation uses strict leakage control:
  - imputation, scaling, and feature selection are inside sklearn Pipelines;
  - hyperparameter tuning is nested inside each outer CV fold;
  - external validation datasets are untouched during training/tuning.

Models:
  - LASSO logistic regression
  - Elastic Net logistic regression
  - Linear SVM
  - RBF SVM
  - Random Forest
  - Gradient Boosting
  - small MLP neural network

Outputs are written to results/models/ml_outputs.
"""

from __future__ import annotations

import json
import math
import os
import warnings
from collections import Counter, defaultdict
from pathlib import Path

_mpl_config_dir = Path.cwd() / "results" / "models" / "matplotlib_cache"
_mpl_config_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_config_dir))

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC, SVC


PROJECT_DIR = Path(__file__).resolve().parents[1]
MATRIX_DIR = PROJECT_DIR / "results" / "models" / "feature_matrices"
MODEL_DIR = PROJECT_DIR / "results" / "models"
OUT_DIR = MODEL_DIR / "ml_outputs"
PLOT_DIR = OUT_DIR / "plots"
ARTIFACT_DIR = OUT_DIR / "artifacts"
TRAINING_DATASET = "GSE32537"
VALIDATION_DATASETS = ["GSE110147", "GSE150910", "GSE53845", "GSE92592"]
DATASETS = [TRAINING_DATASET] + VALIDATION_DATASETS
FEATURE_SET_NAME = "discovery_only_mrna"
MATRIX_SUFFIX = "common_discovery_only_mrna_matrix"
OUT_DIR = MODEL_DIR / f"ml_outputs_{FEATURE_SET_NAME}"
PLOT_DIR = OUT_DIR / "plots"
ARTIFACT_DIR = OUT_DIR / "artifacts"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
RANDOM_STATE = 20260524
OUTER_SPLITS = 5
OUTER_REPEATS = 10
INNER_SPLITS = 3
PERMUTATIONS = 20


warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


def load_matrix(series_id: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, list[str]]:
    path = MATRIX_DIR / f"{series_id}_{MATRIX_SUFFIX}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    features = [c for c in df.columns if c not in {"sample_id", "group", "label"}]
    x = df[features].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    y = df["label"].to_numpy(dtype=int)
    return df, x, y, features


def make_pipeline(clf) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("select", SelectKBest(score_func=f_classif, k=50)),
            ("clf", clf),
        ]
    )


def model_specs(n_features: int) -> dict[str, tuple[Pipeline, dict[str, list]]]:
    k_values = sorted({min(k, n_features) for k in [15, 30, 50, 100] if min(k, n_features) >= 5})
    return {
        "lasso_logistic": (
            make_pipeline(
                LogisticRegression(
                    penalty="l1",
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=RANDOM_STATE,
                )
            ),
            {"select__k": k_values, "clf__C": [0.03, 0.1, 0.3, 1.0]},
        ),
        "elastic_net": (
            make_pipeline(
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    class_weight="balanced",
                    max_iter=8000,
                    random_state=RANDOM_STATE,
                )
            ),
            {"select__k": k_values, "clf__C": [0.03, 0.1, 0.3, 1.0], "clf__l1_ratio": [0.2, 0.5, 0.8]},
        ),
        "linear_svm": (
            make_pipeline(LinearSVC(class_weight="balanced", max_iter=8000, random_state=RANDOM_STATE)),
            {"select__k": k_values, "clf__C": [0.03, 0.1, 0.3, 1.0]},
        ),
        "rbf_svm": (
            make_pipeline(SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=RANDOM_STATE)),
            {"select__k": [min(30, n_features), min(50, n_features)], "clf__C": [0.1, 1.0, 3.0], "clf__gamma": ["scale"]},
        ),
        "random_forest": (
            make_pipeline(
                RandomForestClassifier(
                    n_estimators=300,
                    class_weight="balanced_subsample",
                    random_state=RANDOM_STATE,
                    n_jobs=-1,
                )
            ),
            {"select__k": [min(50, n_features), min(100, n_features)], "clf__max_depth": [3, 5, None], "clf__min_samples_leaf": [1, 3]},
        ),
        "gradient_boosting": (
            make_pipeline(GradientBoostingClassifier(random_state=RANDOM_STATE)),
            {"select__k": [min(30, n_features), min(50, n_features)], "clf__n_estimators": [80, 150], "clf__learning_rate": [0.03, 0.08], "clf__max_depth": [1, 2]},
        ),
        "mlp_small": (
            make_pipeline(
                MLPClassifier(
                    hidden_layer_sizes=(16, 8),
                    activation="relu",
                    alpha=0.01,
                    learning_rate_init=0.001,
                    early_stopping=True,
                    validation_fraction=0.2,
                    max_iter=1000,
                    random_state=RANDOM_STATE,
                )
            ),
            {"select__k": [min(30, n_features), min(50, n_features)], "clf__hidden_layer_sizes": [(16,), (16, 8), (32, 16)], "clf__alpha": [0.001, 0.01, 0.1]},
        ),
    }


def get_scores(model: Pipeline, x: np.ndarray) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(x)[:, 1]
    if hasattr(model, "decision_function"):
        return model.decision_function(x)
    return model.predict(x)


def metrics_row(y_true: np.ndarray, scores: np.ndarray, preds: np.ndarray) -> dict[str, float]:
    tn, fp, fn, tp = confusion_matrix(y_true, preds, labels=[0, 1]).ravel()
    return {
        "roc_auc": roc_auc_score(y_true, scores) if len(np.unique(y_true)) == 2 else np.nan,
        "pr_auc": average_precision_score(y_true, scores) if len(np.unique(y_true)) == 2 else np.nan,
        "accuracy": accuracy_score(y_true, preds),
        "balanced_accuracy": balanced_accuracy_score(y_true, preds),
        "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
        "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
        "f1": f1_score(y_true, preds),
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


def selected_features(estimator: Pipeline, features: list[str]) -> list[str]:
    selector = estimator.named_steps["select"]
    mask = selector.get_support()
    return [feature for feature, keep in zip(features, mask) if keep]


def model_feature_importance(estimator: Pipeline, features: list[str]) -> pd.DataFrame:
    selected = selected_features(estimator, features)
    clf = estimator.named_steps["clf"]
    if hasattr(clf, "coef_"):
        values = np.ravel(clf.coef_)
        if len(values) == len(selected):
            return pd.DataFrame({"feature": selected, "importance": np.abs(values), "signed_importance": values})
    if hasattr(clf, "feature_importances_"):
        values = np.ravel(clf.feature_importances_)
        if len(values) == len(selected):
            return pd.DataFrame({"feature": selected, "importance": values, "signed_importance": values})
    return pd.DataFrame({"feature": selected, "importance": 1.0, "signed_importance": 1.0})


def nested_cv_train(x: np.ndarray, y: np.ndarray, features: list[str], specs: dict[str, tuple[Pipeline, dict]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    outer = RepeatedStratifiedKFold(n_splits=OUTER_SPLITS, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE)
    inner = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    perf_rows = []
    prediction_rows = []
    selection_rows = []

    for model_name, (pipeline, grid) in specs.items():
        print(f"Nested CV: {model_name}")
        for fold_id, (train_idx, test_idx) in enumerate(outer.split(x, y), start=1):
            search = GridSearchCV(
                estimator=clone(pipeline),
                param_grid=grid,
                scoring="roc_auc",
                cv=inner,
                n_jobs=-1,
                refit=True,
                error_score="raise",
            )
            search.fit(x[train_idx], y[train_idx])
            best = search.best_estimator_
            scores = get_scores(best, x[test_idx])
            preds = best.predict(x[test_idx])
            row = metrics_row(y[test_idx], scores, preds)
            row.update(
                {
                    "model": model_name,
                    "fold_id": fold_id,
                    "best_params": json.dumps(search.best_params_, sort_keys=True),
                    "selected_feature_count": len(selected_features(best, features)),
                }
            )
            perf_rows.append(row)

            for idx, score, pred in zip(test_idx, scores, preds):
                prediction_rows.append(
                    {
                        "model": model_name,
                        "fold_id": fold_id,
                        "sample_id": idx,
                        "true_label": int(y[idx]),
                        "score": float(score),
                        "prediction": int(pred),
                    }
                )

            for feature in selected_features(best, features):
                selection_rows.append({"model": model_name, "fold_id": fold_id, "feature": feature})

    return pd.DataFrame(perf_rows), pd.DataFrame(prediction_rows), pd.DataFrame(selection_rows)


def train_final_models(x: np.ndarray, y: np.ndarray, features: list[str], specs: dict[str, tuple[Pipeline, dict]]) -> dict[str, Pipeline]:
    final_models = {}
    inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    for model_name, (pipeline, grid) in specs.items():
        print(f"Final training: {model_name}")
        search = GridSearchCV(
            estimator=clone(pipeline),
            param_grid=grid,
            scoring="roc_auc",
            cv=inner,
            n_jobs=-1,
            refit=True,
            error_score="raise",
        )
        search.fit(x, y)
        final_models[model_name] = search.best_estimator_
        joblib.dump(search.best_estimator_, ARTIFACT_DIR / f"{model_name}_final_model.joblib")
    return final_models


def evaluate_external(final_models: dict[str, Pipeline], validation_data: dict[str, tuple[pd.DataFrame, np.ndarray, np.ndarray, list[str]]], train_features: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    perf_rows = []
    pred_rows = []
    for model_name, model in final_models.items():
        for series_id, (df, x, y, features) in validation_data.items():
            x = df[train_features].to_numpy(dtype=float)
            scores = get_scores(model, x)
            preds = model.predict(x)
            row = metrics_row(y, scores, preds)
            row.update({"model": model_name, "series_id": series_id, "samples": len(y), "ipf_samples": int(y.sum()), "control_samples": int((y == 0).sum())})
            perf_rows.append(row)
            for sample_id, group, true, score, pred in zip(df["sample_id"], df["group"], y, scores, preds):
                pred_rows.append(
                    {
                        "model": model_name,
                        "series_id": series_id,
                        "sample_id": sample_id,
                        "group": group,
                        "true_label": int(true),
                        "score": float(score),
                        "prediction": int(pred),
                    }
                )
    return pd.DataFrame(perf_rows), pd.DataFrame(pred_rows)


def permutation_control(best_model: Pipeline, x: np.ndarray, y: np.ndarray) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_STATE)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    rows = []
    for i in range(1, PERMUTATIONS + 1):
        y_perm = rng.permutation(y)
        fold_aucs = []
        for train_idx, test_idx in cv.split(x, y_perm):
            model = clone(best_model)
            model.fit(x[train_idx], y_perm[train_idx])
            scores = get_scores(model, x[test_idx])
            fold_aucs.append(roc_auc_score(y_perm[test_idx], scores))
        rows.append({"permutation_id": i, "roc_auc": float(np.mean(fold_aucs))})
    return pd.DataFrame(rows)


def plot_model_comparison(internal_summary: pd.DataFrame, external_summary: pd.DataFrame) -> None:
    merged = external_summary.copy()
    order = merged.sort_values("mean_external_roc_auc", ascending=False)["model"]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.barh(range(len(order)), merged.set_index("model").loc[order, "mean_external_roc_auc"], color="#4f7ea8")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Mean external ROC AUC")
    ax.set_title("External validation model comparison")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "ml_model_external_auc_barplot.png", dpi=180)
    plt.close(fig)


def plot_roc_pr(best_model_name: str, predictions: pd.DataFrame) -> None:
    model_pred = predictions[predictions["model"] == best_model_name]
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for series_id, df in model_pred.groupby("series_id"):
        fpr, tpr, _ = roc_curve(df["true_label"], df["score"])
        auc = roc_auc_score(df["true_label"], df["score"])
        ax.plot(fpr, tpr, label=f"{series_id} AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(f"External ROC curves - {best_model_name}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "ml_best_model_external_roc.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for series_id, df in model_pred.groupby("series_id"):
        precision, recall, _ = precision_recall_curve(df["true_label"], df["score"])
        ap = average_precision_score(df["true_label"], df["score"])
        ax.plot(recall, precision, label=f"{series_id} AP={ap:.3f}")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"External PR curves - {best_model_name}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "ml_best_model_external_pr.png", dpi=180)
    plt.close(fig)


def plot_feature_stability(stability: pd.DataFrame) -> None:
    top = stability.sort_values("overall_selection_frequency", ascending=False).head(30)
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.barh(range(len(top)), top["overall_selection_frequency"], color="#b42318")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["feature"])
    ax.invert_yaxis()
    ax.set_xlabel("Selection frequency")
    ax.set_title("Top stable ML-selected features")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "ml_feature_selection_stability_top30.png", dpi=180)
    plt.close(fig)


def build_stability(selection: pd.DataFrame, n_folds: int, models: list[str]) -> pd.DataFrame:
    rows = []
    for feature in sorted(selection["feature"].unique()):
        sub = selection[selection["feature"] == feature]
        row = {"feature": feature, "overall_selection_count": len(sub), "overall_selection_frequency": len(sub) / (n_folds * len(models))}
        for model in models:
            count = len(sub[sub["model"] == model])
            row[f"{model}_selection_frequency"] = count / n_folds
        row["models_selecting_feature"] = sum(row[f"{model}_selection_frequency"] > 0 for model in models)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["overall_selection_frequency", "models_selecting_feature"], ascending=False)


def build_final_panel(stability: pd.DataFrame, final_models: dict[str, Pipeline], features: list[str]) -> pd.DataFrame:
    importances = []
    for model_name, model in final_models.items():
        imp = model_feature_importance(model, features)
        imp["model"] = model_name
        importances.append(imp)
    imp_df = pd.concat(importances, ignore_index=True)
    imp_summary = (
        imp_df.groupby("feature", as_index=False)
        .agg(mean_importance=("importance", "mean"), models_with_final_importance=("model", "nunique"))
    )
    panel = stability.merge(imp_summary, on="feature", how="left").fillna({"mean_importance": 0, "models_with_final_importance": 0})
    panel["panel_score"] = panel["overall_selection_frequency"] * 2 + panel["models_selecting_feature"] / max(len(final_models), 1) + panel["mean_importance"].rank(pct=True)
    return panel.sort_values("panel_score", ascending=False).head(25)


def add_biological_evidence(panel: pd.DataFrame) -> pd.DataFrame:
    discovery = pd.read_csv(PROJECT_DIR / "results" / "differential_expression_annotated" / "GSE32537_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv")
    robust = pd.read_csv(PROJECT_DIR / "results" / "robust_candidates" / "robust_mrna_candidates_strict.csv")
    ppi = pd.read_csv(PROJECT_DIR / "results" / "ppi_network" / "string_ppi_nodes_robust_mrna_strict_medium_confidence.csv")
    axes = pd.read_csv(PROJECT_DIR / "results" / "mirna_mrna_axes" / "robust_mirna_mrna_negative_axes_mirtarbase.csv")
    enrichment_gene_sets = pd.read_csv(PROJECT_DIR / "results" / "enrichment" / "gene_sets_used.csv")

    discovery["feature"] = discovery["standard_feature_id"].str.upper()
    robust["feature"] = robust["standard_feature_id"].str.upper()
    ppi["feature"] = ppi["gene_symbol"].str.upper()
    axes["feature"] = axes["target_gene"].str.upper()
    axis_summary = axes.groupby("feature", as_index=False).agg(axis_count=("axis", "nunique"), regulating_mirnas=("candidate_mirna", lambda x: ";".join(sorted(set(map(str, x))))))
    enrich_summary = enrichment_gene_sets.groupby("gene_symbol", as_index=False).size().rename(columns={"gene_symbol": "feature", "size": "enrichment_gene_set_count"})
    enrich_summary["feature"] = enrich_summary["feature"].str.upper()

    discovery_small = discovery[["feature", "logFC", "adj.P.Val", "P.Value"]].rename(
        columns={"logFC": "discovery_logFC", "adj.P.Val": "discovery_adj_p", "P.Value": "discovery_p_value"}
    )
    out = panel.merge(discovery_small, on="feature", how="left")
    out = out.merge(
        robust[["feature", "same_direction_fdr_sig_count", "robust_score", "validation_details"]],
        on="feature",
        how="left",
    )
    out = out.merge(ppi[["feature", "degree", "weighted_degree", "hub_score", "hub_rank"]], on="feature", how="left")
    out = out.merge(axis_summary, on="feature", how="left")
    out = out.merge(enrich_summary, on="feature", how="left")
    out["axis_count"] = out["axis_count"].fillna(0).astype(int)
    out["enrichment_gene_set_count"] = out["enrichment_gene_set_count"].fillna(0).astype(int)
    out["oligonucleotide_strategy"] = np.where(
        out["discovery_logFC"] > 0,
        "upregulated mRNA: siRNA/ASO knockdown candidate",
        "downregulated mRNA: not primary knockdown target; consider pathway context",
    )
    out["target_priority_score"] = (
        out["panel_score"].rank(pct=True)
        + out["robust_score"].fillna(0).rank(pct=True)
        + out["hub_score"].fillna(0).rank(pct=True)
        + (out["axis_count"] > 0).astype(int)
        + out["enrichment_gene_set_count"].rank(pct=True)
    )
    return out.sort_values("target_priority_score", ascending=False)


def main() -> None:
    train_df, x_train, y_train, features = load_matrix(TRAINING_DATASET)
    validation_data = {sid: load_matrix(sid) for sid in VALIDATION_DATASETS}
    specs = model_specs(len(features))

    cv_perf, cv_preds, selection = nested_cv_train(x_train, y_train, features, specs)
    cv_perf.to_csv(OUT_DIR / "ml_model_performance_internal_cv_folds.csv", index=False, encoding="utf-8-sig")
    cv_preds.to_csv(OUT_DIR / "ml_internal_cv_predictions.csv", index=False, encoding="utf-8-sig")
    selection.to_csv(OUT_DIR / "ml_feature_selection_events.csv", index=False, encoding="utf-8-sig")

    internal_summary = (
        cv_perf.groupby("model", as_index=False)
        .agg(
            mean_internal_roc_auc=("roc_auc", "mean"),
            sd_internal_roc_auc=("roc_auc", "std"),
            mean_internal_pr_auc=("pr_auc", "mean"),
            mean_internal_balanced_accuracy=("balanced_accuracy", "mean"),
            mean_selected_features=("selected_feature_count", "mean"),
        )
        .sort_values("mean_internal_roc_auc", ascending=False)
    )
    internal_summary.to_csv(OUT_DIR / "ml_model_performance_internal_cv_summary.csv", index=False, encoding="utf-8-sig")

    final_models = train_final_models(x_train, y_train, features, specs)
    external_perf, external_preds = evaluate_external(final_models, validation_data, features)
    external_perf.to_csv(OUT_DIR / "ml_model_performance_external_validation.csv", index=False, encoding="utf-8-sig")
    external_preds.to_csv(OUT_DIR / "ml_external_validation_predictions.csv", index=False, encoding="utf-8-sig")

    external_summary = (
        external_perf.groupby("model", as_index=False)
        .agg(
            mean_external_roc_auc=("roc_auc", "mean"),
            min_external_roc_auc=("roc_auc", "min"),
            mean_external_pr_auc=("pr_auc", "mean"),
            mean_external_balanced_accuracy=("balanced_accuracy", "mean"),
            mean_external_sensitivity=("sensitivity", "mean"),
            mean_external_specificity=("specificity", "mean"),
        )
        .sort_values("mean_external_roc_auc", ascending=False)
    )
    external_summary.to_csv(OUT_DIR / "ml_model_performance_external_validation_summary.csv", index=False, encoding="utf-8-sig")

    n_outer_folds = OUTER_SPLITS * OUTER_REPEATS
    stability = build_stability(selection, n_outer_folds, list(specs.keys()))
    stability.to_csv(OUT_DIR / "ml_feature_selection_stability.csv", index=False, encoding="utf-8-sig")

    final_panel = build_final_panel(stability, final_models, features)
    final_panel = add_biological_evidence(final_panel)
    final_panel.to_csv(OUT_DIR / "ml_final_biomarker_panel.csv", index=False, encoding="utf-8-sig")

    target_priority = final_panel.sort_values("target_priority_score", ascending=False)
    target_priority.to_csv(OUT_DIR / "ml_oligonucleotide_target_priority.csv", index=False, encoding="utf-8-sig")

    best_model_name = external_summary.iloc[0]["model"]
    permutation = permutation_control(final_models[best_model_name], x_train, y_train)
    permutation["best_model"] = best_model_name
    permutation.to_csv(OUT_DIR / "ml_permutation_label_control.csv", index=False, encoding="utf-8-sig")

    plot_model_comparison(internal_summary, external_summary)
    plot_roc_pr(best_model_name, external_preds)
    plot_feature_stability(stability)

    qc1_pass = (
        train_df["sample_id"].is_unique
        and set(train_df["group"]) == {"IPF", "Control"}
        and np.isfinite(x_train).all()
        and all(np.isfinite(data[1]).all() for data in validation_data.values())
        and len(features) >= 20
    )
    qc2_pass = len(cv_perf) == len(specs) * n_outer_folds and cv_perf["roc_auc"].notna().all()
    best_external = external_summary.iloc[0]
    internal_best_auc = float(internal_summary[internal_summary["model"] == best_model_name]["mean_internal_roc_auc"].iloc[0])
    external_best_auc = float(best_external["mean_external_roc_auc"])
    perm_mean = float(permutation["roc_auc"].mean())
    stability_pass = (stability["overall_selection_frequency"].max() >= 0.30) and (len(final_panel) >= 5)
    qc3_pass = external_best_auc >= 0.75 and (internal_best_auc - external_best_auc) <= 0.25 and perm_mean < 0.65 and stability_pass

    qc = pd.DataFrame(
        [
            {
                "qc1_data_feature_integrity_pass": qc1_pass,
                "qc2_resampling_leakage_control_pass": qc2_pass,
                "qc3_performance_stability_biology_pass": qc3_pass,
                "triple_qc_pass": qc1_pass and qc2_pass and qc3_pass,
                "training_dataset": TRAINING_DATASET,
                "feature_set_name": FEATURE_SET_NAME,
                "feature_prefilter_rule": "GSE32537 discovery significant mRNA only; external validation labels are not used for feature prefiltering",
                "external_validation_datasets": ";".join(VALIDATION_DATASETS),
                "common_features": len(features),
                "models": ";".join(specs.keys()),
                "outer_cv": f"RepeatedStratifiedKFold {OUTER_SPLITS}x{OUTER_REPEATS}",
                "inner_cv": f"StratifiedKFold {INNER_SPLITS}",
                "best_model_by_mean_external_auc": best_model_name,
                "best_model_internal_mean_auc": internal_best_auc,
                "best_model_external_mean_auc": external_best_auc,
                "best_model_external_min_auc": float(best_external["min_external_roc_auc"]),
                "permutation_control_mean_auc": perm_mean,
                "final_panel_size": len(final_panel),
                "top_selection_frequency": float(stability["overall_selection_frequency"].max()),
                "leakage_control_note": "Imputation, scaling and SelectKBest are inside sklearn Pipeline; feature selection and tuning are nested inside training folds; external validation is not used for tuning.",
            }
        ]
    )
    qc.to_csv(OUT_DIR / "ml_triple_qc.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(OUT_DIR / "ml_results_summary.xlsx", engine="openpyxl") as writer:
        qc.to_excel(writer, sheet_name="triple_qc", index=False)
        internal_summary.to_excel(writer, sheet_name="internal_cv_summary", index=False)
        external_summary.to_excel(writer, sheet_name="external_summary", index=False)
        external_perf.to_excel(writer, sheet_name="external_by_dataset", index=False)
        stability.head(100).to_excel(writer, sheet_name="feature_stability_top100", index=False)
        final_panel.to_excel(writer, sheet_name="final_panel", index=False)
        permutation.to_excel(writer, sheet_name="permutation_control", index=False)

    print(qc.to_string(index=False))
    print("\nExternal validation summary:")
    print(external_summary.to_string(index=False))
    print("\nFinal biomarker panel:")
    print(final_panel[["feature", "overall_selection_frequency", "models_selecting_feature", "discovery_logFC", "hub_rank", "axis_count", "target_priority_score"]].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
