#!/usr/bin/env python
"""
Add targeted boundary analyses for the submission manuscript.

This script deliberately avoids adding a new classifier family. It uses
regularized final-panel refits and target-set tests to clarify two manuscript
boundaries: high external AUC values and the narrow exact-miRNA result layer.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


PROJECT_DIR = Path(__file__).resolve().parents[1]
MODEL_DIR = PROJECT_DIR / "results" / "models"
FEATURE_DIR = MODEL_DIR / "feature_matrices"
ML_DIR = MODEL_DIR / "ml_outputs_discovery_only_mrna"
ROBUST_DIR = PROJECT_DIR / "results" / "robust_candidates"
DE_DIR = PROJECT_DIR / "results" / "differential_expression_annotated"
META_DIR = PROJECT_DIR / "metadata"
FIG_DIR = PROJECT_DIR / "manuscript" / "figures"
AUC_OUT = PROJECT_DIR / "results" / "ml_stress_tests"
MIRNA_OUT = PROJECT_DIR / "results" / "mirna_program_support"
MIRTARBASE_PATH = PROJECT_DIR / "data_external" / "miRTarBase" / "hsa_MTI_miRTarBase_2025_v10.csv"
EXPR_DIR = PROJECT_DIR / "data_processed" / "expression"
FEATURE_MAP_PATH = PROJECT_DIR / "results" / "feature_annotation" / "feature_annotation_map.csv"

for path in [FIG_DIR, AUC_OUT, MIRNA_OUT]:
    path.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260527)


def normalize_gene(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().upper()
    return "" if text in {"", "NA", "NAN"} else text


def normalize_mirna(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return text


def mirna_base(value: object) -> str:
    text = normalize_mirna(value).lower()
    for suffix in ["-3p", "-5p"]:
        if text.endswith(suffix):
            return text[: -len(suffix)]
    return text


def sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.clip(z, -35, 35)
    return 1.0 / (1.0 + np.exp(-z))


def roc_auc_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    pos = y == 1
    neg = y == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(score)
    ranks = np.empty_like(order, dtype=float)
    sorted_scores = score[order]
    start = 0
    rank_values = np.arange(1, len(score) + 1, dtype=float)
    while start < len(score):
        end = start + 1
        while end < len(score) and sorted_scores[end] == sorted_scores[start]:
            end += 1
        ranks[order[start:end]] = rank_values[start:end].mean()
        start = end
    rank_sum_pos = ranks[pos].sum()
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def average_precision_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-score)
    y_sorted = y[order]
    tp = np.cumsum(y_sorted == 1)
    fp = np.cumsum(y_sorted == 0)
    precision = tp / np.maximum(tp + fp, 1)
    return float((precision * (y_sorted == 1)).sum() / y.sum())


def balanced_accuracy(y: np.ndarray, pred: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    pred = np.asarray(pred).astype(int)
    tp = ((y == 1) & (pred == 1)).sum()
    fn = ((y == 1) & (pred == 0)).sum()
    tn = ((y == 0) & (pred == 0)).sum()
    fp = ((y == 0) & (pred == 1)).sum()
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    return float((sens + spec) / 2.0)


def bh_adjust(pvalues: list[float]) -> list[float]:
    arr = np.asarray(pvalues, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    valid = np.isfinite(arr)
    vals = arr[valid]
    if len(vals) == 0:
        return out.tolist()
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty_like(vals)
    tmp[order] = np.minimum(adj, 1.0)
    out[valid] = tmp
    return out.tolist()


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def read_feature_matrix(series_id: str, common: bool = True) -> pd.DataFrame:
    suffix = "common_discovery_only_mrna_matrix.csv" if common else "gene_level_matrix.csv"
    path = FEATURE_DIR / f"{series_id}_{suffix}"
    return pd.read_csv(path)


def fit_regularized_logistic(X: np.ndarray, y: np.ndarray, l2: float = 0.08, lr: float = 0.04, n_iter: int = 2500) -> tuple[np.ndarray, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, p = X.shape
    w = np.zeros(p, dtype=float)
    b = math.log((y.mean() + 1e-4) / (1 - y.mean() + 1e-4))
    pos_weight = n / (2.0 * max(y.sum(), 1.0))
    neg_weight = n / (2.0 * max((1.0 - y).sum(), 1.0))
    weights = np.where(y == 1, pos_weight, neg_weight)
    for i in range(n_iter):
        pred = sigmoid(X @ w + b)
        err = (pred - y) * weights
        grad_w = X.T @ err / n + l2 * w
        grad_b = err.mean()
        step = lr / math.sqrt(1.0 + i / 400.0)
        w -= step * grad_w
        b -= step * grad_b
    return w, float(b)


def prepare_train_scaler(train: pd.DataFrame, features: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X = train[features].apply(pd.to_numeric, errors="coerce")
    med = X.median(axis=0).to_numpy(dtype=float)
    X = X.fillna(pd.Series(med, index=features))
    mean = X.mean(axis=0).to_numpy(dtype=float)
    sd = X.std(axis=0, ddof=0).replace(0, 1.0).to_numpy(dtype=float)
    Xs = (X.to_numpy(dtype=float) - mean) / sd
    y = train["label"].to_numpy(dtype=int)
    return Xs, y, med, np.vstack([mean, sd])


def transform_matrix(df: pd.DataFrame, features: list[str], med: np.ndarray, mean_sd: np.ndarray) -> np.ndarray:
    mean, sd = mean_sd
    X = df[features].apply(pd.to_numeric, errors="coerce")
    X = X.fillna(pd.Series(med, index=features))
    return (X.to_numpy(dtype=float) - mean) / sd


def evaluate_panel(features: list[str], validation_sets: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = read_feature_matrix("GSE32537")
    Xs, y, med, mean_sd = prepare_train_scaler(train, features)
    w, b = fit_regularized_logistic(Xs, y)
    rows = []
    preds = []
    for series_id in validation_sets:
        val = read_feature_matrix(series_id)
        Xv = transform_matrix(val, features, med, mean_sd)
        score = sigmoid(Xv @ w + b)
        pred = (score >= 0.5).astype(int)
        yv = val["label"].to_numpy(dtype=int)
        rows.append(
            {
                "series_id": series_id,
                "samples": len(yv),
                "ipf_samples": int(yv.sum()),
                "control_samples": int((yv == 0).sum()),
                "roc_auc": roc_auc_score(yv, score),
                "pr_auc": average_precision_score(yv, score),
                "balanced_accuracy": balanced_accuracy(yv, pred),
            }
        )
        tmp = val[["sample_id", "group", "label"]].copy()
        tmp["series_id"] = series_id
        tmp["score"] = score
        tmp["prediction"] = pred
        preds.append(tmp)
    return pd.DataFrame(rows), pd.concat(preds, ignore_index=True)


def collapse_gse110147_all_samples() -> pd.DataFrame:
    expr = pd.read_csv(EXPR_DIR / "GSE110147_expression_in_annotation_order.csv.gz")
    expr = expr.rename(columns={expr.columns[0]: "feature_id"})
    expr["feature_id"] = expr["feature_id"].astype(str)
    fmap = pd.read_csv(FEATURE_MAP_PATH)
    mapping = fmap[(fmap["series_id"] == "GSE110147") & (fmap["is_annotated"].astype(str).str.lower().isin(["true", "1", "yes"]))].copy()
    mapping["feature_id"] = mapping["feature_id"].astype(str)
    mapping["standard_feature_id"] = mapping["standard_feature_id"].map(normalize_gene)
    mapping = mapping[mapping["standard_feature_id"] != ""].drop_duplicates("feature_id", keep="first")
    merged = expr.merge(mapping[["feature_id", "standard_feature_id"]], on="feature_id", how="inner")
    sample_cols = [c for c in merged.columns if c not in {"feature_id", "standard_feature_id"}]
    merged[sample_cols] = merged[sample_cols].apply(pd.to_numeric, errors="coerce")
    merged["row_variance"] = merged[sample_cols].var(axis=1, skipna=True)
    merged = merged.sort_values(["standard_feature_id", "row_variance"], ascending=[True, False])
    collapsed = merged.drop_duplicates("standard_feature_id", keep="first")
    out = collapsed[["standard_feature_id"] + sample_cols].set_index("standard_feature_id").T
    out.index.name = "sample_id"
    return out.reset_index()


def score_gse110147_disease_controls(panel: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = read_feature_matrix("GSE32537")
    Xs, y, med, mean_sd = prepare_train_scaler(train, panel)
    w, b = fit_regularized_logistic(Xs, y)
    all_gse = collapse_gse110147_all_samples()
    audit = pd.read_csv(META_DIR / "all_bulk_mirna_sample_label_audit.csv")
    audit = audit[audit["series_id"].eq("GSE110147")].copy()
    audit["stress_group"] = "Other excluded"
    audit.loc[audit["final_analysis_label"].eq("Control"), "stress_group"] = "Normal control"
    audit.loc[audit["final_analysis_label"].eq("IPF"), "stress_group"] = "IPF"
    audit.loc[audit["original_disease_label"].str.contains("subgroup=NSIP", na=False), "stress_group"] = "NSIP"
    audit.loc[audit["original_disease_label"].str.contains("subgroup=IPF-NSIP", na=False), "stress_group"] = "Mixed IPF-NSIP"
    samples = audit[["sample_id", "stress_group", "included_excluded", "original_disease_label"]].merge(all_gse, on="sample_id", how="inner")
    X = transform_matrix(samples, panel, med, mean_sd)
    samples["refit_final_panel_score"] = sigmoid(X @ w + b)
    order = ["Normal control", "NSIP", "Mixed IPF-NSIP", "IPF"]
    rows = []
    for group in order:
        vals = samples.loc[samples["stress_group"].eq(group), "refit_final_panel_score"].dropna().to_numpy(dtype=float)
        if len(vals) == 0:
            continue
        rows.append(
            {
                "stress_group": group,
                "samples": len(vals),
                "mean_score": float(np.mean(vals)),
                "median_score": float(np.median(vals)),
                "q1_score": float(np.quantile(vals, 0.25)),
                "q3_score": float(np.quantile(vals, 0.75)),
                "min_score": float(np.min(vals)),
                "max_score": float(np.max(vals)),
                "interpretation": "final-panel refit disease-state score; not an IPF-specific diagnostic test",
            }
        )
    return samples[["sample_id", "stress_group", "included_excluded", "original_disease_label", "refit_final_panel_score"]], pd.DataFrame(rows)


def matched_random_panel_baseline(panel: list[str], n_iter: int = 500) -> tuple[pd.DataFrame, pd.DataFrame]:
    common_features = set(read_feature_matrix("GSE32537").columns) - {"sample_id", "group", "label"}
    de = pd.read_csv(DE_DIR / "GSE32537_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv")
    de["gene"] = de["standard_feature_id"].map(normalize_gene)
    de = de[de["gene"].isin(common_features)].drop_duplicates("gene")
    support = pd.read_csv(ROBUST_DIR / "mrna_discovery_validation_summary.csv")
    support["gene"] = support["standard_feature_id"].map(normalize_gene)
    support = support.drop_duplicates("gene")
    candidate_df = de[["gene", "logFC", "adj.P.Val"]].rename(columns={"logFC": "discovery_logFC", "adj.P.Val": "discovery_adj_p"})
    candidate_df = candidate_df.merge(
        support[["gene", "same_direction_fdr_sig_count"]],
        on="gene",
        how="left",
    )
    candidate_df["same_direction_fdr_sig_count"] = candidate_df["same_direction_fdr_sig_count"].fillna(0)
    final_panel = pd.read_csv(ML_DIR / "ml_final_biomarker_panel.csv")
    final_panel["gene"] = final_panel["feature"].map(normalize_gene)
    observed = final_panel[final_panel["gene"].isin(panel)].copy()
    observed = observed[["gene", "discovery_logFC", "discovery_adj_p", "same_direction_fdr_sig_count"]].drop_duplicates("gene")
    if observed["gene"].nunique() != len(panel):
        missing = sorted(set(panel) - set(observed["gene"]))
        raise ValueError(f"Panel features missing from final panel table: {missing}")
    missing_from_candidates = sorted(set(panel) - set(candidate_df["gene"]))
    if missing_from_candidates:
        candidate_df = pd.concat(
            [
                candidate_df,
                observed[observed["gene"].isin(missing_from_candidates)].rename(columns={"same_direction_fdr_sig_count": "same_direction_fdr_sig_count"}),
            ],
            ignore_index=True,
        )
    candidate_df = candidate_df.dropna(subset=["discovery_logFC", "discovery_adj_p"]).drop_duplicates("gene")
    candidate_df["discovery_direction"] = np.where(candidate_df["discovery_logFC"] >= 0, "up", "down")
    candidate_df["abs_logfc"] = candidate_df["discovery_logFC"].abs()
    candidate_df["neglog10_fdr"] = -np.log10(candidate_df["discovery_adj_p"].clip(lower=1e-300))
    candidate_df["abs_bin"] = pd.qcut(candidate_df["abs_logfc"].rank(method="first"), q=4, labels=False, duplicates="drop")
    candidate_df["fdr_bin"] = pd.qcut(candidate_df["neglog10_fdr"].rank(method="first"), q=4, labels=False, duplicates="drop")
    candidate_df["support"] = candidate_df["same_direction_fdr_sig_count"].fillna(0).astype(int)
    observed = observed.merge(
        candidate_df[["gene", "discovery_direction", "abs_bin", "fdr_bin", "support"]],
        on="gene",
        how="left",
    )

    validation_sets = ["GSE110147", "GSE150910", "GSE53845", "GSE92592"]
    obs_metrics, _ = evaluate_panel(panel, validation_sets)
    obs_summary = {
        "panel_type": "observed_final_panel_regularized_refit",
        "iteration": -1,
        "mean_roc_auc": float(obs_metrics["roc_auc"].mean()),
        "mean_pr_auc": float(obs_metrics["pr_auc"].mean()),
        "mean_balanced_accuracy": float(obs_metrics["balanced_accuracy"].mean()),
        "genes": ";".join(panel),
    }

    rows = [obs_summary]
    genes = candidate_df["gene"].tolist()
    for iteration in range(n_iter):
        selected: list[str] = []
        used: set[str] = set()
        for _, row in observed.sample(frac=1.0, random_state=20260527 + iteration).iterrows():
            masks = [
                (candidate_df["discovery_direction"].eq(row["discovery_direction"]) & candidate_df["support"].eq(row["support"]) & candidate_df["abs_bin"].eq(row["abs_bin"]) & candidate_df["fdr_bin"].eq(row["fdr_bin"])),
                (candidate_df["discovery_direction"].eq(row["discovery_direction"]) & candidate_df["support"].eq(row["support"]) & candidate_df["abs_bin"].eq(row["abs_bin"])),
                (candidate_df["discovery_direction"].eq(row["discovery_direction"]) & candidate_df["support"].eq(row["support"])),
                (candidate_df["discovery_direction"].eq(row["discovery_direction"])),
            ]
            candidate_genes: list[str] = []
            for mask in masks:
                pool = [g for g in candidate_df.loc[mask, "gene"].tolist() if g not in used and g not in panel]
                if pool:
                    candidate_genes = pool
                    break
            if not candidate_genes:
                candidate_genes = [g for g in genes if g not in used and g not in panel]
            chosen = str(RNG.choice(candidate_genes))
            selected.append(chosen)
            used.add(chosen)
        metrics, _ = evaluate_panel(selected, validation_sets)
        rows.append(
            {
                "panel_type": "matched_random_discovery_feature_panel",
                "iteration": iteration,
                "mean_roc_auc": float(metrics["roc_auc"].mean()),
                "mean_pr_auc": float(metrics["pr_auc"].mean()),
                "mean_balanced_accuracy": float(metrics["balanced_accuracy"].mean()),
                "genes": ";".join(selected),
            }
        )
    out = pd.DataFrame(rows)
    rand = out[out["panel_type"].eq("matched_random_discovery_feature_panel")]
    summary = pd.DataFrame(
        [
            {
                "metric": metric,
                "observed_value": obs_summary[metric],
                "matched_random_mean": float(rand[metric].mean()),
                "matched_random_sd": float(rand[metric].std(ddof=1)),
                "matched_random_95th_percentile": float(rand[metric].quantile(0.95)),
                "observed_percentile_within_random": float((rand[metric] <= obs_summary[metric]).mean()),
                "interpretation": "observed final-panel refit is compared with common discovery features matched for direction, validation support, and discovery-effect strata",
            }
            for metric in ["mean_roc_auc", "mean_pr_auc", "mean_balanced_accuracy"]
        ]
    )
    return out, summary


def nonperfect_and_cohort_adjusted() -> tuple[pd.DataFrame, pd.DataFrame]:
    perf = pd.read_csv(ML_DIR / "ml_model_performance_external_validation.csv")
    enet = perf[perf["model"].eq("elastic_net")].copy()
    nonperfect = enet[enet["roc_auc"] < 0.999999].copy()
    nonperfect_summary = pd.DataFrame(
        [
            {
                "stress_test": "exclude_validation_cohorts_with_perfect_roc_auc",
                "retained_cohorts": ";".join(nonperfect["series_id"].tolist()),
                "excluded_cohorts": ";".join(enet.loc[enet["roc_auc"] >= 0.999999, "series_id"].tolist()),
                "mean_roc_auc": float(nonperfect["roc_auc"].mean()),
                "mean_pr_auc": float(nonperfect["pr_auc"].mean()),
                "mean_balanced_accuracy": float(nonperfect["balanced_accuracy"].mean()),
                "interpretation": "summary after removing validation cohorts with perfect ROC AUC; used to contextualize high external discrimination",
            }
        ]
    )

    pred = pd.read_csv(ML_DIR / "ml_external_validation_predictions.csv")
    pred = pred[pred["model"].eq("elastic_net")].copy()
    pred["score_clip"] = pred["score"].clip(1e-5, 1 - 1e-5)
    pred["logit_score"] = np.log(pred["score_clip"] / (1 - pred["score_clip"]))
    cohorts = sorted(pred["series_id"].unique().tolist())
    X_cols = [np.ones(len(pred)), pred["logit_score"].to_numpy(dtype=float)]
    coef_names = ["intercept", "elastic_net_logit_score"]
    for c in cohorts[1:]:
        X_cols.append(pred["series_id"].eq(c).astype(float).to_numpy())
        coef_names.append(f"cohort_{c}")
    X = np.vstack(X_cols).T
    y = pred["true_label"].to_numpy(dtype=float)
    beta = np.zeros(X.shape[1], dtype=float)
    ridge = 1e-4
    for _ in range(100):
        p = sigmoid(X @ beta)
        W = p * (1 - p)
        grad = X.T @ (p - y) + ridge * np.r_[0.0, beta[1:]]
        H = X.T @ (X * W[:, None]) + np.diag(np.r_[0.0, np.repeat(ridge, X.shape[1] - 1)])
        try:
            step = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(H) @ grad
        beta -= step
        if np.max(np.abs(step)) < 1e-7:
            break
    p = sigmoid(X @ beta)
    W = p * (1 - p)
    H = X.T @ (X * W[:, None]) + np.diag(np.r_[0.0, np.repeat(ridge, X.shape[1] - 1)])
    cov = np.linalg.pinv(H)
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    z = beta / np.where(se == 0, np.nan, se)
    pvals = [float(math.erfc(abs(v) / math.sqrt(2))) if np.isfinite(v) else float("nan") for v in z]
    assoc = pd.DataFrame(
        {
            "term": coef_names,
            "coefficient": beta,
            "standard_error": se,
            "wald_z": z,
            "wald_p_value": pvals,
            "model": "disease_status ~ elastic_net_logit_score + validation_cohort_fixed_effects",
            "interpretation": "association test for score after validation-cohort adjustment; not a clinical calibration model",
        }
    )
    return nonperfect_summary, assoc


def draw_auc_stress_figure(disease_summary: pd.DataFrame, disease_scores: pd.DataFrame, matched: pd.DataFrame, matched_summary: pd.DataFrame) -> None:
    width, height = 2500, 1600
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title = load_font(52, True)
    font = load_font(34)
    small = load_font(28)
    tiny = load_font(24)
    red = (181, 93, 96)
    grey = (115, 123, 135)
    blue = (73, 118, 160)
    draw.text((80, 45), "Figure 3. Disease-boundary stress tests for the fibrotic ILD state score", fill=(0, 0, 0), font=title)
    draw.text((80, 105), "A. GSE110147 excluded ILD samples scored with final-panel refit. B. Matched random discovery-feature baseline.", fill=(70, 70, 70), font=small)

    # Panel A
    ax1 = (110, 230, 1120, 1320)
    draw.text((ax1[0], 170), "A  Fibrotic-disease-state score distribution", fill=(0, 0, 0), font=font)
    groups = ["Normal control", "NSIP", "Mixed IPF-NSIP", "IPF"]
    max_score = 1.0
    for tick in np.linspace(0, 1, 6):
        y = int(ax1[3] - tick * (ax1[3] - ax1[1]))
        draw.line((ax1[0], y, ax1[2], y), fill=(232, 237, 244), width=2)
        draw.text((ax1[0] - 55, y - 12), f"{tick:.1f}", fill=(80, 80, 80), font=tiny)
    for i, group in enumerate(groups):
        x = int(ax1[0] + (i + 0.5) * (ax1[2] - ax1[0]) / len(groups))
        vals = disease_scores.loc[disease_scores["stress_group"].eq(group), "refit_final_panel_score"].to_numpy(dtype=float)
        if len(vals) == 0:
            continue
        q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        yq1 = int(ax1[3] - q1 * (ax1[3] - ax1[1]))
        ymed = int(ax1[3] - med * (ax1[3] - ax1[1]))
        yq3 = int(ax1[3] - q3 * (ax1[3] - ax1[1]))
        color = red if group == "IPF" else (170, 115, 80) if "NSIP" in group else grey
        draw.rectangle((x - 45, yq3, x + 45, yq1), outline=color, width=4)
        draw.line((x - 55, ymed, x + 55, ymed), fill=color, width=5)
        for j, v in enumerate(vals):
            jitter = ((j % 9) - 4) * 7
            y = int(ax1[3] - v * (ax1[3] - ax1[1]))
            draw.ellipse((x + jitter - 5, y - 5, x + jitter + 5, y + 5), fill=color)
        label = group.replace(" ", "\n")
        draw.multiline_text((x - 85, ax1[3] + 25), label, fill=(30, 30, 30), font=tiny, align="center")
        draw.text((x - 35, ax1[1] - 35), f"n={len(vals)}", fill=(70, 70, 70), font=tiny)
    draw.text((ax1[0] - 80, ax1[1] - 45), "Score", fill=(0, 0, 0), font=tiny)

    # Panel B
    ax2 = (1350, 230, 2380, 1320)
    draw.text((ax2[0], 170), "B  Matched random-panel baseline", fill=(0, 0, 0), font=font)
    rand = matched[matched["panel_type"].eq("matched_random_discovery_feature_panel")]["mean_roc_auc"].to_numpy(dtype=float)
    observed = float(matched.loc[matched["panel_type"].eq("observed_final_panel_regularized_refit"), "mean_roc_auc"].iloc[0])
    bins = np.linspace(max(0.4, rand.min() - 0.02), min(1.0, max(rand.max(), observed) + 0.02), 24)
    counts, edges = np.histogram(rand, bins=bins)
    max_count = max(counts.max(), 1)
    for i, count in enumerate(counts):
        x0 = int(ax2[0] + (edges[i] - bins[0]) / (bins[-1] - bins[0]) * (ax2[2] - ax2[0]))
        x1 = int(ax2[0] + (edges[i + 1] - bins[0]) / (bins[-1] - bins[0]) * (ax2[2] - ax2[0]))
        y0 = int(ax2[3] - count / max_count * (ax2[3] - ax2[1]))
        draw.rectangle((x0, y0, x1 - 2, ax2[3]), fill=(193, 209, 225), outline=(255, 255, 255))
    obs_x = int(ax2[0] + (observed - bins[0]) / (bins[-1] - bins[0]) * (ax2[2] - ax2[0]))
    draw.line((obs_x, ax2[1], obs_x, ax2[3]), fill=red, width=6)
    q95 = float(rand.quantile(0.95) if isinstance(rand, pd.Series) else np.quantile(rand, 0.95))
    q95_x = int(ax2[0] + (q95 - bins[0]) / (bins[-1] - bins[0]) * (ax2[2] - ax2[0]))
    draw.line((q95_x, ax2[1], q95_x, ax2[3]), fill=blue, width=4)
    for tick in np.linspace(round(bins[0], 2), round(bins[-1], 2), 5):
        x = int(ax2[0] + (tick - bins[0]) / (bins[-1] - bins[0]) * (ax2[2] - ax2[0]))
        draw.text((x - 28, ax2[3] + 25), f"{tick:.2f}", fill=(50, 50, 50), font=tiny)
    draw.text((ax2[0] + 260, ax2[3] + 75), "Mean external ROC AUC", fill=(0, 0, 0), font=small)
    draw.text((ax2[0] + 25, ax2[1] + 25), f"Observed refit={observed:.3f}", fill=red, font=small)
    draw.text((ax2[0] + 25, ax2[1] + 70), f"Matched random 95th={q95:.3f}", fill=blue, font=small)
    image.save(FIG_DIR / "Additional_Figure_S13_auc_stress_tests.png")
    image.save(AUC_OUT / "additional_figure_s13_auc_stress_tests.pdf", "PDF", resolution=300.0)


def hypergeom_sf(k: int, K: int, n: int, N: int) -> float:
    if any(v < 0 for v in [k, K, n, N]) or K > N or n > N:
        return float("nan")
    max_i = min(K, n)
    if k > max_i:
        return 0.0
    denom = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
    logs = []
    for i in range(k, max_i + 1):
        if K - i > N - n:
            continue
        logp = (
            math.lgamma(K + 1)
            - math.lgamma(i + 1)
            - math.lgamma(K - i + 1)
            + math.lgamma(N - K + 1)
            - math.lgamma(n - i + 1)
            - math.lgamma(N - K - n + i + 1)
            - denom
        )
        logs.append(logp)
    if not logs:
        return 0.0
    m = max(logs)
    return float(min(1.0, math.exp(m) * sum(math.exp(x - m) for x in logs)))


def odds_ratio(k: int, K: int, n: int, N: int) -> float:
    a = k + 0.5
    b = n - k + 0.5
    c = K - k + 0.5
    d = N - K - n + k + 0.5
    return float((a * d) / (b * c))


def build_mirtarbase_targets(robust_mirnas: pd.DataFrame) -> dict[tuple[str, str], set[str]]:
    mti = pd.read_csv(MIRTARBASE_PATH)
    mti = mti[mti["Species (miRNA)"].eq("hsa") & mti["Species (Target Gene)"].eq("hsa")].copy()
    mti = mti[mti["Support Type"].astype(str).str.startswith("Functional MTI", na=False)]
    mti["target_gene_norm"] = mti["Target Gene"].map(normalize_gene)
    mti = mti[mti["target_gene_norm"] != ""]
    mti["mirna_norm"] = mti["miRNA"].map(normalize_mirna)
    mti["mirna_base"] = mti["miRNA"].map(mirna_base)
    targets: dict[tuple[str, str], set[str]] = {}
    for mirna in robust_mirnas["standard_feature_id"].map(normalize_mirna):
        exact = set(mti.loc[mti["mirna_norm"].eq(mirna), "target_gene_norm"])
        base = set(mti.loc[mti["mirna_base"].eq(mirna_base(mirna)), "target_gene_norm"])
        targets[(mirna, "miRTarBase_exact")] = exact
        targets[(mirna, "miRTarBase_arm_recoverable")] = base
    return targets


def mirna_program_support() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    robust_mirna = pd.read_csv(ROBUST_DIR / "robust_mirna_candidates_strict.csv")
    robust_mrna = pd.read_csv(ROBUST_DIR / "robust_mrna_candidates_strict.csv")
    robust_mirna["standard_feature_id"] = robust_mirna["standard_feature_id"].map(normalize_mirna)
    robust_mrna["gene"] = robust_mrna["standard_feature_id"].map(normalize_gene)
    down_mirnas = robust_mirna[robust_mirna["discovery_direction"].eq("down")].copy()
    up_genes = set(robust_mrna.loc[robust_mrna["discovery_direction"].eq("up"), "gene"].dropna())
    gene_universe = set(read_feature_matrix("GSE32537", common=False).columns) - {"sample_id", "group", "label"}
    up_genes = up_genes & gene_universe
    targets = build_mirtarbase_targets(down_mirnas)
    rows = []
    for (mirna, source), target_set in targets.items():
        t = set(target_set) & gene_universe
        k = len(t & up_genes)
        K = len(t)
        n = len(up_genes)
        N = len(gene_universe)
        p = hypergeom_sf(k, K, n, N) if K > 0 else float("nan")
        rows.append(
            {
                "candidate_mirna": mirna,
                "mirna_direction": "down",
                "target_source": source,
                "target_genes_in_background": K,
                "upregulated_robust_mrna_targets": k,
                "upregulated_robust_mrna_foreground": n,
                "gene_background": N,
                "odds_ratio": odds_ratio(k, K, n, N) if K > 0 else float("nan"),
                "fisher_hypergeom_p": p,
                "target_genes_overlapping_upregulated_robust_mRNAs": ";".join(sorted(t & up_genes)),
            }
        )
    enrich = pd.DataFrame(rows)
    enrich["fdr_bh"] = bh_adjust(enrich["fisher_hypergeom_p"].tolist())
    enrich["interpretation"] = np.where(
        enrich["fdr_bh"].le(0.1),
        "target-set-level support for derepression-like signal among IPF-upregulated robust mRNAs",
        "no FDR-supported target-set enrichment in this background; retain as hypothesis-generating",
    )

    de = pd.read_csv(DE_DIR / "GSE32537_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv")
    de["gene"] = de["standard_feature_id"].map(normalize_gene)
    de = de[de["gene"].isin(gene_universe)].drop_duplicates("gene")
    h375_targets = targets.get(("hsa-miR-375", "miRTarBase_arm_recoverable"), set()) & set(de["gene"])
    de["is_hsa_miR_375_target"] = de["gene"].isin(h375_targets)
    target_vals = de.loc[de["is_hsa_miR_375_target"], "logFC"].dropna().to_numpy(dtype=float)
    bg_vals = de.loc[~de["is_hsa_miR_375_target"], "logFC"].dropna().to_numpy(dtype=float)
    observed_diff = float(np.mean(target_vals) - np.mean(bg_vals)) if len(target_vals) and len(bg_vals) else float("nan")
    combined = np.r_[target_vals, bg_vals]
    n_t = len(target_vals)
    perm_diffs = []
    for _ in range(5000):
        shuffled = RNG.permutation(combined)
        perm_diffs.append(float(np.mean(shuffled[:n_t]) - np.mean(shuffled[n_t:])))
    perm_diffs = np.asarray(perm_diffs)
    p_emp = float((np.sum(perm_diffs >= observed_diff) + 1) / (len(perm_diffs) + 1)) if np.isfinite(observed_diff) else float("nan")
    repression = pd.DataFrame(
        [
            {
                "candidate_mirna": "hsa-miR-375",
                "target_source": "miRTarBase_arm_recoverable",
                "targets_in_discovery_mrna_background": n_t,
                "mean_target_discovery_logFC": float(np.mean(target_vals)) if len(target_vals) else float("nan"),
                "mean_non_target_discovery_logFC": float(np.mean(bg_vals)) if len(bg_vals) else float("nan"),
                "target_minus_non_target_mean_logFC": observed_diff,
                "permutation_iterations": 5000,
                "one_sided_permutation_p_target_greater_than_background": p_emp,
                "targets_upregulated_fraction": float((target_vals > 0).mean()) if len(target_vals) else float("nan"),
                "interpretation": "program-level derepression support if targets show higher discovery logFC than non-target genes",
            }
        ]
    )

    axes = pd.read_csv(PROJECT_DIR / "results" / "submission_enhancements" / "mirna_mrna_axes_evidence_graded.csv")
    relaxed = axes.copy()
    relaxed["sensitivity_layer"] = "relaxed_validated_axis_sensitivity"
    relaxed["sensitivity_interpretation"] = np.where(
        relaxed["match_type"].eq("exact"),
        "exact mature/source-name axes retained in main text",
        "arm-recoverable miRTarBase axes retained as exploratory sensitivity only",
    )

    audit = pd.DataFrame(
        [
            {
                "analysis": "paired GSE32537 mRNA and GSE32538 miRNA inverse correlation",
                "status": "not_attempted",
                "reason": "No explicit shared donor identifier was available in curated metadata; GEO sample accessions differ between mRNA and miRNA series, so ordinal matching was not used.",
                "manuscript_interpretation": "miRNA-mRNA relationships are interpreted at cohort/target-set level rather than paired donor correlation level.",
            }
        ]
    )
    qc = pd.DataFrame(
        [
            {
                "qc_item": "robust_downregulated_mirnas_tested",
                "value": down_mirnas["standard_feature_id"].nunique(),
                "pass": down_mirnas["standard_feature_id"].nunique() > 0,
            },
            {
                "qc_item": "mirtarbase_functional_rows_used",
                "value": sum(len(v) for v in targets.values()),
                "pass": sum(len(v) for v in targets.values()) > 0,
            },
            {
                "qc_item": "hsa_mir_375_targets_in_background",
                "value": int(n_t),
                "pass": int(n_t) > 0,
            },
            {
                "qc_item": "paired_correlation_attempted_only_if_metadata_supported",
                "value": "not_attempted_no_explicit_shared_donor_id",
                "pass": True,
            },
        ]
    )
    return enrich, repression, relaxed, audit, qc


def draw_mirna_program_figure(enrich: pd.DataFrame, repression: pd.DataFrame) -> None:
    width, height = 2500, 1450
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title = load_font(52, True)
    font = load_font(34)
    small = load_font(28)
    tiny = load_font(24)
    red = (181, 93, 96)
    blue = (73, 118, 160)
    grey = (115, 123, 135)
    draw.text((80, 45), "Figure 4. miRNA evidence-gate stress tests", fill=(0, 0, 0), font=title)
    draw.text((80, 105), "A. Target-set enrichment among IPF-upregulated robust mRNAs. B. hsa-miR-375 target release-like score.", fill=(70, 70, 70), font=small)

    top = enrich[enrich["target_source"].eq("miRTarBase_arm_recoverable")].copy()
    top = top.sort_values(["fdr_bh", "odds_ratio"], ascending=[True, False]).head(8)
    ax1 = (130, 230, 1180, 1180)
    draw.text((ax1[0], 170), "A  Robust downregulated miRNAs", fill=(0, 0, 0), font=font)
    values = top["odds_ratio"].replace([np.inf, -np.inf], np.nan).fillna(0).to_numpy(dtype=float)
    max_v = max(1.6, float(np.nanmax(values)) if len(values) else 1.6)
    ref_x = int(ax1[0] + 300 + (1.0 / max_v) * (ax1[2] - ax1[0] - 420))
    draw.line((ref_x, ax1[1] - 15, ref_x, ax1[3] - 180), fill=(160, 160, 160), width=3)
    draw.text((ref_x - 32, ax1[3] - 165), "OR=1", fill=(90, 90, 90), font=tiny)
    for i, (_, row) in enumerate(top.iterrows()):
        y = int(ax1[1] + i * 105)
        val = float(row["odds_ratio"]) if pd.notna(row["odds_ratio"]) else 0.0
        bar_w = int((val / max_v) * (ax1[2] - ax1[0] - 420))
        draw.text((ax1[0], y + 12), str(row["candidate_mirna"]), fill=(0, 0, 0), font=small)
        color = red if row["candidate_mirna"] == "hsa-miR-375" else blue
        if bar_w > 0:
            draw.rectangle((ax1[0] + 300, y + 8, ax1[0] + 300 + bar_w, y + 56), fill=color)
        fdr_text = "NA" if pd.isna(row["fdr_bh"]) else f"{row['fdr_bh']:.2g}"
        or_text = "NA" if pd.isna(row["odds_ratio"]) else f"{row['odds_ratio']:.2f}"
        draw.text((ax1[0] + 315 + max(bar_w, 4), y + 12), f"OR={or_text}; FDR={fdr_text}", fill=(40, 40, 40), font=tiny)
    draw.text((ax1[0] + 300, ax1[3] + 15), "Odds ratio for target overlap", fill=(0, 0, 0), font=small)

    ax2 = (1420, 230, 2320, 1180)
    draw.text((ax2[0], 170), "B  hsa-miR-375 target release-like score", fill=(0, 0, 0), font=font)
    row = repression.iloc[0]
    labels = ["Targets", "Non-targets"]
    vals = [float(row["mean_target_discovery_logFC"]), float(row["mean_non_target_discovery_logFC"])]
    min_v = min(vals + [0]) - 0.2
    max_v = max(vals + [0]) + 0.2
    zero_y = int(ax2[3] - (0 - min_v) / (max_v - min_v) * (ax2[3] - ax2[1]))
    draw.line((ax2[0], zero_y, ax2[2], zero_y), fill=(170, 170, 170), width=3)
    for i, (lab, val) in enumerate(zip(labels, vals)):
        x = int(ax2[0] + (i + 0.5) * (ax2[2] - ax2[0]) / 2)
        y = int(ax2[3] - (val - min_v) / (max_v - min_v) * (ax2[3] - ax2[1]))
        color = red if lab == "Targets" else grey
        draw.rectangle((x - 80, min(y, zero_y), x + 80, max(y, zero_y)), fill=color)
        draw.text((x - 70, ax2[3] + 25), lab, fill=(0, 0, 0), font=small)
        draw.text((x - 65, y - 40), f"{val:.2f}", fill=(0, 0, 0), font=small)
    draw.text((ax2[0] + 25, ax2[1] + 35), f"Permutation p={row['one_sided_permutation_p_target_greater_than_background']:.3g}", fill=(40, 40, 40), font=small)
    draw.text((ax2[0] + 25, ax2[1] + 80), f"n targets={int(row['targets_in_discovery_mrna_background'])}", fill=(40, 40, 40), font=small)
    draw.text((ax2[0] + 60, ax2[3] + 85), "Mean discovery mRNA logFC", fill=(0, 0, 0), font=small)

    image.save(FIG_DIR / "Additional_Figure_S14_mirna_program_stress_test.png")
    image.save(MIRNA_OUT / "additional_figure_s14_mirna_program_stress_test.pdf", "PDF", resolution=300.0)


def main() -> None:
    panel = pd.read_csv(ML_DIR / "ml_final_biomarker_panel.csv")["feature"].map(normalize_gene).tolist()
    disease_scores, disease_summary = score_gse110147_disease_controls(panel)
    disease_scores.to_csv(AUC_OUT / "gse110147_excluded_ild_disease_control_scores.csv", index=False, encoding="utf-8-sig")
    disease_summary.to_csv(AUC_OUT / "gse110147_excluded_ild_disease_control_summary.csv", index=False, encoding="utf-8-sig")

    matched, matched_summary = matched_random_panel_baseline(panel, n_iter=500)
    matched.to_csv(AUC_OUT / "matched_random_discovery_feature_panel_baseline.csv", index=False, encoding="utf-8-sig")
    matched_summary.to_csv(AUC_OUT / "matched_random_discovery_feature_panel_summary.csv", index=False, encoding="utf-8-sig")

    nonperfect, cohort_adj = nonperfect_and_cohort_adjusted()
    nonperfect.to_csv(AUC_OUT / "nonperfect_cohort_stress_summary.csv", index=False, encoding="utf-8-sig")
    cohort_adj.to_csv(AUC_OUT / "cohort_adjusted_external_logistic_association.csv", index=False, encoding="utf-8-sig")
    auc_qc = pd.DataFrame(
        [
            {"qc_item": "disease_control_scores_include_excluded_nsip_mixed_samples", "value": disease_scores["stress_group"].isin(["NSIP", "Mixed IPF-NSIP"]).sum(), "pass": disease_scores["stress_group"].isin(["NSIP", "Mixed IPF-NSIP"]).sum() == 15},
            {"qc_item": "matched_random_iterations", "value": int((matched["panel_type"] == "matched_random_discovery_feature_panel").sum()), "pass": int((matched["panel_type"] == "matched_random_discovery_feature_panel").sum()) == 500},
            {"qc_item": "nonperfect_cohort_summary_created", "value": str(nonperfect["retained_cohorts"].iloc[0]), "pass": True},
            {"qc_item": "cohort_adjusted_score_term_present", "value": "elastic_net_logit_score" in cohort_adj["term"].tolist(), "pass": "elastic_net_logit_score" in cohort_adj["term"].tolist()},
        ]
    )
    auc_qc.to_csv(AUC_OUT / "auc_stress_test_qc.csv", index=False, encoding="utf-8-sig")
    draw_auc_stress_figure(disease_summary, disease_scores, matched, matched_summary)
    pd.DataFrame(
        [
            {
                "figure": "Figure 3",
                "png": "manuscript/figures/Additional_Figure_S13_auc_stress_tests.png",
                "pdf": "results/ml_stress_tests/additional_figure_s13_auc_stress_tests.pdf",
                "description": "Disease-control score distribution for excluded GSE110147 NSIP/mixed samples and matched random discovery-feature panel baseline.",
            }
        ]
    ).to_csv(AUC_OUT / "figure_s13_manifest.csv", index=False, encoding="utf-8-sig")

    enrich, repression, relaxed, audit, mirna_qc = mirna_program_support()
    enrich.to_csv(MIRNA_OUT / "robust_mirna_target_set_enrichment.csv", index=False, encoding="utf-8-sig")
    repression.to_csv(MIRNA_OUT / "hsa_mir_375_target_repression_release_score.csv", index=False, encoding="utf-8-sig")
    relaxed.to_csv(MIRNA_OUT / "relaxed_mirna_axis_sensitivity_set.csv", index=False, encoding="utf-8-sig")
    audit.to_csv(MIRNA_OUT / "paired_mirna_mrna_correlation_audit.csv", index=False, encoding="utf-8-sig")
    mirna_qc.to_csv(MIRNA_OUT / "mirna_program_support_qc.csv", index=False, encoding="utf-8-sig")
    draw_mirna_program_figure(enrich, repression)
    pd.DataFrame(
        [
            {
                "figure": "Figure 4",
                "png": "manuscript/figures/Additional_Figure_S14_mirna_program_stress_test.png",
                "pdf": "results/mirna_program_support/additional_figure_s14_mirna_program_stress_test.pdf",
                "description": "Validated miRTarBase target-set stress test for robust downregulated miRNAs and hsa-miR-375 target release-like score.",
            }
        ]
    ).to_csv(MIRNA_OUT / "figure_s14_manifest.csv", index=False, encoding="utf-8-sig")

    print(AUC_OUT)
    print(MIRNA_OUT)


if __name__ == "__main__":
    main()
