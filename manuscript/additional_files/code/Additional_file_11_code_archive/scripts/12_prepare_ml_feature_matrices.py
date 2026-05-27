#!/usr/bin/env python
"""
Prepare gene-level machine-learning feature matrices for IPF vs Control.

Inputs:
  - data_processed/expression/*_expression_in_annotation_order.csv.gz
  - metadata/all_bulk_mirna_annotation.csv
  - results/feature_annotation/feature_annotation_map.csv
  - results/robust_candidates/robust_mrna_candidates_strict.csv

Outputs:
  - results/models/feature_matrices/{GSE}_gene_level_matrix.csv
  - results/models/feature_matrices/{GSE}_robust_mrna_matrix.csv
  - results/models/feature_matrices/{GSE}_common_robust_mrna_matrix.csv
  - results/models/ml_feature_matrix_qc.csv
  - results/models/common_robust_mrna_features.txt
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
EXPRESSION_DIR = PROJECT_DIR / "data_processed" / "expression"
METADATA_PATH = PROJECT_DIR / "metadata" / "all_bulk_mirna_annotation.csv"
FEATURE_MAP_PATH = PROJECT_DIR / "results" / "feature_annotation" / "feature_annotation_map.csv"
ROBUST_MRNA_PATH = PROJECT_DIR / "results" / "robust_candidates" / "robust_mrna_candidates_strict.csv"
DISCOVERY_SIG_PATH = PROJECT_DIR / "results" / "differential_expression_annotated" / "GSE32537_IPF_vs_Control_de_significant_annotated_gene_or_mirna_level_fdr0.05_logfc1.csv"
OUT_DIR = PROJECT_DIR / "results" / "models" / "feature_matrices"
MODELS_DIR = PROJECT_DIR / "results" / "models"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRAINING_DATASET = "GSE32537"
VALIDATION_DATASETS = ["GSE110147", "GSE150910", "GSE53845", "GSE92592"]
DATASETS = [TRAINING_DATASET] + VALIDATION_DATASETS


def normalize_gene(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.upper()
    return text if text and text != "NA" else ""


def read_expression(series_id: str) -> pd.DataFrame:
    path = EXPRESSION_DIR / f"{series_id}_expression_in_annotation_order.csv.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    expr = pd.read_csv(path)
    feature_col = expr.columns[0]
    expr = expr.rename(columns={feature_col: "feature_id"})
    expr["feature_id"] = expr["feature_id"].astype(str)
    return expr


def collapse_to_gene_level(expr: pd.DataFrame, fmap: pd.DataFrame, series_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = fmap[(fmap["series_id"] == series_id) & (fmap["is_annotated"].astype(str).str.lower().isin(["true", "1", "yes"]))].copy()
    mapping["feature_id"] = mapping["feature_id"].astype(str)
    mapping["standard_feature_id"] = mapping["standard_feature_id"].map(normalize_gene)
    mapping = mapping[mapping["standard_feature_id"] != ""]
    mapping = mapping.drop_duplicates(subset=["feature_id"], keep="first")

    merged = expr.merge(mapping[["feature_id", "standard_feature_id"]], on="feature_id", how="inner")
    sample_cols = [c for c in merged.columns if c not in {"feature_id", "standard_feature_id"}]
    values = merged[sample_cols].apply(pd.to_numeric, errors="coerce")
    merged[sample_cols] = values
    merged["row_variance"] = merged[sample_cols].var(axis=1, skipna=True)
    merged = merged.sort_values(["standard_feature_id", "row_variance"], ascending=[True, False])
    collapsed = merged.drop_duplicates(subset=["standard_feature_id"], keep="first").copy()

    collapse_qc = pd.DataFrame(
        [
            {
                "series_id": series_id,
                "expression_features": expr.shape[0],
                "mapped_feature_rows": merged.shape[0],
                "unique_gene_symbols": collapsed["standard_feature_id"].nunique(),
                "multi_probe_genes": int(merged["standard_feature_id"].duplicated().sum()),
                "collapse_rule": "highest within-dataset variance per gene symbol",
            }
        ]
    )

    gene_expr = collapsed[["standard_feature_id"] + sample_cols].copy()
    gene_expr = gene_expr.set_index("standard_feature_id").T
    gene_expr.index.name = "sample_id"
    gene_expr = gene_expr.reset_index()
    return gene_expr, collapse_qc


def add_labels(gene_expr: pd.DataFrame, annotation: pd.DataFrame, series_id: str) -> pd.DataFrame:
    pheno = annotation[
        (annotation["series_id"] == series_id)
        & (annotation["include"].astype(str).str.lower() == "yes")
        & (annotation["group"].isin(["IPF", "Control"]))
    ][["sample_id", "group"]].copy()
    pheno = pheno.drop_duplicates(subset=["sample_id"], keep="first")
    out = pheno.merge(gene_expr, on="sample_id", how="inner")
    out["label"] = (out["group"] == "IPF").astype(int)
    cols = ["sample_id", "group", "label"] + [c for c in out.columns if c not in {"sample_id", "group", "label"}]
    return out[cols]


def write_feature_list(path: Path, genes: list[str]) -> None:
    path.write_text("\n".join(genes) + "\n", encoding="utf-8")


def main() -> None:
    annotation = pd.read_csv(METADATA_PATH)
    fmap = pd.read_csv(FEATURE_MAP_PATH)
    robust = pd.read_csv(ROBUST_MRNA_PATH)
    robust_genes = sorted({normalize_gene(g) for g in robust["standard_feature_id"] if normalize_gene(g)})
    discovery_sig = pd.read_csv(DISCOVERY_SIG_PATH)
    discovery_genes = sorted({normalize_gene(g) for g in discovery_sig["standard_feature_id"] if normalize_gene(g)})

    qc_rows = []
    collapse_qcs = []
    available_robust_by_dataset = {}
    available_discovery_by_dataset = {}

    for series_id in DATASETS:
        expr = read_expression(series_id)
        gene_expr, collapse_qc = collapse_to_gene_level(expr, fmap, series_id)
        labeled = add_labels(gene_expr, annotation, series_id)

        gene_cols = [c for c in labeled.columns if c not in {"sample_id", "group", "label"}]
        available_robust = sorted(set(gene_cols).intersection(robust_genes))
        available_discovery = sorted(set(gene_cols).intersection(discovery_genes))
        available_robust_by_dataset[series_id] = set(available_robust)
        available_discovery_by_dataset[series_id] = set(available_discovery)

        robust_matrix = labeled[["sample_id", "group", "label"] + available_robust].copy()
        discovery_matrix = labeled[["sample_id", "group", "label"] + available_discovery].copy()

        gene_level_path = OUT_DIR / f"{series_id}_gene_level_matrix.csv"
        robust_path = OUT_DIR / f"{series_id}_robust_mrna_matrix.csv"
        discovery_path = OUT_DIR / f"{series_id}_discovery_only_mrna_matrix.csv"
        labeled.to_csv(gene_level_path, index=False, encoding="utf-8-sig")
        robust_matrix.to_csv(robust_path, index=False, encoding="utf-8-sig")
        discovery_matrix.to_csv(discovery_path, index=False, encoding="utf-8-sig")

        numeric = robust_matrix[available_robust].apply(pd.to_numeric, errors="coerce")
        qc_rows.append(
            {
                "series_id": series_id,
                "role": "training" if series_id == TRAINING_DATASET else "external_validation",
                "samples": robust_matrix.shape[0],
                "ipf_samples": int((robust_matrix["group"] == "IPF").sum()),
                "control_samples": int((robust_matrix["group"] == "Control").sum()),
                "gene_level_features": len(gene_cols),
                "robust_mrna_features_available": len(available_robust),
                "robust_mrna_features_missing": len(set(robust_genes) - set(available_robust)),
                "discovery_only_features_available": len(available_discovery),
                "discovery_only_features_missing": len(set(discovery_genes) - set(available_discovery)),
                "numeric_missing_fraction": float(numeric.isna().to_numpy().mean()) if numeric.size else np.nan,
                "finite_fraction": float(np.isfinite(numeric.to_numpy(dtype=float)).mean()) if numeric.size else np.nan,
                "feature_matrix_path": str(robust_path),
                "qc1_labels_pass": robust_matrix["group"].isin(["IPF", "Control"]).all()
                and robust_matrix["label"].isin([0, 1]).all()
                and robust_matrix["sample_id"].is_unique
                and robust_matrix["group"].nunique() == 2,
                "qc2_features_pass": len(available_robust) >= 20,
                "qc3_numeric_pass": numeric.size > 0 and numeric.isna().to_numpy().mean() == 0,
            }
        )
        collapse_qcs.append(collapse_qc)

    common_features = sorted(set.intersection(*available_robust_by_dataset.values()))
    common_discovery_features = sorted(set.intersection(*available_discovery_by_dataset.values()))
    write_feature_list(MODELS_DIR / "common_robust_mrna_features.txt", common_features)
    write_feature_list(MODELS_DIR / "common_discovery_only_mrna_features.txt", common_discovery_features)

    for series_id in DATASETS:
        robust_path = OUT_DIR / f"{series_id}_robust_mrna_matrix.csv"
        mat = pd.read_csv(robust_path)
        common_mat = mat[["sample_id", "group", "label"] + common_features].copy()
        common_path = OUT_DIR / f"{series_id}_common_robust_mrna_matrix.csv"
        common_mat.to_csv(common_path, index=False, encoding="utf-8-sig")
        discovery_path = OUT_DIR / f"{series_id}_discovery_only_mrna_matrix.csv"
        discovery_mat = pd.read_csv(discovery_path)
        common_discovery_mat = discovery_mat[["sample_id", "group", "label"] + common_discovery_features].copy()
        common_discovery_path = OUT_DIR / f"{series_id}_common_discovery_only_mrna_matrix.csv"
        common_discovery_mat.to_csv(common_discovery_path, index=False, encoding="utf-8-sig")

    qc = pd.DataFrame(qc_rows)
    qc["common_robust_features_across_all_datasets"] = len(common_features)
    qc["common_discovery_only_features_across_all_datasets"] = len(common_discovery_features)
    qc["triple_qc_pass"] = qc["qc1_labels_pass"] & qc["qc2_features_pass"] & qc["qc3_numeric_pass"] & (len(common_features) >= 20) & (len(common_discovery_features) >= 20)
    qc.to_csv(MODELS_DIR / "ml_feature_matrix_qc.csv", index=False, encoding="utf-8-sig")
    pd.concat(collapse_qcs, ignore_index=True).to_csv(MODELS_DIR / "ml_gene_collapse_qc.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(MODELS_DIR / "ml_feature_matrix_qc.xlsx", engine="openpyxl") as writer:
        qc.to_excel(writer, sheet_name="feature_matrix_qc", index=False)
        pd.concat(collapse_qcs, ignore_index=True).to_excel(writer, sheet_name="gene_collapse_qc", index=False)
        pd.DataFrame({"common_feature": common_features}).to_excel(writer, sheet_name="common_features", index=False)
        pd.DataFrame({"common_discovery_only_feature": common_discovery_features}).to_excel(writer, sheet_name="common_discovery_only", index=False)

    print(qc.to_string(index=False))
    print(f"\nCommon robust mRNA features across all datasets: {len(common_features)}")
    print(f"Common discovery-only mRNA features across all datasets: {len(common_discovery_features)}")


if __name__ == "__main__":
    main()
