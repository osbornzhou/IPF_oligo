#!/usr/bin/env python
"""
Integrate robust IPF miRNA and mRNA candidates with miRTarBase targets.

The script keeps candidate axes where:
  - miRNA is a robust IPF candidate
  - mRNA is a robust IPF candidate
  - miRTarBase supports the miRNA-target relationship
  - miRNA and mRNA show opposite discovery directions

Because some array annotations lack mature arm suffixes, matching is performed
with two levels:
  - exact: hsa-miR-30d-5p == hsa-miR-30d-5p
  - arm_agnostic: hsa-miR-30d == hsa-miR-30d-5p or hsa-miR-30d-3p
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
ROBUST_DIR = PROJECT_DIR / "results" / "robust_candidates"
MIRTARBASE_PATH = PROJECT_DIR / "data_external" / "miRTarBase" / "hsa_MTI_miRTarBase_2025_v10.csv"
OUTPUT_DIR = PROJECT_DIR / "results" / "mirna_mrna_axes"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def normalize_mirna(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.replace("hsa-mir-", "hsa-miR-")
    text = re.sub(r"_st$", "", text)
    return text


def mirna_core(value: object) -> str:
    text = normalize_mirna(value)
    text = re.sub(r"-(3p|5p)$", "", text)
    return text


def direction(value: float) -> str:
    if value > 0:
        return "up"
    if value < 0:
        return "down"
    return "zero"


def support_rank(support_type: object, experiments: object) -> int:
    support = "" if pd.isna(support_type) else str(support_type)
    exp = "" if pd.isna(experiments) else str(experiments)
    rank = 0
    if "Functional MTI" in support:
        rank += 10
    if "Non-Functional" in support:
        rank -= 2
    strong_terms = [
        "Luciferase reporter assay",
        "Reporter assay",
        "Western blot",
        "qRT-PCR",
        "RT-PCR",
    ]
    rank += sum(2 for term in strong_terms if term.lower() in exp.lower())
    if "CLIP" in exp.upper() or "PAR-CLIP" in exp.upper():
        rank += 1
    return rank


def read_candidates() -> tuple[pd.DataFrame, pd.DataFrame]:
    mrna = pd.read_csv(ROBUST_DIR / "robust_mrna_candidates_strict.csv")
    mirna = pd.read_csv(ROBUST_DIR / "robust_mirna_candidates_strict.csv")

    mrna["target_gene"] = mrna["standard_feature_id"].astype(str).str.upper()
    mrna["mrna_direction"] = mrna["discovery_logFC"].apply(direction)
    mirna["mirna_name_norm"] = mirna["standard_feature_id"].map(normalize_mirna)
    mirna["mirna_core"] = mirna["mirna_name_norm"].map(mirna_core)
    mirna["mirna_direction"] = mirna["discovery_logFC"].apply(direction)
    return mrna, mirna


def read_mirtarbase() -> pd.DataFrame:
    if not MIRTARBASE_PATH.exists():
        raise FileNotFoundError(
            f"miRTarBase file not found: {MIRTARBASE_PATH}. "
            "Run the download step or place hsa_MTI.csv there."
        )
    cols = [
        "miRTarBase ID",
        "miRNA",
        "Species (miRNA)",
        "Target Gene",
        "Target Gene (Entrez ID)",
        "Species (Target Gene)",
        "Experiments",
        "Support Type",
        "References (PMID)",
    ]
    df = pd.read_csv(MIRTARBASE_PATH, usecols=cols)
    df = df[
        (df["Species (miRNA)"].astype(str) == "hsa")
        & (df["Species (Target Gene)"].astype(str) == "hsa")
    ].copy()
    df["miRNA_norm"] = df["miRNA"].map(normalize_mirna)
    df["miRNA_core"] = df["miRNA_norm"].map(mirna_core)
    df["target_gene"] = df["Target Gene"].astype(str).str.upper()
    df["support_rank"] = [
        support_rank(support, experiments)
        for support, experiments in zip(df["Support Type"], df["Experiments"])
    ]
    return df


def collapse_targets(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["candidate_mirna", "target_gene", "match_type"]
    agg = (
        df.groupby(group_cols, as_index=False)
        .agg(
            mirtarbase_mirnas=("miRNA_norm", lambda x: ";".join(sorted(set(map(str, x))))),
            mirtarbase_ids=("miRTarBase ID", lambda x: ";".join(sorted(set(map(str, x))))),
            support_types=("Support Type", lambda x: ";".join(sorted(set(map(str, x))))),
            experiments=("Experiments", lambda x: ";".join(sorted(set(map(str, x))))),
            pmids=("References (PMID)", lambda x: ";".join(sorted(set(str(v).replace(".0", "") for v in x if pd.notna(v))))),
            evidence_count=("miRTarBase ID", "size"),
            best_support_rank=("support_rank", "max"),
        )
    )
    return agg


def main() -> None:
    mrna, mirna = read_candidates()
    mti = read_mirtarbase()

    exact = mti.merge(
        mirna[["standard_feature_id", "mirna_name_norm", "mirna_core", "discovery_logFC", "discovery_adj_p", "mirna_direction"]],
        left_on="miRNA_norm",
        right_on="mirna_name_norm",
        how="inner",
    )
    exact["candidate_mirna"] = exact["standard_feature_id"]
    exact["match_type"] = "exact"

    arm_agnostic = mti.merge(
        mirna[["standard_feature_id", "mirna_name_norm", "mirna_core", "discovery_logFC", "discovery_adj_p", "mirna_direction"]],
        left_on="miRNA_core",
        right_on="mirna_core",
        how="inner",
    )
    arm_agnostic = arm_agnostic[arm_agnostic["miRNA_norm"] != arm_agnostic["mirna_name_norm"]].copy()
    arm_agnostic["candidate_mirna"] = arm_agnostic["standard_feature_id"]
    arm_agnostic["match_type"] = "arm_agnostic"

    matched_mti = pd.concat([exact, arm_agnostic], ignore_index=True)
    if matched_mti.empty:
        raise RuntimeError("No miRTarBase interactions matched robust miRNA candidates.")

    target_hits = matched_mti.merge(
        mrna[
            [
                "standard_feature_id",
                "target_gene",
                "discovery_logFC",
                "discovery_adj_p",
                "mrna_direction",
                "same_direction_fdr_sig_count",
                "robust_score",
            ]
        ],
        on="target_gene",
        how="inner",
        suffixes=("_mirna", "_mrna"),
    )

    target_hits["opposite_direction"] = target_hits["mirna_direction"] != target_hits["mrna_direction"]
    target_hits = target_hits[target_hits["opposite_direction"]].copy()

    collapsed = collapse_targets(target_hits)
    axes = collapsed.merge(
        mirna[
            [
                "standard_feature_id",
                "discovery_logFC",
                "discovery_adj_p",
                "mirna_direction",
                "same_direction_fdr_sig_count",
                "robust_score",
                "validation_details",
            ]
        ],
        left_on="candidate_mirna",
        right_on="standard_feature_id",
        how="left",
    ).rename(
        columns={
            "discovery_logFC": "mirna_discovery_logFC",
            "discovery_adj_p": "mirna_discovery_adj_p",
            "same_direction_fdr_sig_count": "mirna_validation_fdr_support_count",
            "robust_score": "mirna_robust_score",
            "validation_details": "mirna_validation_details",
        }
    )
    axes = axes.merge(
        mrna[
            [
                "standard_feature_id",
                "target_gene",
                "discovery_logFC",
                "discovery_adj_p",
                "mrna_direction",
                "same_direction_fdr_sig_count",
                "robust_score",
                "validation_details",
            ]
        ],
        on="target_gene",
        how="left",
    ).rename(
        columns={
            "standard_feature_id_y": "target_gene_symbol",
            "discovery_logFC": "mrna_discovery_logFC",
            "discovery_adj_p": "mrna_discovery_adj_p",
            "same_direction_fdr_sig_count": "mrna_validation_fdr_support_count",
            "robust_score": "mrna_robust_score",
            "validation_details": "mrna_validation_details",
        }
    )
    axes = axes.drop(columns=[c for c in ["standard_feature_id_x"] if c in axes.columns])

    axes["axis"] = axes["candidate_mirna"] + " -> " + axes["target_gene"]
    axes["axis_score"] = (
        axes["best_support_rank"]
        + axes["evidence_count"].clip(upper=20)
        + axes["mirna_robust_score"].fillna(0)
        + axes["mrna_robust_score"].fillna(0)
    )
    axes["match_priority"] = axes["match_type"].map({"exact": 0, "arm_agnostic": 1}).fillna(9)
    axes = axes.sort_values(
        ["match_priority", "axis_score", "best_support_rank", "evidence_count"],
        ascending=[True, False, False, False],
    )

    axes.to_csv(OUTPUT_DIR / "robust_mirna_mrna_negative_axes_mirtarbase.csv", index=False, encoding="utf-8-sig")

    exact_axes = axes[axes["match_type"] == "exact"].copy()
    exact_axes.to_csv(OUTPUT_DIR / "robust_mirna_mrna_negative_axes_mirtarbase_exact.csv", index=False, encoding="utf-8-sig")

    top_axes = axes.head(100).copy()
    top_axes.to_csv(OUTPUT_DIR / "top100_robust_mirna_mrna_axes.csv", index=False, encoding="utf-8-sig")

    mirna_summary = (
        axes.groupby("candidate_mirna", as_index=False)
        .agg(
            target_count=("target_gene", "nunique"),
            exact_target_count=("match_type", lambda x: int((x == "exact").sum())),
            mean_axis_score=("axis_score", "mean"),
            top_targets=("target_gene", lambda x: ";".join(list(dict.fromkeys(x.astype(str)))[:20])),
        )
        .sort_values(["target_count", "mean_axis_score"], ascending=[False, False])
    )
    mirna_summary.to_csv(OUTPUT_DIR / "mirna_axis_summary.csv", index=False, encoding="utf-8-sig")

    gene_summary = (
        axes.groupby("target_gene", as_index=False)
        .agg(
            regulator_mirna_count=("candidate_mirna", "nunique"),
            regulators=("candidate_mirna", lambda x: ";".join(sorted(set(map(str, x))))),
            mean_axis_score=("axis_score", "mean"),
        )
        .sort_values(["regulator_mirna_count", "mean_axis_score"], ascending=[False, False])
    )
    gene_summary.to_csv(OUTPUT_DIR / "target_gene_axis_summary.csv", index=False, encoding="utf-8-sig")

    qc = pd.DataFrame(
        [
            {
                "robust_mirna_candidates": len(mirna),
                "robust_mrna_candidates": len(mrna),
                "mirtarbase_human_rows": len(mti),
                "matched_mti_rows_for_candidate_mirnas": len(matched_mti),
                "negative_direction_axes": len(axes),
                "exact_negative_direction_axes": len(exact_axes),
                "unique_candidate_mirnas_in_axes": axes["candidate_mirna"].nunique(),
                "unique_target_genes_in_axes": axes["target_gene"].nunique(),
                "target_database": "miRTarBase 2025 v10 hsa_MTI.csv",
                "negative_direction_rule": "miRNA and mRNA discovery logFC have opposite signs",
            }
        ]
    )
    qc.to_csv(OUTPUT_DIR / "mirna_mrna_axis_qc.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(OUTPUT_DIR / "mirna_mrna_axes_summary.xlsx", engine="openpyxl") as writer:
        qc.to_excel(writer, sheet_name="qc", index=False)
        top_axes.drop(columns=["match_priority"]).to_excel(writer, sheet_name="top100_axes", index=False)
        axes.drop(columns=["match_priority"]).to_excel(writer, sheet_name="all_axes", index=False)
        mirna_summary.to_excel(writer, sheet_name="mirna_summary", index=False)
        gene_summary.to_excel(writer, sheet_name="target_gene_summary", index=False)

    print(qc.to_string(index=False))
    print("\nTop axes:")
    print(top_axes[["axis", "match_type", "axis_score", "mirna_discovery_logFC", "mrna_discovery_logFC", "support_types"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
