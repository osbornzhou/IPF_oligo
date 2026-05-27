#!/usr/bin/env python
"""
Select robust IPF mRNA and miRNA candidates by discovery-validation agreement.

Strategy:
  - Discovery mRNA dataset: GSE32537
  - Discovery miRNA dataset: GSE32538
  - Validation mRNA datasets: GSE110147, GSE150910, GSE53845, GSE92592
  - Validation miRNA datasets: GSE21394, GSE27430

Discovery candidates must pass FDR < 0.05 and |logFC| >= 1.
Validation support is counted only when the validation logFC has the same
direction as discovery.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
ANNOTATED_DIR = PROJECT_DIR / "results" / "differential_expression_annotated"
OUTPUT_DIR = PROJECT_DIR / "results" / "robust_candidates"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DISCOVERY = {
    "bulk mRNA": "GSE32537",
    "miRNA": "GSE32538",
}

VALIDATION = {
    "bulk mRNA": ["GSE110147", "GSE150910", "GSE53845", "GSE92592"],
    "miRNA": ["GSE21394", "GSE27430"],
}

FDR_CUTOFF = 0.05
LOGFC_CUTOFF = 1.0
NOMINAL_P_CUTOFF = 0.05


def load_de(series_id: str) -> pd.DataFrame:
    path = ANNOTATED_DIR / f"{series_id}_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    df = df[df["is_annotated"].astype(str).str.lower().isin(["true", "1", "yes"])]
    df["standard_feature_id"] = df["standard_feature_id"].astype(str)
    df = df[df["standard_feature_id"].notna() & (df["standard_feature_id"] != "")]
    df["direction"] = df["logFC"].apply(lambda x: "up" if x > 0 else ("down" if x < 0 else "zero"))
    df["fdr_sig"] = (df["adj.P.Val"] < FDR_CUTOFF) & (df["logFC"].abs() >= LOGFC_CUTOFF)
    df["nominal_sig"] = (df["P.Value"] < NOMINAL_P_CUTOFF) & (df["logFC"].abs() >= LOGFC_CUTOFF)
    return df


def signed_direction(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def select_for_data_type(data_type: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    discovery_id = DISCOVERY[data_type]
    validation_ids = VALIDATION[data_type]

    discovery = load_de(discovery_id)
    validations = {sid: load_de(sid) for sid in validation_ids}

    discovery_sig = discovery[
        (discovery["adj.P.Val"] < FDR_CUTOFF)
        & (discovery["logFC"].abs() >= LOGFC_CUTOFF)
    ].copy()

    rows = []
    matrix_rows = []
    for _, drow in discovery_sig.iterrows():
        feature = drow["standard_feature_id"]
        discovery_direction = signed_direction(drow["logFC"])

        validation_present = 0
        same_direction = 0
        opposite_direction = 0
        same_direction_fdr = 0
        same_direction_nominal = 0
        validation_details = []

        out = {
            "data_type": data_type,
            "standard_feature_id": feature,
            "gene_symbol": drow.get("gene_symbol", ""),
            "mirna_name": drow.get("mirna_name", ""),
            "discovery_series_id": discovery_id,
            "discovery_logFC": drow["logFC"],
            "discovery_adj_p": drow["adj.P.Val"],
            "discovery_p_value": drow["P.Value"],
            "discovery_direction": drow["direction"],
        }

        matrix = {
            "data_type": data_type,
            "standard_feature_id": feature,
            f"{discovery_id}_logFC": drow["logFC"],
            f"{discovery_id}_adj_p": drow["adj.P.Val"],
            f"{discovery_id}_direction": drow["direction"],
            f"{discovery_id}_fdr_sig": True,
        }

        for sid, vdf in validations.items():
            vmatch = vdf[vdf["standard_feature_id"] == feature]
            if vmatch.empty:
                out[f"{sid}_present"] = False
                out[f"{sid}_logFC"] = math.nan
                out[f"{sid}_adj_p"] = math.nan
                out[f"{sid}_p_value"] = math.nan
                out[f"{sid}_same_direction"] = False
                out[f"{sid}_same_direction_fdr_sig"] = False
                out[f"{sid}_same_direction_nominal_sig"] = False
                matrix[f"{sid}_logFC"] = math.nan
                matrix[f"{sid}_adj_p"] = math.nan
                matrix[f"{sid}_direction"] = "missing"
                matrix[f"{sid}_fdr_sig"] = False
                continue

            # Gene/miRNA-level files are already collapsed; keep first row if duplicated.
            vrow = vmatch.iloc[0]
            validation_present += 1
            v_direction = signed_direction(vrow["logFC"])
            same = v_direction == discovery_direction and v_direction != 0
            opposite = v_direction == -discovery_direction and v_direction != 0
            fdr_same = same and bool(vrow["fdr_sig"])
            nominal_same = same and bool(vrow["nominal_sig"])

            same_direction += int(same)
            opposite_direction += int(opposite)
            same_direction_fdr += int(fdr_same)
            same_direction_nominal += int(nominal_same)

            detail = f"{sid}:{vrow['logFC']:.3g}/FDR={vrow['adj.P.Val']:.3g}"
            validation_details.append(detail)

            out[f"{sid}_present"] = True
            out[f"{sid}_logFC"] = vrow["logFC"]
            out[f"{sid}_adj_p"] = vrow["adj.P.Val"]
            out[f"{sid}_p_value"] = vrow["P.Value"]
            out[f"{sid}_same_direction"] = same
            out[f"{sid}_same_direction_fdr_sig"] = fdr_same
            out[f"{sid}_same_direction_nominal_sig"] = nominal_same

            matrix[f"{sid}_logFC"] = vrow["logFC"]
            matrix[f"{sid}_adj_p"] = vrow["adj.P.Val"]
            matrix[f"{sid}_direction"] = vrow["direction"]
            matrix[f"{sid}_fdr_sig"] = bool(vrow["fdr_sig"])

        validation_available = max(validation_present, 1)
        direction_consistency_fraction = same_direction / validation_available

        if data_type == "bulk mRNA":
            strict = same_direction_fdr >= 2
            relaxed = (same_direction_fdr >= 1 and same_direction_nominal >= 2) or same_direction_fdr >= 2
        else:
            strict = same_direction_fdr >= 1
            relaxed = same_direction_fdr >= 1 or same_direction_nominal >= 1

        score = (
            same_direction_fdr * 4
            + same_direction_nominal * 2
            + same_direction
            - opposite_direction * 2
            + direction_consistency_fraction
        )

        out.update(
            {
                "validation_datasets_total": len(validation_ids),
                "validation_datasets_present": validation_present,
                "same_direction_validation_count": same_direction,
                "opposite_direction_validation_count": opposite_direction,
                "same_direction_fdr_sig_count": same_direction_fdr,
                "same_direction_nominal_sig_count": same_direction_nominal,
                "direction_consistency_fraction": round(direction_consistency_fraction, 6),
                "robust_strict": strict,
                "robust_relaxed": relaxed,
                "robust_score": round(score, 6),
                "validation_details": "; ".join(validation_details),
            }
        )
        matrix.update(
            {
                "validation_datasets_present": validation_present,
                "same_direction_validation_count": same_direction,
                "same_direction_fdr_sig_count": same_direction_fdr,
                "robust_strict": strict,
                "robust_relaxed": relaxed,
            }
        )
        rows.append(out)
        matrix_rows.append(matrix)

    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values(
            ["robust_strict", "robust_relaxed", "same_direction_fdr_sig_count", "robust_score", "discovery_adj_p"],
            ascending=[False, False, False, False, True],
        )

    matrix = pd.DataFrame(matrix_rows)
    strict_df = summary[summary["robust_strict"]].copy() if not summary.empty else summary.copy()
    relaxed_df = summary[summary["robust_relaxed"]].copy() if not summary.empty else summary.copy()
    return summary, strict_df, relaxed_df, matrix


def main() -> None:
    qc_rows = []
    all_summary = []
    all_matrix = []

    for data_type in ["bulk mRNA", "miRNA"]:
        summary, strict_df, relaxed_df, matrix = select_for_data_type(data_type)
        label = "mrna" if data_type == "bulk mRNA" else "mirna"

        summary.to_csv(OUTPUT_DIR / f"{label}_discovery_validation_summary.csv", index=False, encoding="utf-8-sig")
        strict_df.to_csv(OUTPUT_DIR / f"robust_{label}_candidates_strict.csv", index=False, encoding="utf-8-sig")
        relaxed_df.to_csv(OUTPUT_DIR / f"robust_{label}_candidates_relaxed.csv", index=False, encoding="utf-8-sig")
        matrix.to_csv(OUTPUT_DIR / f"{label}_direction_matrix.csv", index=False, encoding="utf-8-sig")

        all_summary.append(summary)
        all_matrix.append(matrix)

        qc_rows.append(
            {
                "data_type": data_type,
                "discovery_dataset": DISCOVERY[data_type],
                "validation_datasets": ";".join(VALIDATION[data_type]),
                "discovery_significant_features": len(summary),
                "robust_strict_candidates": len(strict_df),
                "robust_relaxed_candidates": len(relaxed_df),
                "strict_rule": "mRNA: >=2 same-direction validation datasets with FDR<0.05 and |logFC|>=1; miRNA: >=1",
                "relaxed_rule": "mRNA: >=1 FDR same-direction plus >=2 nominal same-direction, or strict; miRNA: >=1 FDR or nominal same-direction",
            }
        )

    pd.concat(all_summary, ignore_index=True).to_csv(
        OUTPUT_DIR / "all_discovery_validation_summary.csv", index=False, encoding="utf-8-sig"
    )
    pd.concat(all_matrix, ignore_index=True).to_csv(
        OUTPUT_DIR / "all_direction_matrix.csv", index=False, encoding="utf-8-sig"
    )

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUTPUT_DIR / "robust_candidate_qc.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(OUTPUT_DIR / "robust_candidates_summary.xlsx", engine="openpyxl") as writer:
        qc.to_excel(writer, sheet_name="qc", index=False)
        pd.concat(all_summary, ignore_index=True).to_excel(writer, sheet_name="all_summary", index=False)
        pd.read_csv(OUTPUT_DIR / "robust_mrna_candidates_strict.csv").to_excel(writer, sheet_name="mrna_strict", index=False)
        pd.read_csv(OUTPUT_DIR / "robust_mrna_candidates_relaxed.csv").to_excel(writer, sheet_name="mrna_relaxed", index=False)
        pd.read_csv(OUTPUT_DIR / "robust_mirna_candidates_strict.csv").to_excel(writer, sheet_name="mirna_strict", index=False)
        pd.read_csv(OUTPUT_DIR / "robust_mirna_candidates_relaxed.csv").to_excel(writer, sheet_name="mirna_relaxed", index=False)

    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
