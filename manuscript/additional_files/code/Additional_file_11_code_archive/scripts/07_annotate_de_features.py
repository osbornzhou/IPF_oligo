#!/usr/bin/env python
"""
Annotate differential-expression feature IDs to gene symbols or miRNA names.

Inputs:
  - results/differential_expression/*_de_results.csv
  - GEO family.soft.gz platform annotation files under data_raw/GEO

Outputs:
  - results/feature_annotation/feature_annotation_map.csv
  - results/feature_annotation/feature_annotation_qc.csv
  - results/differential_expression_annotated/*_de_results_annotated.csv
  - results/differential_expression_annotated/*_de_results_annotated_gene_or_mirna_level.csv
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
GEO_DIR = PROJECT_DIR / "data_raw" / "GEO"
DE_DIR = PROJECT_DIR / "results" / "differential_expression"
ANNOTATED_DIR = PROJECT_DIR / "results" / "differential_expression_annotated"
FEATURE_ANNOTATION_DIR = PROJECT_DIR / "results" / "feature_annotation"

ANNOTATED_DIR.mkdir(parents=True, exist_ok=True)
FEATURE_ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)


SERIES_PLATFORM = {
    "GSE110147": "GPL6244",
    "GSE150910": "RNAseq_gene_symbol",
    "GSE21394": "GPL8936",
    "GSE27430": "GPL8227",
    "GSE32537": "GPL6244",
    "GSE32538": "GPL8786",
    "GSE53845": "GPL6480",
    "GSE92592": "RNAseq_gene_symbol_chr_suffix",
}


def clean_text(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "---", "null", "none"}:
        return ""
    return text


def split_unique(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        value = clean_text(value)
        if not value or value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def normalize_mirna_name(feature_id: str, fallback: str = "") -> str:
    fid = clean_text(feature_id)
    name = clean_text(fallback) or fid
    name = re.sub(r"_st$", "", name)
    name = re.sub(r"_x_st$", "", name)
    name = name.replace("hsa-mir-", "hsa-miR-")
    return name


def parse_gene_assignment(gene_assignment: str) -> tuple[str, str, str]:
    """Parse Affymetrix GPL6244 gene_assignment field."""
    text = clean_text(gene_assignment)
    if not text:
        return "", "", ""
    symbols = []
    entrez_ids = []
    descriptions = []
    for item in text.split("///"):
        parts = [clean_text(p) for p in item.split("//")]
        if len(parts) >= 2:
            symbol = parts[1]
            if symbol and not symbol.startswith("ENST") and not symbol.startswith("NONHS"):
                symbols.append(symbol)
        if len(parts) >= 3:
            descriptions.append(parts[2])
        if len(parts) >= 5:
            entrez_ids.append(parts[4])
    symbols = split_unique(symbols)
    entrez_ids = split_unique(entrez_ids)
    descriptions = split_unique(descriptions)
    return (
        symbols[0] if symbols else "",
        ";".join(symbols),
        ";".join(entrez_ids),
    )


def read_platform_table(series_id: str) -> tuple[str, pd.DataFrame]:
    soft_path = GEO_DIR / series_id / "soft" / f"{series_id}_family.soft.gz"
    if not soft_path.exists():
        return "", pd.DataFrame()

    platform_id = ""
    rows = []
    header = None
    in_table = False
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^PLATFORM"):
                platform_id = line.split("=", 1)[1].strip()
                continue
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(handle).rstrip("\n").split("\t")
                continue
            if in_table:
                if line.startswith("!platform_table_end"):
                    break
                values = line.split("\t")
                if header and len(values) < len(header):
                    values.extend([""] * (len(header) - len(values)))
                rows.append(dict(zip(header or [], values)))

    df = pd.DataFrame(rows)
    return platform_id, df


def build_platform_mapping(series_id: str, feature_ids: set[str]) -> pd.DataFrame:
    platform_declared = SERIES_PLATFORM.get(series_id, "")

    if platform_declared == "RNAseq_gene_symbol":
        records = []
        for fid in sorted(feature_ids):
            symbol = clean_text(fid)
            records.append(
                {
                    "series_id": series_id,
                    "platform": platform_declared,
                    "feature_id": fid,
                    "standard_feature_id": symbol,
                    "gene_symbol": symbol,
                    "mirna_name": "",
                    "entrez_id": "",
                    "all_mapped_symbols": symbol,
                    "annotation_source": "feature_id_gene_symbol",
                    "is_annotated": bool(symbol),
                    "is_multi_mapped": False,
                }
            )
        return pd.DataFrame(records)

    if platform_declared == "RNAseq_gene_symbol_chr_suffix":
        records = []
        for fid in sorted(feature_ids):
            symbol = re.sub(r"\.chr[0-9XYM]+$", "", clean_text(fid), flags=re.IGNORECASE)
            records.append(
                {
                    "series_id": series_id,
                    "platform": platform_declared,
                    "feature_id": fid,
                    "standard_feature_id": symbol,
                    "gene_symbol": symbol,
                    "mirna_name": "",
                    "entrez_id": "",
                    "all_mapped_symbols": symbol,
                    "annotation_source": "feature_id_gene_symbol_chr_suffix_removed",
                    "is_annotated": bool(symbol),
                    "is_multi_mapped": False,
                }
            )
        return pd.DataFrame(records)

    platform_id, platform_df = read_platform_table(series_id)
    platform = platform_id or platform_declared
    if platform_df.empty or "ID" not in platform_df.columns:
        return pd.DataFrame(
            {
                "series_id": [series_id for _ in feature_ids],
                "platform": [platform for _ in feature_ids],
                "feature_id": sorted(feature_ids),
                "standard_feature_id": ["" for _ in feature_ids],
                "gene_symbol": ["" for _ in feature_ids],
                "mirna_name": ["" for _ in feature_ids],
                "entrez_id": ["" for _ in feature_ids],
                "all_mapped_symbols": ["" for _ in feature_ids],
                "annotation_source": ["missing_platform_table" for _ in feature_ids],
                "is_annotated": [False for _ in feature_ids],
                "is_multi_mapped": [False for _ in feature_ids],
            }
        )

    platform_df = platform_df.drop_duplicates(subset=["ID"], keep="first")
    platform_df = platform_df[platform_df["ID"].isin(feature_ids)].copy()
    records = []

    for row in platform_df.to_dict(orient="records"):
        fid = clean_text(row.get("ID"))
        gene_symbol = ""
        mirna_name = ""
        entrez_id = ""
        all_symbols = ""
        source = f"{platform}_platform_table"

        if platform == "GPL6244":
            gene_symbol, all_symbols, entrez_id = parse_gene_assignment(row.get("gene_assignment", ""))
        elif platform == "GPL6480":
            gene_symbol = clean_text(row.get("GENE_SYMBOL"))
            entrez_id = clean_text(row.get("GENE"))
            all_symbols = gene_symbol
        elif platform == "GPL8936":
            mirna_name = normalize_mirna_name(fid, row.get("miRNA_ID", ""))
            all_symbols = mirna_name
        elif platform == "GPL8227":
            mirna_name = normalize_mirna_name(fid, row.get("miRNA_ID", ""))
            all_symbols = mirna_name
        elif platform == "GPL8786":
            mirna_name = normalize_mirna_name(fid, fid)
            all_symbols = clean_text(row.get("miRNA_ID_LIST")) or mirna_name
        else:
            gene_symbol = clean_text(row.get("GENE_SYMBOL"))
            mirna_name = clean_text(row.get("miRNA_ID"))
            all_symbols = gene_symbol or mirna_name

        standard = mirna_name or gene_symbol
        records.append(
            {
                "series_id": series_id,
                "platform": platform,
                "feature_id": fid,
                "standard_feature_id": standard,
                "gene_symbol": gene_symbol,
                "mirna_name": mirna_name,
                "entrez_id": entrez_id,
                "all_mapped_symbols": all_symbols,
                "annotation_source": source,
                "is_annotated": bool(standard),
                "is_multi_mapped": ";" in all_symbols,
            }
        )

    mapped = pd.DataFrame(records)
    missing_ids = sorted(feature_ids - set(mapped["feature_id"] if not mapped.empty else []))
    if missing_ids:
        missing = pd.DataFrame(
            {
                "series_id": [series_id for _ in missing_ids],
                "platform": [platform for _ in missing_ids],
                "feature_id": missing_ids,
                "standard_feature_id": ["" for _ in missing_ids],
                "gene_symbol": ["" for _ in missing_ids],
                "mirna_name": ["" for _ in missing_ids],
                "entrez_id": ["" for _ in missing_ids],
                "all_mapped_symbols": ["" for _ in missing_ids],
                "annotation_source": ["feature_not_found_in_platform_table" for _ in missing_ids],
                "is_annotated": [False for _ in missing_ids],
                "is_multi_mapped": [False for _ in missing_ids],
            }
        )
        mapped = pd.concat([mapped, missing], ignore_index=True)
    return mapped


def collapse_to_standard_feature(df: pd.DataFrame) -> pd.DataFrame:
    keep = df[df["is_annotated"] & df["standard_feature_id"].astype(bool)].copy()
    if keep.empty:
        return keep
    keep["abs_logFC"] = keep["logFC"].abs()
    keep = keep.sort_values(
        ["series_id", "standard_feature_id", "adj.P.Val", "abs_logFC"],
        ascending=[True, True, True, False],
    )
    keep = keep.drop_duplicates(subset=["series_id", "standard_feature_id"], keep="first")
    keep = keep.drop(columns=["abs_logFC"])
    return keep


def main() -> None:
    result_paths = sorted(DE_DIR.glob("GSE*_IPF_vs_Control_de_results.csv"))
    if not result_paths:
        raise FileNotFoundError(f"No differential-expression result files found in {DE_DIR}")

    all_mappings = []
    qc_rows = []

    for result_path in result_paths:
        match = re.match(r"(GSE\d+)_IPF_vs_Control_de_results\.csv", result_path.name)
        if not match:
            continue
        series_id = match.group(1)
        de = pd.read_csv(result_path)
        de["series_id"] = de["series_id"].astype(str)
        de["feature_id"] = de["feature_id"].astype(str)
        feature_ids = set(de["feature_id"].astype(str))
        mapping = build_platform_mapping(series_id, feature_ids)
        all_mappings.append(mapping)

        annotated = de.merge(mapping, on=["series_id", "feature_id"], how="left")
        for col in [
            "platform",
            "standard_feature_id",
            "gene_symbol",
            "mirna_name",
            "entrez_id",
            "all_mapped_symbols",
            "annotation_source",
        ]:
            annotated[col] = annotated[col].fillna("")
        annotated["is_annotated"] = annotated["is_annotated"].fillna(False).astype(bool)
        annotated["is_multi_mapped"] = annotated["is_multi_mapped"].fillna(False).astype(bool)

        annotated_path = ANNOTATED_DIR / f"{series_id}_IPF_vs_Control_de_results_annotated.csv"
        annotated_sig_path = ANNOTATED_DIR / f"{series_id}_IPF_vs_Control_de_significant_annotated_fdr0.05_logfc1.csv"
        collapsed_path = ANNOTATED_DIR / f"{series_id}_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv"
        collapsed_sig_path = ANNOTATED_DIR / f"{series_id}_IPF_vs_Control_de_significant_annotated_gene_or_mirna_level_fdr0.05_logfc1.csv"

        annotated.to_csv(annotated_path, index=False, encoding="utf-8-sig")
        annotated[
            (annotated["adj.P.Val"] < 0.05)
            & (annotated["logFC"].abs() >= 1)
            & (annotated["is_annotated"])
        ].to_csv(annotated_sig_path, index=False, encoding="utf-8-sig")

        collapsed = collapse_to_standard_feature(annotated)
        collapsed.to_csv(collapsed_path, index=False, encoding="utf-8-sig")
        collapsed[
            (collapsed["adj.P.Val"] < 0.05) & (collapsed["logFC"].abs() >= 1)
        ].to_csv(collapsed_sig_path, index=False, encoding="utf-8-sig")

        qc_rows.append(
            {
                "series_id": series_id,
                "platform": mapping["platform"].iloc[0] if not mapping.empty else "",
                "input_features": len(feature_ids),
                "annotated_features": int(annotated["is_annotated"].sum()),
                "unannotated_features": int((~annotated["is_annotated"]).sum()),
                "annotation_rate": round(float(annotated["is_annotated"].mean()), 6),
                "multi_mapped_features": int(annotated["is_multi_mapped"].sum()),
                "unique_standard_features": int(annotated.loc[annotated["is_annotated"], "standard_feature_id"].nunique()),
                "significant_annotated_features_fdr0.05_logfc1": int(
                    (
                        (annotated["adj.P.Val"] < 0.05)
                        & (annotated["logFC"].abs() >= 1)
                        & (annotated["is_annotated"])
                    ).sum()
                ),
                "collapsed_significant_features_fdr0.05_logfc1": int(
                    ((collapsed["adj.P.Val"] < 0.05) & (collapsed["logFC"].abs() >= 1)).sum()
                ),
                "annotated_result_path": str(annotated_path),
                "collapsed_result_path": str(collapsed_path),
            }
        )

    feature_map = pd.concat(all_mappings, ignore_index=True)
    feature_map.to_csv(FEATURE_ANNOTATION_DIR / "feature_annotation_map.csv", index=False, encoding="utf-8-sig")

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(FEATURE_ANNOTATION_DIR / "feature_annotation_qc.csv", index=False, encoding="utf-8-sig")

    significant_paths = sorted(ANNOTATED_DIR.glob("*_de_significant_annotated_gene_or_mirna_level_fdr0.05_logfc1.csv"))
    combined = []
    for path in significant_paths:
        df = pd.read_csv(path)
        if not df.empty:
            combined.append(df)
    if combined:
        pd.concat(combined, ignore_index=True).to_csv(
            ANNOTATED_DIR / "all_significant_annotated_gene_or_mirna_level_fdr0.05_logfc1.csv",
            index=False,
            encoding="utf-8-sig",
        )
    else:
        pd.DataFrame().to_csv(
            ANNOTATED_DIR / "all_significant_annotated_gene_or_mirna_level_fdr0.05_logfc1.csv",
            index=False,
            encoding="utf-8-sig",
        )

    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
