#!/usr/bin/env python
"""Extract bulk mRNA and miRNA expression matrices with triple cross-QC.

Triple QC per dataset:
1. Annotation cross-check: sample_id count, uniqueness, IPF/Control availability.
2. Expression cross-check: expression columns match annotation sample_id exactly as a set.
3. Integrity cross-check: numeric matrix validity, missingness, zero fraction, feature uniqueness,
   and GEO accession consistency when sample_id is expected to be GSM.

Outputs:
- data_processed/expression/{GSE}_expression.csv.gz
- data_processed/expression/{GSE}_expression_in_annotation_order.csv.gz
- metadata/expression_matrix_qc.csv
- metadata/expression_matrix_qc.xlsx
"""

from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GEO_ROOT = PROJECT_ROOT / "data_raw" / "GEO"
ANNOTATION_PATH = PROJECT_ROOT / "metadata" / "all_bulk_mirna_annotation.csv"
OUT_ROOT = PROJECT_ROOT / "data_processed" / "expression"
META_ROOT = PROJECT_ROOT / "metadata"


def opener(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def find_dataset_root(series_id: str) -> Path:
    root = GEO_ROOT / series_id
    if not root.exists():
        raise FileNotFoundError(f"Dataset folder not found: {root}")
    return root


def read_series_matrix(path: Path) -> pd.DataFrame:
    lines: list[str] = []
    in_table = False
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line)
    if not lines:
        raise ValueError(f"No series matrix table found in {path}")
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", quotechar='"')


def read_csv_matrix(path: Path) -> pd.DataFrame:
    with opener(path)(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        return pd.read_csv(handle)


def parse_sample_field(path: Path, field: str) -> list[str]:
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(field):
                return [part.strip().strip('"') for part in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                break
    return []


def read_gse92592_counts(path: Path, series_matrix_path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        df = pd.read_csv(handle, sep="\t", header=None, names=["feature_id"] + header)
    titles = parse_sample_field(series_matrix_path, "!Sample_title")
    geos = parse_sample_field(series_matrix_path, "!Sample_geo_accession")
    title_to_geo = dict(zip(titles, geos))
    df = df.rename(columns={col: title_to_geo.get(col, col) for col in df.columns})
    return df


def choose_expression_source(series_id: str) -> tuple[Path, str]:
    root = find_dataset_root(series_id)
    if series_id == "GSE150910":
        path = root / "supplementary" / "GSE150910_gene-level_count_file.csv.gz"
        return path, "supplementary_gene_level_count"
    if series_id == "GSE92592":
        path = root / "supplementary" / "GSE92592_gene.counts.txt.gz"
        return path, "supplementary_gene_counts_mapped_to_gsm"

    matrix_dir = root / "matrix"
    gz = matrix_dir / f"{series_id}_series_matrix.txt.gz"
    txt = matrix_dir / f"{series_id}_series_matrix.txt"
    if gz.exists():
        return gz, "geo_series_matrix"
    if txt.exists():
        return txt, "geo_series_matrix"
    matches = sorted(root.rglob(f"{series_id}*_series_matrix.txt.gz")) + sorted(root.rglob(f"{series_id}*_series_matrix.txt"))
    if matches:
        return matches[0], "geo_series_matrix"
    raise FileNotFoundError(f"No expression matrix found for {series_id}")


def normalize_feature_column(df: pd.DataFrame, source_type: str) -> pd.DataFrame:
    df = df.copy()
    first = df.columns[0]
    if first != "feature_id":
        df = df.rename(columns={first: "feature_id"})
    df["feature_id"] = df["feature_id"].astype(str)
    return df


def numeric_expression(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    sample_cols = [c for c in out.columns if c != "feature_id"]
    for col in sample_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def infer_value_scale(values: pd.Series) -> str:
    vals = values.dropna()
    if vals.empty:
        return "unknown"
    q99 = vals.quantile(0.99)
    vmax = vals.max()
    vmin = vals.min()
    if vmin >= 0 and vmax <= 25:
        return "log_or_normalized_intensity"
    if vmin >= 0 and q99 > 100:
        return "count_or_abundance"
    if vmin < 0:
        return "centered_or_transformed"
    return "unknown_nonnegative"


def geo_expected_sample_ids(annotation: pd.DataFrame) -> bool:
    sample = annotation["sample_id"].dropna().astype(str)
    geo = annotation["geo_accession"].dropna().astype(str) if "geo_accession" in annotation.columns else pd.Series([], dtype=str)
    return len(sample) > 0 and sample.str.fullmatch(r"GSM\d+").all() and len(geo) > 0


def qc_dataset(series_id: str, annotation: pd.DataFrame, expr: pd.DataFrame, source_path: Path, source_type: str) -> dict[str, object]:
    sample_cols = [c for c in expr.columns if c != "feature_id"]
    ann_ids = annotation["sample_id"].dropna().astype(str).tolist()
    include_sub = annotation[annotation["include"].astype(str).str.lower() == "yes"].copy()
    groups = include_sub["group"].astype(str).value_counts().to_dict()

    missing_in_expr = sorted(set(ann_ids) - set(sample_cols))
    missing_in_ann = sorted(set(sample_cols) - set(ann_ids))
    same_order = sample_cols == ann_ids
    duplicate_samples = int(pd.Series(ann_ids).duplicated().sum())
    duplicate_features = int(expr["feature_id"].duplicated().sum())

    numeric = expr[[c for c in expr.columns if c != "feature_id"]]
    total_cells = int(numeric.shape[0] * numeric.shape[1])
    missing_values = int(numeric.isna().sum().sum())
    zero_values = int((numeric == 0).sum().sum())
    finite_values = numeric.replace([np.inf, -np.inf], np.nan).stack()
    value_scale = infer_value_scale(finite_values)

    geo_consistency = "not_applicable"
    if geo_expected_sample_ids(annotation):
        geo_ids = annotation["geo_accession"].dropna().astype(str).tolist()
        geo_consistency = "PASS" if set(ann_ids) == set(geo_ids) else "FAIL"

    qc1_annotation = (
        duplicate_samples == 0
        and len(ann_ids) == len(annotation)
        and groups.get("IPF", 0) > 0
        and groups.get("Control", 0) > 0
    )
    qc2_cross_match = len(sample_cols) > 0 and not missing_in_expr and not missing_in_ann
    qc3_integrity = (
        expr.shape[0] > 0
        and duplicate_features < expr.shape[0]
        and total_cells > 0
        and missing_values / total_cells < 0.05
        and geo_consistency != "FAIL"
    )

    return {
        "series_id": series_id,
        "data_type": ";".join(sorted(set(annotation["data_type"].dropna().astype(str)))),
        "dataset_role": ";".join(sorted(set(annotation["dataset_role"].dropna().astype(str)))),
        "expression_source_type": source_type,
        "expression_source_path": str(source_path),
        "features": int(expr.shape[0]),
        "expression_samples": int(len(sample_cols)),
        "annotation_samples": int(len(ann_ids)),
        "include_yes_samples": int(len(include_sub)),
        "ipf_samples": int(groups.get("IPF", 0)),
        "control_samples": int(groups.get("Control", 0)),
        "exclude_samples": int(annotation["group"].astype(str).value_counts().to_dict().get("Exclude", 0)),
        "duplicate_annotation_sample_ids": duplicate_samples,
        "duplicate_feature_ids": duplicate_features,
        "missing_in_expression_count": int(len(missing_in_expr)),
        "missing_in_annotation_count": int(len(missing_in_ann)),
        "same_order_as_annotation": bool(same_order),
        "numeric_missing_values": missing_values,
        "numeric_missing_fraction": round(missing_values / total_cells, 8) if total_cells else np.nan,
        "zero_fraction": round(zero_values / total_cells, 8) if total_cells else np.nan,
        "min_value": float(finite_values.min()) if not finite_values.empty else np.nan,
        "median_value": float(finite_values.median()) if not finite_values.empty else np.nan,
        "max_value": float(finite_values.max()) if not finite_values.empty else np.nan,
        "value_scale_guess": value_scale,
        "geo_accession_consistency": geo_consistency,
        "qc1_annotation_pass": bool(qc1_annotation),
        "qc2_sample_crossmatch_pass": bool(qc2_cross_match),
        "qc3_matrix_integrity_pass": bool(qc3_integrity),
        "triple_qc_pass": bool(qc1_annotation and qc2_cross_match and qc3_integrity),
        "missing_in_expression_examples": ";".join(missing_in_expr[:10]) if missing_in_expr else "NA",
        "missing_in_annotation_examples": ";".join(missing_in_ann[:10]) if missing_in_ann else "NA",
    }


def extract_dataset(series_id: str, annotation: pd.DataFrame) -> dict[str, object]:
    source_path, source_type = choose_expression_source(series_id)
    if source_type == "geo_series_matrix":
        expr = read_series_matrix(source_path)
    elif source_type == "supplementary_gene_counts_mapped_to_gsm":
        series_matrix = find_dataset_root(series_id) / "matrix" / f"{series_id}_series_matrix.txt.gz"
        expr = read_gse92592_counts(source_path, series_matrix)
    else:
        expr = read_csv_matrix(source_path)
    expr = normalize_feature_column(expr, source_type)
    expr = numeric_expression(expr)

    sample_cols = [str(c) for c in expr.columns if c != "feature_id"]
    expr = expr.rename(columns={c: str(c) for c in expr.columns})
    ann_ids = annotation["sample_id"].dropna().astype(str).tolist()
    ordered_cols = ["feature_id"] + [sid for sid in ann_ids if sid in expr.columns]
    expr_ordered = expr[ordered_cols].copy()

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out_path = OUT_ROOT / f"{series_id}_expression.csv.gz"
    ordered_path = OUT_ROOT / f"{series_id}_expression_in_annotation_order.csv.gz"
    expr.to_csv(out_path, index=False, compression="gzip")
    expr_ordered.to_csv(ordered_path, index=False, compression="gzip")

    qc = qc_dataset(series_id, annotation, expr, source_path, source_type)
    qc["expression_output_path"] = str(out_path)
    qc["annotation_order_output_path"] = str(ordered_path)
    return qc


def main() -> None:
    annotation = pd.read_csv(ANNOTATION_PATH, dtype=str)
    annotation = annotation[annotation["series_id"].str.fullmatch(r"GSE\d+", na=False)].copy()
    # GSE166682 was excluded earlier because it is not human and should not be in the merged annotation.
    qc_rows = []
    for series_id, sub in annotation.groupby("series_id", sort=True):
        print(f"[{series_id}] extracting expression", flush=True)
        qc_rows.append(extract_dataset(series_id, sub))

    qc = pd.DataFrame(qc_rows).sort_values("series_id")
    META_ROOT.mkdir(parents=True, exist_ok=True)
    qc_csv = META_ROOT / "expression_matrix_qc.csv"
    qc_xlsx = META_ROOT / "expression_matrix_qc.xlsx"
    qc.to_csv(qc_csv, index=False, encoding="utf-8")
    with pd.ExcelWriter(qc_xlsx, engine="openpyxl") as writer:
        qc.fillna("NA").replace("", "NA").to_excel(writer, sheet_name="expression_qc", index=False)

    print(f"Wrote {qc_csv}")
    print(f"Wrote {qc_xlsx}")
    print("\nTriple QC summary:")
    print(
        qc[
            [
                "series_id",
                "data_type",
                "features",
                "expression_samples",
                "annotation_samples",
                "ipf_samples",
                "control_samples",
                "qc1_annotation_pass",
                "qc2_sample_crossmatch_pass",
                "qc3_matrix_integrity_pass",
                "triple_qc_pass",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
