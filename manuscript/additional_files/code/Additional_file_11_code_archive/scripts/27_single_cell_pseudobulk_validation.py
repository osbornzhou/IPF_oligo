#!/usr/bin/env python
"""
Donor-aware pseudobulk validation for core prioritized IPF candidates.

This script streams the same sparse single-cell matrices used for localization
and aggregates target-gene counts at donor/sample x broad-celltype resolution.
It is intentionally limited to core candidates to provide a review-ready
pseudobulk sensitivity layer without turning the manuscript into a full
single-cell differential-expression paper.
"""

from __future__ import annotations

import gzip
import math
import os
from collections import defaultdict
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
GEO_DIR = PROJECT_DIR / "data_raw" / "GEO"
OUT_DIR = PROJECT_DIR / "results" / "single_cell_pseudobulk"
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_DIR / "results" / "_mpl_config"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


CORE_GENES = ["SPP1", "COL14A1", "POSTN", "COL1A1", "COL3A1", "GPX3", "CD24", "PTGFRN"]
EXCLUDE_LABEL_PATTERNS = ["MULTIPLET", "OUTLIER", "CELLCYCLE", "MT-TRNAS"]


def normalize_gene(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().upper()
    return text if text and text != "NA" else ""


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", errors="replace") as handle:
        return [line.rstrip("\n").strip('"') for line in handle]


def parse_mtx_header(handle) -> tuple[int, int, int]:
    for line in handle:
        text = line.strip()
        if not text or text.startswith("%"):
            continue
        parts = text.split()
        if len(parts) >= 3:
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise RuntimeError("Matrix Market header not found")


def load_gse135893() -> dict:
    supp = GEO_DIR / "GSE135893" / "supplementary"
    genes = [normalize_gene(x) for x in read_gzip_lines(supp / "GSE135893_genes.tsv.gz")]
    barcodes = read_gzip_lines(supp / "GSE135893_barcodes.tsv.gz")
    meta = pd.read_csv(supp / "GSE135893_IPF_metadata.csv.gz")
    first_col = meta.columns[0]
    meta = meta.rename(columns={first_col: "cell_id"})
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["donor_id"] = meta["Sample_Name"].astype(str)
    meta["group"] = meta["Diagnosis"].astype(str).replace({"control": "Control", "IPF": "IPF"})
    meta["broad_celltype"] = meta["population"].astype(str)
    meta["fine_celltype"] = meta["celltype"].astype(str)
    meta["library_size"] = pd.to_numeric(meta["nCount_RNA"], errors="coerce")
    return {
        "series_id": "GSE135893",
        "matrix_path": supp / "GSE135893_matrix.mtx.gz",
        "genes": genes,
        "barcodes": barcodes,
        "metadata": meta,
    }


def load_gse136831() -> dict:
    supp = GEO_DIR / "GSE136831" / "supplementary"
    gene_df = pd.read_csv(supp / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    genes = [normalize_gene(x) for x in gene_df["HGNC_EnsemblAlt_GeneID"]]
    barcodes = read_gzip_lines(supp / "GSE136831_AllCells.cellBarcodes.txt.gz")
    meta = pd.read_csv(supp / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    meta = meta.rename(columns={"CellBarcode_Identity": "cell_id"})
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["donor_id"] = meta["Subject_Identity"].astype(str)
    meta["group"] = meta["Disease_Identity"].astype(str)
    meta["broad_celltype"] = meta["CellType_Category"].astype(str)
    meta["fine_celltype"] = meta["Subclass_Cell_Identity"].astype(str)
    meta["library_size"] = pd.to_numeric(meta["nUMI"], errors="coerce")
    return {
        "series_id": "GSE136831",
        "matrix_path": supp / "GSE136831_RawCounts_Sparse.mtx.gz",
        "genes": genes,
        "barcodes": barcodes,
        "metadata": meta,
    }


def is_informative(row: pd.Series) -> bool:
    label = f"{row['broad_celltype']} {row['fine_celltype']}".upper()
    return not any(pattern in label for pattern in EXCLUDE_LABEL_PATTERNS)


def stream_pseudobulk(ds: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    series_id = ds["series_id"]
    meta = ds["metadata"].copy()
    meta = meta[meta["group"].isin(["IPF", "Control"])].copy()
    meta = meta[meta.apply(is_informative, axis=1)].copy()
    meta = meta.drop_duplicates(subset=["cell_id"], keep="first")
    meta_lookup = meta.set_index("cell_id")

    barcode_meta = []
    cell_counts = defaultdict(int)
    lib_sums = defaultdict(float)
    for barcode in ds["barcodes"]:
        if barcode not in meta_lookup.index:
            barcode_meta.append(None)
            continue
        row = meta_lookup.loc[barcode]
        lib = float(row["library_size"]) if pd.notna(row["library_size"]) and float(row["library_size"]) > 0 else math.nan
        if not math.isfinite(lib):
            barcode_meta.append(None)
            continue
        item = (str(row["donor_id"]), str(row["group"]), str(row["broad_celltype"]), lib)
        barcode_meta.append(item)
        key = item[:3]
        cell_counts[key] += 1
        lib_sums[key] += lib

    row_to_gene = {}
    target_found = set()
    for idx, gene in enumerate(ds["genes"], start=1):
        if gene in CORE_GENES:
            row_to_gene[idx] = gene
            target_found.add(gene)

    raw_sums = defaultdict(float)
    entries_seen = 0
    target_entries = 0
    with gzip.open(ds["matrix_path"], "rt", errors="replace") as handle:
        n_rows, n_cols, nnz_header = parse_mtx_header(handle)
        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            row_idx = int(parts[0])
            col_idx = int(parts[1])
            value = float(parts[2])
            entries_seen += 1
            gene = row_to_gene.get(row_idx)
            if gene is None:
                continue
            cell_meta = barcode_meta[col_idx - 1] if 0 <= col_idx - 1 < len(barcode_meta) else None
            if cell_meta is None:
                continue
            donor, group, broad, _lib = cell_meta
            raw_sums[(donor, group, broad, gene)] += value
            target_entries += 1

    rows = []
    for donor, group, broad in sorted(cell_counts):
        denom = lib_sums[(donor, group, broad)]
        for gene in sorted(target_found):
            raw = raw_sums.get((donor, group, broad, gene), 0.0)
            norm = raw / denom * 10000.0 if denom > 0 else np.nan
            rows.append(
                {
                    "series_id": series_id,
                    "donor_id": donor,
                    "group": group,
                    "broad_celltype": broad,
                    "gene_symbol": gene,
                    "cells": cell_counts[(donor, group, broad)],
                    "library_sum": denom,
                    "raw_count_sum": raw,
                    "norm_count_per_10k": norm,
                    "log1p_norm_count_per_10k": math.log1p(norm) if pd.notna(norm) else np.nan,
                }
            )

    donor_df = pd.DataFrame(rows)
    qc = pd.DataFrame(
        [
            {
                "series_id": series_id,
                "barcodes": len(ds["barcodes"]),
                "metadata_ipf_control_cells_after_label_filter": len(meta),
                "donors": meta["donor_id"].nunique(),
                "ipf_donors": meta.loc[meta["group"].eq("IPF"), "donor_id"].nunique(),
                "control_donors": meta.loc[meta["group"].eq("Control"), "donor_id"].nunique(),
                "mtx_rows": n_rows,
                "mtx_cols": n_cols,
                "mtx_nnz_header": nnz_header,
                "mtx_entries_streamed": entries_seen,
                "target_genes_requested": len(CORE_GENES),
                "target_genes_found": len(target_found),
                "target_entries_streamed": target_entries,
                "donor_celltype_gene_rows": len(donor_df),
                "qc1_donor_metadata_pass": meta["donor_id"].nunique() >= 4 and meta["group"].nunique() == 2,
                "qc2_matrix_integrity_pass": n_cols == len(ds["barcodes"]) and n_rows == len(ds["genes"]),
                "qc3_target_recovery_pass": len(target_found) >= len(CORE_GENES) - 1,
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_donor_metadata_pass"] & qc["qc2_matrix_integrity_pass"] & qc["qc3_target_recovery_pass"]
    return donor_df, qc


def test_pseudobulk(donor_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (series_id, broad, gene), sub in donor_df.groupby(["series_id", "broad_celltype", "gene_symbol"]):
        ipf = sub.loc[sub["group"].eq("IPF"), "log1p_norm_count_per_10k"].dropna()
        ctrl = sub.loc[sub["group"].eq("Control"), "log1p_norm_count_per_10k"].dropna()
        if len(ipf) < 2 or len(ctrl) < 2:
            continue
        stat = ttest_ind(ipf, ctrl, equal_var=False)
        rows.append(
            {
                "series_id": series_id,
                "broad_celltype": broad,
                "gene_symbol": gene,
                "ipf_donors": len(ipf),
                "control_donors": len(ctrl),
                "ipf_mean_log1p_norm": float(ipf.mean()),
                "control_mean_log1p_norm": float(ctrl.mean()),
                "ipf_minus_control_log1p_norm": float(ipf.mean() - ctrl.mean()),
                "welch_t": float(stat.statistic) if pd.notna(stat.statistic) else np.nan,
                "p_value": float(stat.pvalue) if pd.notna(stat.pvalue) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh"] = bh_fdr(out["p_value"])
        out = out.sort_values(["fdr_bh", "p_value", "series_id", "broad_celltype", "gene_symbol"])
    return out


def bh_fdr(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(out, index=p_values.index)
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    valid_idx = np.where(valid)[0]
    out[valid_idx[order]] = adj
    return pd.Series(out, index=p_values.index)


def summarize_core_contexts(stats: pd.DataFrame) -> pd.DataFrame:
    desired = {
        "SPP1": ["Immune", "Myeloid"],
        "COL14A1": ["Mesenchymal", "Stromal"],
        "POSTN": ["Mesenchymal", "Stromal", "Endothelial"],
        "COL1A1": ["Mesenchymal", "Stromal"],
        "COL3A1": ["Mesenchymal", "Stromal"],
        "GPX3": ["Mesenchymal", "Stromal"],
        "CD24": ["Epithelial"],
        "PTGFRN": ["Mesenchymal", "Stromal"],
    }
    rows = []
    for gene, celltypes in desired.items():
        sub = stats[(stats["gene_symbol"].eq(gene)) & (stats["broad_celltype"].isin(celltypes))].copy()
        if sub.empty:
            continue
        sub["abs_delta"] = sub["ipf_minus_control_log1p_norm"].abs()
        rows.append(sub.sort_values(["fdr_bh", "abs_delta"], ascending=[True, False]).head(1))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_summary(summary: pd.DataFrame) -> None:
    if summary.empty:
        return
    plot = summary.sort_values("ipf_minus_control_log1p_norm")
    colors = np.where(plot["ipf_minus_control_log1p_norm"] >= 0, "#B55D60", "#4C78A8")
    labels = plot["gene_symbol"] + " (" + plot["broad_celltype"] + ", " + plot["series_id"] + ")"
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    ax.barh(labels, plot["ipf_minus_control_log1p_norm"], color=colors)
    ax.axvline(0, color="#666666", lw=0.8)
    ax.set_xlabel("IPF-control donor-level pseudobulk log1p normalized expression")
    ax.set_title("Donor-aware pseudobulk validation of core candidates")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "single_cell_pseudobulk_core_candidates.png", dpi=300)
    fig.savefig(OUT_DIR / "single_cell_pseudobulk_core_candidates.pdf")
    plt.close(fig)


def main() -> None:
    frames = []
    qcs = []
    for loader in [load_gse135893, load_gse136831]:
        donor_df, qc = stream_pseudobulk(loader())
        frames.append(donor_df)
        qcs.append(qc)
    donor_df = pd.concat(frames, ignore_index=True)
    qc = pd.concat(qcs, ignore_index=True)
    stats = test_pseudobulk(donor_df)
    summary = summarize_core_contexts(stats)
    donor_df.to_csv(OUT_DIR / "single_cell_pseudobulk_donor_celltype_expression.csv", index=False, encoding="utf-8-sig")
    stats.to_csv(OUT_DIR / "single_cell_pseudobulk_differential_summary.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT_DIR / "single_cell_pseudobulk_core_candidate_summary.csv", index=False, encoding="utf-8-sig")
    qc.to_csv(OUT_DIR / "single_cell_pseudobulk_qc.csv", index=False, encoding="utf-8-sig")
    plot_summary(summary)
    print(OUT_DIR)
    print(qc.to_string(index=False))
    print(summary[["series_id", "broad_celltype", "gene_symbol", "ipf_minus_control_log1p_norm", "ipf_donors", "control_donors", "fdr_bh"]].to_string(index=False))


if __name__ == "__main__":
    main()
