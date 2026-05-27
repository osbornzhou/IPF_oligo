#!/usr/bin/env python
"""
Single-cell validation of ML panel, PPI hubs, and miRNA-mRNA axis targets.

This script streams raw sparse matrices and extracts only the target genes.
It computes cell-type-level normalized expression summaries without loading
the full sparse matrices into memory.

Triple QC:
  QC1 input/metadata alignment
  QC2 sparse matrix integrity and target gene recovery
  QC3 cell-type-level result validity
"""

from __future__ import annotations

import gzip
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import matplotlib

_mpl_config_dir = Path.cwd() / "results" / "single_cell_validation" / "matplotlib_cache"
_mpl_config_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_config_dir))
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "single_cell_validation"
PLOT_DIR = OUT_DIR / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

GEO_DIR = PROJECT_DIR / "data_raw" / "GEO"


def normalize_gene(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.upper()
    return text if text and text != "NA" else ""


def read_target_genes() -> pd.DataFrame:
    records = []

    panel = pd.read_csv(PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna" / "ml_final_biomarker_panel.csv")
    for _, row in panel.iterrows():
        records.append({"gene_symbol": normalize_gene(row["feature"]), "source": "ml_final_panel", "source_rank": int(row.name) + 1})

    priority = pd.read_csv(PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna" / "ml_oligonucleotide_target_priority.csv")
    for _, row in priority.head(50).iterrows():
        records.append({"gene_symbol": normalize_gene(row["feature"]), "source": "ml_target_priority_top50", "source_rank": int(row.name) + 1})

    ppi = pd.read_csv(PROJECT_DIR / "results" / "ppi_network" / "string_ppi_hub_genes_robust_mrna_strict_medium_confidence.csv")
    for _, row in ppi.head(30).iterrows():
        records.append({"gene_symbol": normalize_gene(row["gene_symbol"]), "source": "ppi_hub_top30", "source_rank": int(row.name) + 1})

    axes = pd.read_csv(PROJECT_DIR / "results" / "mirna_mrna_axes" / "robust_mirna_mrna_negative_axes_mirtarbase.csv")
    for _, row in axes.iterrows():
        records.append({"gene_symbol": normalize_gene(row["target_gene"]), "source": "mirna_mrna_axis_target", "source_rank": int(row.name) + 1})

    target = pd.DataFrame(records)
    target = target[target["gene_symbol"] != ""].copy()
    collapsed = (
        target.groupby("gene_symbol", as_index=False)
        .agg(
            sources=("source", lambda x: ";".join(sorted(set(x)))),
            source_count=("source", "nunique"),
            best_source_rank=("source_rank", "min"),
        )
        .sort_values(["source_count", "best_source_rank", "gene_symbol"], ascending=[False, True, True])
    )
    collapsed.to_csv(OUT_DIR / "single_cell_target_genes.csv", index=False, encoding="utf-8-sig")
    return collapsed


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", errors="replace") as handle:
        return [line.rstrip("\n").strip('"') for line in handle]


def load_gse135893(target_genes: set[str]) -> dict:
    supp = GEO_DIR / "GSE135893" / "supplementary"
    genes = [normalize_gene(x) for x in read_gzip_lines(supp / "GSE135893_genes.tsv.gz")]
    barcodes = read_gzip_lines(supp / "GSE135893_barcodes.tsv.gz")
    meta = pd.read_csv(supp / "GSE135893_IPF_metadata.csv.gz")
    first_col = meta.columns[0]
    meta = meta.rename(columns={first_col: "cell_id"})
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["group"] = meta["Diagnosis"].astype(str).replace({"control": "Control", "IPF": "IPF"})
    meta["fine_celltype"] = meta["celltype"].astype(str)
    meta["broad_celltype"] = meta["population"].astype(str)
    meta["library_size"] = pd.to_numeric(meta["nCount_RNA"], errors="coerce")
    return {
        "series_id": "GSE135893",
        "matrix_path": supp / "GSE135893_matrix.mtx.gz",
        "genes": genes,
        "barcodes": barcodes,
        "metadata": meta,
        "target_genes": target_genes,
    }


def load_gse136831(target_genes: set[str]) -> dict:
    supp = GEO_DIR / "GSE136831" / "supplementary"
    gene_df = pd.read_csv(supp / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    genes = [normalize_gene(x) for x in gene_df["HGNC_EnsemblAlt_GeneID"]]
    barcodes = read_gzip_lines(supp / "GSE136831_AllCells.cellBarcodes.txt.gz")
    meta = pd.read_csv(supp / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    meta = meta.rename(columns={"CellBarcode_Identity": "cell_id"})
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["group"] = meta["Disease_Identity"].astype(str)
    meta["fine_celltype"] = meta["Subclass_Cell_Identity"].astype(str)
    meta["broad_celltype"] = meta["CellType_Category"].astype(str)
    meta["library_size"] = pd.to_numeric(meta["nUMI"], errors="coerce")
    return {
        "series_id": "GSE136831",
        "matrix_path": supp / "GSE136831_RawCounts_Sparse.mtx.gz",
        "genes": genes,
        "barcodes": barcodes,
        "metadata": meta,
        "target_genes": target_genes,
    }


def parse_mtx_header(handle) -> tuple[int, int, int]:
    for line in handle:
        text = line.strip()
        if not text or text.startswith("%"):
            continue
        parts = text.split()
        if len(parts) >= 3:
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise RuntimeError("Matrix Market header not found")


def stream_dataset(ds: dict) -> tuple[pd.DataFrame, dict]:
    series_id = ds["series_id"]
    genes = ds["genes"]
    barcodes = ds["barcodes"]
    meta = ds["metadata"]
    target_genes = ds["target_genes"]

    meta = meta[meta["group"].isin(["IPF", "Control"])].copy()
    meta = meta.drop_duplicates(subset=["cell_id"], keep="first")
    meta_lookup = meta.set_index("cell_id")

    barcode_to_meta = []
    valid_cells = 0
    group_celltype_counts = defaultdict(int)
    for barcode in barcodes:
        if barcode in meta_lookup.index:
            row = meta_lookup.loc[barcode]
            group = str(row["group"])
            broad = str(row["broad_celltype"])
            fine = str(row["fine_celltype"])
            lib = float(row["library_size"]) if pd.notna(row["library_size"]) and float(row["library_size"]) > 0 else math.nan
            barcode_to_meta.append((group, broad, fine, lib))
            valid_cells += 1
            group_celltype_counts[(group, broad, fine)] += 1
        else:
            barcode_to_meta.append(None)

    row_to_gene = {}
    target_found = set()
    for idx, gene in enumerate(genes, start=1):
        if gene in target_genes:
            row_to_gene[idx] = gene
            target_found.add(gene)

    sums = defaultdict(float)
    detected = defaultdict(int)
    raw_sums = defaultdict(float)
    entries_seen = 0
    target_entries = 0
    negative_values = 0
    start = time.time()

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
            if value < 0:
                negative_values += 1
            gene = row_to_gene.get(row_idx)
            if gene is None:
                continue
            cell_meta = barcode_to_meta[col_idx - 1] if 0 <= col_idx - 1 < len(barcode_to_meta) else None
            if cell_meta is None:
                continue
            group, broad, fine, lib = cell_meta
            if not math.isfinite(lib) or lib <= 0:
                continue
            key = (series_id, group, broad, fine, gene)
            norm = value / lib * 10000.0
            sums[key] += norm
            raw_sums[key] += value
            if value > 0:
                detected[key] += 1
            target_entries += 1

    rows = []
    for (sid, group, broad, fine, gene), total_cells in [
        ((series_id, group, broad, fine, gene), count)
        for (group, broad, fine), count in group_celltype_counts.items()
        for gene in sorted(target_found)
    ]:
        key = (sid, group, broad, fine, gene)
        mean_norm = sums.get(key, 0.0) / total_cells
        rows.append(
            {
                "series_id": sid,
                "group": group,
                "broad_celltype": broad,
                "fine_celltype": fine,
                "gene_symbol": gene,
                "cells": total_cells,
                "mean_norm_count_per_10k": mean_norm,
                "log1p_mean_norm": math.log1p(mean_norm),
                "raw_count_sum": raw_sums.get(key, 0.0),
                "detected_cells": detected.get(key, 0),
                "detection_fraction": detected.get(key, 0) / total_cells if total_cells else math.nan,
            }
        )

    qc = {
        "series_id": series_id,
        "genes_in_reference": len(genes),
        "barcodes_in_reference": len(barcodes),
        "metadata_cells_ipf_control": len(meta),
        "barcodes_with_metadata": valid_cells,
        "mtx_rows": n_rows,
        "mtx_cols": n_cols,
        "mtx_nnz_header": nnz_header,
        "mtx_entries_streamed": entries_seen,
        "target_genes_requested": len(target_genes),
        "target_genes_found": len(target_found),
        "target_genes_missing": len(target_genes - target_found),
        "target_entries_streamed": target_entries,
        "negative_values": negative_values,
        "group_count": meta["group"].nunique(),
        "broad_celltype_count": meta["broad_celltype"].nunique(),
        "fine_celltype_count": meta["fine_celltype"].nunique(),
        "runtime_seconds": round(time.time() - start, 2),
        "qc1_metadata_alignment_pass": n_cols == len(barcodes) and valid_cells > 0 and meta["group"].nunique() == 2,
        "qc2_matrix_integrity_pass": n_rows == len(genes) and entries_seen == nnz_header and negative_values == 0 and len(target_found) >= 5,
        "qc3_result_validity_pass": len(rows) > 0 and target_entries > 0,
    }
    qc["triple_qc_pass"] = qc["qc1_metadata_alignment_pass"] and qc["qc2_matrix_integrity_pass"] and qc["qc3_result_validity_pass"]
    return pd.DataFrame(rows), qc


def compute_ipf_control_delta(summary: pd.DataFrame) -> pd.DataFrame:
    keys = ["series_id", "broad_celltype", "fine_celltype", "gene_symbol"]
    ipf = summary[summary["group"] == "IPF"][keys + ["log1p_mean_norm", "detection_fraction", "cells"]].rename(
        columns={
            "log1p_mean_norm": "ipf_log1p_mean_norm",
            "detection_fraction": "ipf_detection_fraction",
            "cells": "ipf_cells",
        }
    )
    control = summary[summary["group"] == "Control"][keys + ["log1p_mean_norm", "detection_fraction", "cells"]].rename(
        columns={
            "log1p_mean_norm": "control_log1p_mean_norm",
            "detection_fraction": "control_detection_fraction",
            "cells": "control_cells",
        }
    )
    delta = ipf.merge(control, on=keys, how="inner")
    delta["ipf_minus_control_log1p_mean_norm"] = delta["ipf_log1p_mean_norm"] - delta["control_log1p_mean_norm"]
    delta["ipf_minus_control_detection_fraction"] = delta["ipf_detection_fraction"] - delta["control_detection_fraction"]
    return delta.sort_values(["gene_symbol", "series_id", "broad_celltype", "fine_celltype"])


def plot_heatmaps(delta: pd.DataFrame, target_genes: pd.DataFrame) -> None:
    top_genes = target_genes.sort_values(["source_count", "best_source_rank"], ascending=[False, True]).head(30)["gene_symbol"].tolist()
    for series_id, sdf in delta.groupby("series_id"):
        broad = (
            sdf[sdf["gene_symbol"].isin(top_genes)]
            .groupby(["gene_symbol", "broad_celltype"], as_index=False)["ipf_minus_control_log1p_mean_norm"]
            .mean()
        )
        if broad.empty:
            continue
        matrix = broad.pivot(index="gene_symbol", columns="broad_celltype", values="ipf_minus_control_log1p_mean_norm").reindex(top_genes)
        matrix = matrix.dropna(how="all")
        fig, ax = plt.subplots(figsize=(max(7, matrix.shape[1] * 1.1), max(6, matrix.shape[0] * 0.28)))
        im = ax.imshow(matrix.fillna(0).to_numpy(), cmap="RdBu_r", aspect="auto", vmin=-np.nanmax(np.abs(matrix.to_numpy())), vmax=np.nanmax(np.abs(matrix.to_numpy())))
        ax.set_xticks(range(matrix.shape[1]))
        ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(matrix.shape[0]))
        ax.set_yticklabels(matrix.index, fontsize=8)
        ax.set_title(f"{series_id} IPF-Control expression delta")
        fig.colorbar(im, ax=ax, label="IPF - Control log1p mean normalized expression")
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"{series_id}_single_cell_ipf_control_delta_heatmap.png", dpi=180)
        plt.close(fig)


def plot_dotplot(summary: pd.DataFrame, target_genes: pd.DataFrame) -> None:
    top_genes = target_genes.sort_values(["source_count", "best_source_rank"], ascending=[False, True]).head(25)["gene_symbol"].tolist()
    for series_id, sdf in summary.groupby("series_id"):
        sdf = sdf[(sdf["gene_symbol"].isin(top_genes)) & (sdf["group"].isin(["IPF", "Control"]))].copy()
        broad = (
            sdf.groupby(["group", "broad_celltype", "gene_symbol"], as_index=False)
            .agg(log1p_mean_norm=("log1p_mean_norm", "mean"), detection_fraction=("detection_fraction", "mean"))
        )
        broad["x"] = broad["broad_celltype"] + " | " + broad["group"]
        x_order = sorted(broad["x"].unique())
        gene_order = [g for g in top_genes if g in set(broad["gene_symbol"])]
        broad["x_idx"] = broad["x"].map({x: i for i, x in enumerate(x_order)})
        broad["y_idx"] = broad["gene_symbol"].map({g: i for i, g in enumerate(gene_order)})
        fig, ax = plt.subplots(figsize=(max(9, len(x_order) * 0.55), max(6, len(gene_order) * 0.32)))
        sc = ax.scatter(
            broad["x_idx"],
            broad["y_idx"],
            s=20 + 180 * broad["detection_fraction"],
            c=broad["log1p_mean_norm"],
            cmap="viridis",
            alpha=0.85,
        )
        ax.set_xticks(range(len(x_order)))
        ax.set_xticklabels(x_order, rotation=70, ha="right", fontsize=7)
        ax.set_yticks(range(len(gene_order)))
        ax.set_yticklabels(gene_order, fontsize=8)
        ax.set_title(f"{series_id} target gene expression by broad cell type")
        fig.colorbar(sc, ax=ax, label="log1p mean normalized expression")
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"{series_id}_single_cell_target_gene_dotplot.png", dpi=180)
        plt.close(fig)


def main() -> None:
    targets = read_target_genes()
    target_set = set(targets["gene_symbol"])

    datasets = [load_gse135893(target_set), load_gse136831(target_set)]
    summaries = []
    qc_rows = []
    for ds in datasets:
        print(f"Streaming {ds['series_id']} sparse matrix for {len(target_set)} target genes...")
        summary, qc = stream_dataset(ds)
        summaries.append(summary)
        qc_rows.append(qc)
        summary.to_csv(OUT_DIR / f"{ds['series_id']}_single_cell_target_expression_summary.csv", index=False, encoding="utf-8-sig")
        print(qc)

    unavailable = {
        "series_id": "GSE122960",
        "genes_in_reference": 0,
        "barcodes_in_reference": 0,
        "metadata_cells_ipf_control": 0,
        "barcodes_with_metadata": 0,
        "mtx_rows": 0,
        "mtx_cols": 0,
        "mtx_nnz_header": 0,
        "mtx_entries_streamed": 0,
        "target_genes_requested": len(target_set),
        "target_genes_found": 0,
        "target_genes_missing": len(target_set),
        "target_entries_streamed": 0,
        "negative_values": 0,
        "group_count": 0,
        "broad_celltype_count": 0,
        "fine_celltype_count": 0,
        "runtime_seconds": 0,
        "qc1_metadata_alignment_pass": False,
        "qc2_matrix_integrity_pass": False,
        "qc3_result_validity_pass": False,
        "triple_qc_pass": False,
        "note": "Cell-level expression matrix not available locally; excluded from single-cell expression validation.",
    }
    qc_rows.append(unavailable)

    all_summary = pd.concat(summaries, ignore_index=True)
    all_summary = all_summary.merge(targets, on="gene_symbol", how="left")
    all_summary.to_csv(OUT_DIR / "single_cell_target_expression_summary_all.csv", index=False, encoding="utf-8-sig")

    delta = compute_ipf_control_delta(all_summary)
    delta = delta.merge(targets, on="gene_symbol", how="left")
    delta.to_csv(OUT_DIR / "single_cell_ipf_control_delta_by_celltype.csv", index=False, encoding="utf-8-sig")

    top_delta = delta.reindex(delta["ipf_minus_control_log1p_mean_norm"].abs().sort_values(ascending=False).index).head(200)
    top_delta.to_csv(OUT_DIR / "single_cell_top_ipf_control_celltype_gene_changes.csv", index=False, encoding="utf-8-sig")

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT_DIR / "single_cell_validation_triple_qc.csv", index=False, encoding="utf-8-sig")

    plot_heatmaps(delta, targets)
    plot_dotplot(all_summary, targets)

    with pd.ExcelWriter(OUT_DIR / "single_cell_validation_summary.xlsx", engine="openpyxl") as writer:
        qc.to_excel(writer, sheet_name="triple_qc", index=False)
        targets.to_excel(writer, sheet_name="target_genes", index=False)
        top_delta.to_excel(writer, sheet_name="top_celltype_changes", index=False)
        delta.to_excel(writer, sheet_name="ipf_control_delta", index=False)

    print("\nQC:")
    print(qc.to_string(index=False))
    print("\nTop cell-type changes:")
    print(top_delta[["series_id", "gene_symbol", "broad_celltype", "fine_celltype", "ipf_minus_control_log1p_mean_norm", "ipf_detection_fraction", "control_detection_fraction"]].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
