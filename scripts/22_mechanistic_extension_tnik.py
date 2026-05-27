#!/usr/bin/env python
"""
Mechanistic extension for the IPF oligonucleotide manuscript.

This script adds a publication-oriented, reproducible layer inspired by
immune-module/WGCNA + single-cell communication + virtual perturbation papers,
but adapted to the data actually available in this IPF project.

Outputs are intentionally named as "marker-based", "WGCNA-like",
"curated ligand-receptor scoring", and "coexpression-neighborhood perturbation-priority proxy" unless the
exact external software package was used.

Triple QC:
  QC1 bulk marker/module inputs are aligned and non-empty
  QC2 single-cell ligand/receptor genes are recovered from sparse matrices
  QC3 TNIK is not over-claimed as a discovery/ML-selected target
"""

from __future__ import annotations

import gzip
import math
import os
import time
from collections import defaultdict
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "mechanistic_extension"
PLOT_DIR = OUT_DIR / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / "_mpl_config"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr


GEO_DIR = PROJECT_DIR / "data_raw" / "GEO"
FEATURE_DIR = PROJECT_DIR / "results" / "models" / "feature_matrices"

COL_IPF = "#B04A3A"
COL_CTRL = "#386FA4"
COL_ACCENT = "#4C7C59"
COL_NEUTRAL = "#6E7781"


MARKER_SETS = {
    "macrophage_myeloid": ["CD68", "LST1", "AIF1", "C1QA", "C1QB", "APOE", "SPP1", "MARCO", "MSR1", "TYROBP", "FCER1G"],
    "t_cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "IL7R", "GZMB"],
    "b_cell": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"],
    "fibroblast_matrix": ["COL1A1", "COL1A2", "COL3A1", "COL14A1", "DCN", "LUM", "POSTN", "ACTA2", "THY1", "ASPN"],
    "epithelial_ciliary": ["EPCAM", "KRT8", "KRT18", "FOXJ1", "TPPP3", "DNAI1", "TEKT1", "RSPH4A", "MUC5B", "SCGB1A1"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "KDR", "ENG", "EMCN"],
    "wnt_tnik_axis": ["TNIK", "CTNNB1", "TCF7L2", "LEF1", "WNT5A", "WNT7B", "FZD1", "FZD2", "FZD5", "LRP5", "LRP6"],
}


LR_PAIRS = [
    ("SPP1", "CD44", "SPP1-CD44", "macrophage-matrix immune remodeling"),
    ("SPP1", "ITGAV", "SPP1-ITGAV", "macrophage-integrin signaling"),
    ("SPP1", "ITGB1", "SPP1-ITGB1", "macrophage-integrin signaling"),
    ("SPP1", "ITGA4", "SPP1-ITGA4", "macrophage-integrin signaling"),
    ("COL1A1", "ITGA1", "COL1A1-ITGA1", "collagen-integrin matrix remodeling"),
    ("COL1A1", "ITGA2", "COL1A1-ITGA2", "collagen-integrin matrix remodeling"),
    ("COL1A1", "ITGB1", "COL1A1-ITGB1", "collagen-integrin matrix remodeling"),
    ("COL3A1", "ITGB1", "COL3A1-ITGB1", "collagen-integrin matrix remodeling"),
    ("COL14A1", "ITGB1", "COL14A1-ITGB1", "collagen-integrin matrix remodeling"),
    ("POSTN", "ITGAV", "POSTN-ITGAV", "periostin-integrin stromal signaling"),
    ("POSTN", "ITGB1", "POSTN-ITGB1", "periostin-integrin stromal signaling"),
    ("MIF", "CD74", "MIF-CD74", "MIF inflammatory signaling"),
    ("MIF", "CXCR4", "MIF-CXCR4", "MIF inflammatory signaling"),
    ("CXCL12", "CXCR4", "CXCL12-CXCR4", "chemokine recruitment"),
    ("CCL2", "CCR2", "CCL2-CCR2", "monocyte recruitment"),
    ("TGFB1", "TGFBR1", "TGFB1-TGFBR1", "fibrotic growth factor signaling"),
    ("TGFB1", "TGFBR2", "TGFB1-TGFBR2", "fibrotic growth factor signaling"),
    ("IL1B", "IL1R1", "IL1B-IL1R1", "inflammatory cytokine signaling"),
    ("TNF", "TNFRSF1A", "TNF-TNFRSF1A", "TNF inflammatory signaling"),
    ("TNFSF12", "TNFRSF12A", "TWEAK-FN14", "injury-remodeling signaling"),
    ("WNT5A", "FZD2", "WNT5A-FZD2", "Wnt/TNIK bridge"),
    ("WNT5A", "FZD5", "WNT5A-FZD5", "Wnt/TNIK bridge"),
    ("WNT7B", "FZD1", "WNT7B-FZD1", "Wnt/TNIK bridge"),
    ("WNT7B", "LRP6", "WNT7B-LRP6", "Wnt/TNIK bridge"),
]


PERTURB_TARGETS = ["SPP1", "POSTN", "COL14A1", "CD24", "PTGFRN", "CDH3", "TNIK"]


def normalize_gene(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().upper()
    return text if text and text != "NA" else ""


def bh_fdr(pvals: list[float]) -> np.ndarray:
    p = np.asarray([1.0 if pd.isna(x) else float(x) for x in pvals], dtype=float)
    order = np.argsort(p)
    ranked = np.empty_like(p)
    n = len(p)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        val = min(prev, p[idx] * n / (i + 1))
        ranked[idx] = val
        prev = val
    return ranked


def read_matrix(series_id: str) -> pd.DataFrame:
    path = FEATURE_DIR / f"{series_id}_gene_level_matrix.csv"
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    df["series_id"] = series_id
    return df


def zscore_frame(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mean = np.nanmean(arr, axis=0)
    sd = np.nanstd(arr, axis=0)
    sd[sd == 0] = np.nan
    out = (arr - mean) / sd
    out = np.nan_to_num(out, nan=0.0)
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def marker_scores() -> pd.DataFrame:
    rows = []
    for series in ["GSE32537", "GSE110147", "GSE150910", "GSE53845", "GSE92592"]:
        df = read_matrix(series)
        meta_cols = ["sample_id", "group", "label", "series_id"]
        genes = [c for c in df.columns if c not in meta_cols]
        expr = df[genes].apply(pd.to_numeric, errors="coerce")
        expr_z = zscore_frame(expr)
        out = df[["sample_id", "group", "label", "series_id"]].copy()
        for name, markers in MARKER_SETS.items():
            present = [g for g in markers if g in expr_z.columns]
            out[f"{name}_score"] = expr_z[present].mean(axis=1) if present else np.nan
            out[f"{name}_markers_found"] = len(present)
        rows.append(out)
    scores = pd.concat(rows, ignore_index=True)
    scores.to_csv(OUT_DIR / "bulk_marker_scores.csv", index=False, encoding="utf-8-sig")
    return scores


def build_modules(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = read_matrix("GSE32537")
    meta_cols = ["sample_id", "group", "label", "series_id"]
    genes = [c for c in df.columns if c not in meta_cols]
    expr = df[genes].apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[:, expr.notna().sum(axis=0) >= max(10, int(0.75 * len(expr)))]
    expr = expr.fillna(expr.median(numeric_only=True))
    var = expr.var(axis=0).sort_values(ascending=False)

    robust = pd.read_csv(PROJECT_DIR / "results/robust_candidates/robust_mrna_candidates_strict.csv")["gene_symbol"].map(normalize_gene).dropna().tolist()
    panel = pd.read_csv(PROJECT_DIR / "results/models/ml_outputs_discovery_only_mrna/ml_final_biomarker_panel.csv")["feature"].map(normalize_gene).dropna().tolist()
    priority = pd.read_csv(PROJECT_DIR / "results/submission_enhancements/final_target_priority_integrated.csv")["gene_symbol"].map(normalize_gene).dropna().head(80).tolist()
    lr_genes = sorted({g for pair in LR_PAIRS for g in pair[:2]})
    forced = set(robust + panel + priority + lr_genes + PERTURB_TARGETS + [g for vals in MARKER_SETS.values() for g in vals])
    forced = {g for g in forced if g in expr.columns}
    selected = list(dict.fromkeys(list(var.head(1400).index) + sorted(forced)))
    if len(selected) > 2200:
        top_set = set(var.head(1800).index) | forced
        selected = [g for g in selected if g in top_set]
    expr_sel = expr[selected]
    expr_z = zscore_frame(expr_sel)

    corr = np.corrcoef(expr_z.to_numpy().T)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    dist = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    module_ids = fcluster(z, t=10, criterion="maxclust")
    membership = pd.DataFrame({"gene_symbol": selected, "module": [f"M{m:02d}" for m in module_ids]})

    gse_scores = scores[scores["series_id"].eq("GSE32537")].set_index("sample_id")
    sample_order = df["sample_id"].astype(str).tolist()
    traits = pd.DataFrame(index=sample_order)
    traits["IPF_status"] = df["label"].astype(float).to_numpy()
    for col in [c for c in gse_scores.columns if c.endswith("_score")]:
        traits[col.replace("_score", "")] = gse_scores.reindex(sample_order)[col].to_numpy()

    module_trait_rows = []
    eigengenes = pd.DataFrame(index=sample_order)
    for module, mdf in membership.groupby("module"):
        mod_genes = mdf["gene_symbol"].tolist()
        mat = expr_z[mod_genes].to_numpy()
        mat = mat - mat.mean(axis=0, keepdims=True)
        u, s, vt = np.linalg.svd(mat, full_matrices=False)
        eig = u[:, 0] * s[0]
        if np.corrcoef(eig, traits["IPF_status"].to_numpy())[0, 1] < 0:
            eig = -eig
        eigengenes[module] = eig
        for trait in traits.columns:
            try:
                r, p = pearsonr(eig, traits[trait].astype(float).to_numpy())
            except Exception:
                r, p = np.nan, np.nan
            module_trait_rows.append(
                {
                    "module": module,
                    "module_size": len(mod_genes),
                    "trait": trait,
                    "pearson_r": r,
                    "p_value": p,
                }
            )
    module_traits = pd.DataFrame(module_trait_rows)
    module_traits["fdr"] = bh_fdr(module_traits["p_value"].tolist())

    gene_rows = []
    for module, mdf in membership.groupby("module"):
        eig = eigengenes[module].to_numpy()
        for gene in mdf["gene_symbol"]:
            try:
                kme = pearsonr(expr_z[gene].to_numpy(), eig)[0]
            except Exception:
                kme = np.nan
            gene_rows.append({"gene_symbol": gene, "module": module, "module_membership_kME": kme})
    membership = pd.DataFrame(gene_rows)
    membership["in_robust_mrna"] = membership["gene_symbol"].isin(set(robust))
    membership["in_ml_panel"] = membership["gene_symbol"].isin(set(panel))
    membership["in_priority_top80"] = membership["gene_symbol"].isin(set(priority))
    membership["is_tnik"] = membership["gene_symbol"].eq("TNIK")
    membership.to_csv(OUT_DIR / "coexpression_gene_module_membership.csv", index=False, encoding="utf-8-sig")
    module_traits.to_csv(OUT_DIR / "coexpression_module_trait_correlations.csv", index=False, encoding="utf-8-sig")

    overlap = (
        membership.groupby("module", as_index=False)
        .agg(
            module_size=("gene_symbol", "size"),
            robust_mrna_count=("in_robust_mrna", "sum"),
            ml_panel_count=("in_ml_panel", "sum"),
            priority_top80_count=("in_priority_top80", "sum"),
            contains_tnik=("is_tnik", "max"),
        )
        .merge(
            module_traits[module_traits["trait"].eq("IPF_status")][["module", "pearson_r", "p_value", "fdr"]].rename(
                columns={"pearson_r": "ipf_module_r", "p_value": "ipf_module_p", "fdr": "ipf_module_fdr"}
            ),
            on="module",
            how="left",
        )
    )
    overlap.to_csv(OUT_DIR / "candidate_module_overlap.csv", index=False, encoding="utf-8-sig")

    plot_module_heatmap(module_traits)
    plot_candidate_module_overlap(overlap)
    return membership, module_traits, overlap


def plot_module_heatmap(module_traits: pd.DataFrame) -> None:
    traits = ["IPF_status", "macrophage_myeloid", "fibroblast_matrix", "epithelial_ciliary", "wnt_tnik_axis", "endothelial", "t_cell", "b_cell"]
    pivot = module_traits[module_traits["trait"].isin(traits)].pivot(index="module", columns="trait", values="pearson_r").sort_index()
    fig, ax = plt.subplots(figsize=(9.5, max(4.5, 0.35 * len(pivot))))
    im = ax.imshow(pivot.fillna(0).to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            if pd.notna(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("WGCNA-like module-trait correlations in GSE32537")
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "module_trait_correlation_heatmap.png", dpi=220)
    fig.savefig(PLOT_DIR / "module_trait_correlation_heatmap.pdf")
    plt.close(fig)


def plot_candidate_module_overlap(overlap: pd.DataFrame) -> None:
    plot_df = overlap.sort_values("robust_mrna_count", ascending=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.barh(plot_df["module"], plot_df["robust_mrna_count"], color=COL_ACCENT, label="robust mRNA")
    ax.barh(plot_df["module"], plot_df["ml_panel_count"], color=COL_IPF, label="ML panel")
    ax.set_xlabel("Candidate genes in module")
    ax.set_title("Candidate concentration across coexpression modules")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "candidate_module_overlap.png", dpi=220)
    fig.savefig(PLOT_DIR / "candidate_module_overlap.pdf")
    plt.close(fig)


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


def stream_sc_expression(target_genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    datasets = [load_gse135893(target_genes), load_gse136831(target_genes)]
    all_rows = []
    qc_rows = []
    for ds in datasets:
        start = time.time()
        series_id = ds["series_id"]
        genes = ds["genes"]
        barcodes = ds["barcodes"]
        meta = ds["metadata"]
        meta = meta[meta["group"].isin(["IPF", "Control"])].drop_duplicates(subset=["cell_id"], keep="first").copy()
        meta_lookup = meta.set_index("cell_id")
        barcode_to_meta = []
        group_ct_counts = defaultdict(int)
        valid_cells = 0
        for barcode in barcodes:
            if barcode in meta_lookup.index:
                row = meta_lookup.loc[barcode]
                group = str(row["group"])
                broad = str(row["broad_celltype"])
                fine = str(row["fine_celltype"])
                lib = float(row["library_size"]) if pd.notna(row["library_size"]) and float(row["library_size"]) > 0 else math.nan
                barcode_to_meta.append((group, broad, fine, lib))
                valid_cells += 1
                group_ct_counts[(group, broad, fine)] += 1
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
        entries_seen = 0
        target_entries = 0
        negative_values = 0
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
                sums[key] += value / lib * 10000.0
                if value > 0:
                    detected[key] += 1
                target_entries += 1

        for (group, broad, fine), cell_count in group_ct_counts.items():
            for gene in sorted(target_found):
                key = (series_id, group, broad, fine, gene)
                mean_norm = sums.get(key, 0.0) / cell_count
                all_rows.append(
                    {
                        "series_id": series_id,
                        "group": group,
                        "broad_celltype": broad,
                        "fine_celltype": fine,
                        "gene_symbol": gene,
                        "cells": cell_count,
                        "mean_norm_count_per_10k": mean_norm,
                        "log1p_mean_norm": math.log1p(mean_norm),
                        "detected_cells": detected.get(key, 0),
                        "detection_fraction": detected.get(key, 0) / cell_count if cell_count else np.nan,
                    }
                )
        qc_rows.append(
            {
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
                "target_entries_streamed": target_entries,
                "negative_values": negative_values,
                "qc1_metadata_alignment_pass": n_cols == len(barcodes) and valid_cells > 0 and meta["group"].nunique() == 2,
                "qc2_matrix_integrity_pass": n_rows == len(genes) and entries_seen == nnz_header and negative_values == 0 and len(target_found) >= 10,
                "qc3_result_validity_pass": target_entries > 0,
                "runtime_seconds": round(time.time() - start, 2),
            }
        )
    summary = pd.DataFrame(all_rows)
    qc = pd.DataFrame(qc_rows)
    qc["triple_qc_pass"] = qc["qc1_metadata_alignment_pass"] & qc["qc2_matrix_integrity_pass"] & qc["qc3_result_validity_pass"]
    summary.to_csv(OUT_DIR / "single_cell_lr_gene_expression_summary.csv", index=False, encoding="utf-8-sig")
    qc.to_csv(OUT_DIR / "single_cell_lr_streaming_qc.csv", index=False, encoding="utf-8-sig")
    return summary, qc


def communication_scores(summary: pd.DataFrame) -> pd.DataFrame:
    excluded = ("multiplet", "outlier", "cellcycle", "cell cycle", "mt-trna", "mt-trnas", "nan")
    summary = summary[
        ~summary["broad_celltype"].astype(str).str.lower().apply(lambda x: any(token in x for token in excluded))
        & ~summary["fine_celltype"].astype(str).str.lower().apply(lambda x: any(token in x for token in excluded))
    ].copy()
    broad = (
        summary.groupby(["series_id", "group", "broad_celltype", "gene_symbol"], as_index=False)
        .agg(log1p_mean_norm=("log1p_mean_norm", "mean"), detection_fraction=("detection_fraction", "mean"), cells=("cells", "sum"))
    )
    lookup = broad.set_index(["series_id", "group", "broad_celltype", "gene_symbol"])
    celltypes = sorted(broad["broad_celltype"].dropna().unique())
    rows = []
    for series_id in sorted(broad["series_id"].unique()):
        for group in ["Control", "IPF"]:
            for sender in celltypes:
                for receiver in celltypes:
                    for ligand, receptor, pair_name, pathway in LR_PAIRS:
                        lig = lookup["log1p_mean_norm"].get((series_id, group, sender, ligand), np.nan)
                        rec = lookup["log1p_mean_norm"].get((series_id, group, receiver, receptor), np.nan)
                        lig_det = lookup["detection_fraction"].get((series_id, group, sender, ligand), np.nan)
                        rec_det = lookup["detection_fraction"].get((series_id, group, receiver, receptor), np.nan)
                        if pd.isna(lig) or pd.isna(rec):
                            continue
                        score = lig * rec
                        rows.append(
                            {
                                "series_id": series_id,
                                "group": group,
                                "sender_celltype": sender,
                                "receiver_celltype": receiver,
                                "ligand": ligand,
                                "receptor": receptor,
                                "pair": pair_name,
                                "pathway": pathway,
                                "ligand_log1p_mean": lig,
                                "receptor_log1p_mean": rec,
                                "ligand_detection_fraction": lig_det,
                                "receptor_detection_fraction": rec_det,
                                "interaction_score": score,
                            }
                        )
    scores = pd.DataFrame(rows)
    keys = ["series_id", "sender_celltype", "receiver_celltype", "ligand", "receptor", "pair", "pathway"]
    ipf = scores[scores["group"].eq("IPF")][keys + ["interaction_score"]].rename(columns={"interaction_score": "ipf_interaction_score"})
    ctrl = scores[scores["group"].eq("Control")][keys + ["interaction_score"]].rename(columns={"interaction_score": "control_interaction_score"})
    delta = ipf.merge(ctrl, on=keys, how="inner")
    delta["ipf_minus_control_interaction_score"] = delta["ipf_interaction_score"] - delta["control_interaction_score"]
    delta["interaction_log2_ratio"] = np.log2((delta["ipf_interaction_score"] + 0.01) / (delta["control_interaction_score"] + 0.01))
    delta = delta.sort_values("ipf_minus_control_interaction_score", ascending=False)
    delta.to_csv(OUT_DIR / "curated_ligand_receptor_interaction_delta.csv", index=False, encoding="utf-8-sig")
    top = delta.head(80)
    top.to_csv(OUT_DIR / "top_curated_ligand_receptor_interactions.csv", index=False, encoding="utf-8-sig")
    plot_lr(top)
    return delta


def plot_lr(top: pd.DataFrame) -> None:
    if top.empty:
        return
    plot_df = top.head(20).iloc[::-1].copy()
    plot_df["label"] = plot_df["pair"] + " | " + plot_df["sender_celltype"] + "->" + plot_df["receiver_celltype"]
    fig, ax = plt.subplots(figsize=(9, 7.2))
    colors = np.where(plot_df["pathway"].str.contains("Wnt", case=False, na=False), "#7B5EA7", COL_IPF)
    ax.barh(plot_df["label"], plot_df["ipf_minus_control_interaction_score"], color=colors)
    ax.set_xlabel("IPF - Control interaction score")
    ax.set_title("Top curated ligand-receptor score changes")
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "top_curated_ligand_receptor_interactions.png", dpi=220)
    fig.savefig(PLOT_DIR / "top_curated_ligand_receptor_interactions.pdf")
    plt.close(fig)


def tnik_evidence(membership: pd.DataFrame, module_traits: pd.DataFrame, sc_summary: pd.DataFrame, lr_delta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for path in sorted((PROJECT_DIR / "results/differential_expression_annotated").glob("*_gene_or_mirna_level.csv")):
        df = pd.read_csv(path)
        if "gene_symbol" not in df.columns:
            continue
        hit = df[df["gene_symbol"].astype(str).str.upper().eq("TNIK")].copy()
        if hit.empty:
            continue
        for _, row in hit.iterrows():
            rows.append(
                {
                    "evidence_layer": "bulk_differential_expression",
                    "series_id": row.get("series_id", ""),
                    "metric": "logFC/adj.P.Val",
                    "value": f"logFC={row.get('logFC', np.nan):.3g}; adj.P={row.get('adj.P.Val', np.nan):.3g}; role={row.get('dataset_role', '')}",
                    "supports_tnik_as_primary_discovery_target": bool(row.get("series_id", "") == "GSE32537" and float(row.get("adj.P.Val", 1.0)) < 0.05 and abs(float(row.get("logFC", 0.0))) >= 1),
                    "interpretation": "TNIK differential-expression evidence by cohort",
                }
            )
    tnik_module = membership.loc[membership["gene_symbol"].eq("TNIK"), "module"]
    if not tnik_module.empty:
        module = tnik_module.iloc[0]
        ipf_r = module_traits[(module_traits["module"].eq(module)) & (module_traits["trait"].eq("IPF_status"))]
        wnt_r = module_traits[(module_traits["module"].eq(module)) & (module_traits["trait"].eq("wnt_tnik_axis"))]
        rows.append(
            {
                "evidence_layer": "coexpression_module",
                "series_id": "GSE32537",
                "metric": "TNIK module",
                "value": f"{module}; IPF r={ipf_r['pearson_r'].iloc[0]:.3f}; Wnt/TNIK marker r={wnt_r['pearson_r'].iloc[0]:.3f}" if not ipf_r.empty and not wnt_r.empty else module,
                "supports_tnik_as_primary_discovery_target": False,
                "interpretation": "TNIK can be discussed as pathway-bridge evidence, not as a primary robust/ML-selected target",
            }
        )
    sc_tnik = sc_summary[sc_summary["gene_symbol"].eq("TNIK")].copy()
    if not sc_tnik.empty:
        top = sc_tnik.sort_values("log1p_mean_norm", ascending=False).head(5)
        rows.append(
            {
                "evidence_layer": "single_cell_expression",
                "series_id": "GSE135893/GSE136831",
                "metric": "top TNIK cell contexts",
                "value": "; ".join(f"{r.series_id}:{r.group}:{r.broad_celltype}={r.log1p_mean_norm:.2f}" for r in top.itertuples()),
                "supports_tnik_as_primary_discovery_target": False,
                "interpretation": "Cellular context available for wet validation design",
            }
        )
    wnt = lr_delta[lr_delta["pathway"].str.contains("Wnt", case=False, na=False)].head(10)
    if not wnt.empty:
        rows.append(
            {
                "evidence_layer": "curated_ligand_receptor_wnt_bridge",
                "series_id": "GSE135893/GSE136831",
                "metric": "top Wnt/TNIK upstream LR deltas",
                "value": "; ".join(f"{r.pair}:{r.sender_celltype}->{r.receiver_celltype} delta={r.ipf_minus_control_interaction_score:.2f}" for r in wnt.itertuples()),
                "supports_tnik_as_primary_discovery_target": False,
                "interpretation": "Wnt ligand-receptor changes can bridge an external TNIK oligonucleotide hypothesis",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "tnik_evidence_summary.csv", index=False, encoding="utf-8-sig")
    plot_tnik_de(out)
    return out


def plot_tnik_de(evidence: pd.DataFrame) -> None:
    bulk = evidence[evidence["evidence_layer"].eq("bulk_differential_expression")].copy()
    rows = []
    for item in bulk.itertuples():
        text = item.value
        try:
            logfc = float(text.split("logFC=")[1].split(";")[0])
            adjp = float(text.split("adj.P=")[1].split(";")[0])
        except Exception:
            continue
        rows.append({"series_id": item.series_id, "logFC": logfc, "neglog10_adj_p": -math.log10(max(adjp, 1e-300))})
    df = pd.DataFrame(rows)
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.axhline(-math.log10(0.05), color=COL_NEUTRAL, lw=0.8, ls="--")
    ax.bar(df["series_id"], df["logFC"], color=np.where(df["logFC"] > 0, COL_IPF, COL_CTRL))
    ax.set_ylabel("TNIK logFC (IPF vs Control)")
    ax.set_title("TNIK is validation-supported but not discovery-selected")
    for i, row in df.iterrows():
        ax.text(i, row["logFC"], f"FDR={10 ** (-row['neglog10_adj_p']):.2g}", ha="center", va="bottom" if row["logFC"] >= 0 else "top", fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "tnik_bulk_evidence.png", dpi=220)
    fig.savefig(PLOT_DIR / "tnik_bulk_evidence.pdf")
    plt.close(fig)


def perturbation_priority_proxy(membership: pd.DataFrame, overlap: pd.DataFrame) -> pd.DataFrame:
    df = read_matrix("GSE32537")
    meta_cols = ["sample_id", "group", "label", "series_id"]
    genes = [c for c in df.columns if c not in meta_cols]
    expr = df[genes].apply(pd.to_numeric, errors="coerce").fillna(df[genes].median(numeric_only=True))
    expr_z = zscore_frame(expr)
    de = pd.read_csv(PROJECT_DIR / "results/differential_expression_annotated/GSE32537_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv")
    de = de.drop_duplicates("gene_symbol").set_index("gene_symbol")
    module_lookup = membership.set_index("gene_symbol")["module"].to_dict()
    ipf_module_r = overlap.set_index("module")["ipf_module_r"].to_dict()
    rows = []
    for target in PERTURB_TARGETS:
        if target not in expr_z.columns:
            continue
        target_vec = expr_z[target].to_numpy()
        corrs = {}
        for gene in expr_z.columns:
            if gene == target:
                continue
            try:
                r = pearsonr(target_vec, expr_z[gene].to_numpy())[0]
            except Exception:
                r = np.nan
            if pd.notna(r):
                corrs[gene] = r
        corr_s = pd.Series(corrs).sort_values(key=lambda x: x.abs(), ascending=False)
        neighborhood = corr_s[corr_s.abs() >= 0.35].head(200)
        target_logfc = float(de.loc[target, "logFC"]) if target in de.index else np.nan
        target_adjp = float(de.loc[target, "adj.P.Val"]) if target in de.index else np.nan
        if pd.isna(target_logfc):
            perturbation_direction = "unknown"
            reversal_score = np.nan
        elif target_logfc > 0:
            perturbation_direction = "knockdown aligns with IPF reversal"
            reversal_score = float((neighborhood.abs() * np.sign(neighborhood) * np.sign(target_logfc)).sum())
        else:
            perturbation_direction = "knockdown may not align; restoration/inhibition-context test needed"
            reversal_score = float((neighborhood.abs() * np.sign(neighborhood) * np.sign(target_logfc)).sum())
        mod = module_lookup.get(target, "")
        rows.append(
            {
                "target_gene": target,
                "target_module": mod,
                "target_module_ipf_r": ipf_module_r.get(mod, np.nan),
                "gse32537_logFC": target_logfc,
                "gse32537_adj_p": target_adjp,
                "coexpression_neighbors_abs_r_ge_0_35": int(len(neighborhood)),
                "top_positive_neighbors": ";".join(neighborhood[neighborhood > 0].head(8).index.tolist()),
                "top_negative_neighbors": ";".join(neighborhood[neighborhood < 0].head(8).index.tolist()),
                "perturbation_priority_proxy_score": reversal_score,
                "oligonucleotide_interpretation": perturbation_direction,
                "claim_strength": "primary_candidate" if target != "TNIK" and pd.notna(target_logfc) and target_logfc > 0 else "bridge_or_context_dependent",
            }
        )
    out = pd.DataFrame(rows).sort_values("perturbation_priority_proxy_score", ascending=False)
    out.to_csv(OUT_DIR / "perturbation_priority_proxy_summary.csv", index=False, encoding="utf-8-sig")
    plot_perturbation_priority(out)
    return out


def plot_perturbation_priority(out: pd.DataFrame) -> None:
    if out.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    colors = np.where(out["target_gene"].eq("TNIK"), "#7B5EA7", np.where(out["gse32537_logFC"] > 0, COL_IPF, COL_CTRL))
    ax.bar(out["target_gene"], out["perturbation_priority_proxy_score"], color=colors)
    ax.axhline(0, color=COL_NEUTRAL, lw=0.8)
    ax.set_ylabel("Network-neighborhood proxy score")
    ax.set_title("Coexpression-neighborhood proxy prioritizes candidates")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "perturbation_priority_proxy.png", dpi=220)
    fig.savefig(PLOT_DIR / "perturbation_priority_proxy.pdf")
    plt.close(fig)


def export_excel() -> None:
    with pd.ExcelWriter(OUT_DIR / "mechanistic_extension_summary.xlsx", engine="openpyxl") as writer:
        for name in [
            "bulk_marker_scores",
            "coexpression_module_trait_correlations",
            "coexpression_gene_module_membership",
            "candidate_module_overlap",
            "single_cell_lr_streaming_qc",
            "single_cell_lr_gene_expression_summary",
            "curated_ligand_receptor_interaction_delta",
            "top_curated_ligand_receptor_interactions",
            "perturbation_priority_proxy_summary",
            "tnik_evidence_summary",
            "mechanistic_extension_qc",
        ]:
            path = OUT_DIR / f"{name}.csv"
            if path.exists():
                pd.read_csv(path).to_excel(writer, sheet_name=name[:31], index=False)


def main() -> None:
    scores = marker_scores()
    membership, module_traits, overlap = build_modules(scores)
    lr_genes = sorted({g for pair in LR_PAIRS for g in pair[:2]} | set(PERTURB_TARGETS) | {"TNIK", "CTNNB1", "TCF7L2", "LEF1"})
    sc_summary_path = OUT_DIR / "single_cell_lr_gene_expression_summary.csv"
    sc_qc_path = OUT_DIR / "single_cell_lr_streaming_qc.csv"
    if sc_summary_path.exists() and sc_qc_path.exists():
        sc_summary = pd.read_csv(sc_summary_path)
        sc_qc = pd.read_csv(sc_qc_path)
    else:
        sc_summary, sc_qc = stream_sc_expression(lr_genes)
    lr_delta = communication_scores(sc_summary)
    perturb = perturbation_priority_proxy(membership, overlap)
    tnik = tnik_evidence(membership, module_traits, sc_summary, lr_delta)

    bulk_discovery_ok = scores[scores["series_id"].eq("GSE32537")].shape[0] > 20 and module_traits["module"].nunique() >= 5
    sc_ok = bool(sc_qc["triple_qc_pass"].all()) and not lr_delta.empty
    tnik_not_overclaimed = not bool(tnik["supports_tnik_as_primary_discovery_target"].fillna(False).any())
    qc = pd.DataFrame(
        [
            {
                "bulk_samples_scored": int(len(scores)),
                "gse32537_modules": int(module_traits["module"].nunique()),
                "lr_interactions_scored": int(len(lr_delta)),
                "perturbation_priority_targets_scored": int(len(perturb)),
                "tnik_evidence_rows": int(len(tnik)),
                "qc1_bulk_module_pass": bulk_discovery_ok,
                "qc2_single_cell_lr_pass": sc_ok,
                "qc3_tnik_not_overclaimed_pass": tnik_not_overclaimed,
            }
        ]
    )
    qc["triple_qc_pass"] = qc["qc1_bulk_module_pass"] & qc["qc2_single_cell_lr_pass"] & qc["qc3_tnik_not_overclaimed_pass"]
    qc.to_csv(OUT_DIR / "mechanistic_extension_qc.csv", index=False, encoding="utf-8-sig")
    export_excel()
    print("Mechanistic extension outputs written to:")
    print(OUT_DIR)
    print(qc.to_string(index=False))
    print("\nTop module overlaps:")
    print(overlap.sort_values(["robust_mrna_count", "ml_panel_count"], ascending=False).head(8).to_string(index=False))
    print("\nTop LR changes:")
    print(lr_delta.head(10)[["series_id", "pair", "sender_celltype", "receiver_celltype", "ipf_minus_control_interaction_score"]].to_string(index=False))
    print("\nTNIK evidence:")
    print(tnik.to_string(index=False))


if __name__ == "__main__":
    main()
