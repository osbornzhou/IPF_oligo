#!/usr/bin/env python
"""Create manuscript-ready single-cell validation summary tables and plots."""

from pathlib import Path
import os

import matplotlib

_mpl_config_dir = Path.cwd() / "results" / "single_cell_validation" / "matplotlib_cache"
_mpl_config_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_config_dir))
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SC_DIR = PROJECT_DIR / "results" / "single_cell_validation"
PLOT_DIR = SC_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def is_clean_celltype(row: pd.Series) -> bool:
    broad = str(row["broad_celltype"])
    fine = str(row["fine_celltype"])
    bad_terms = ["Multiplet", "Outlier", "MT-tRNAs", "CellCycle"]
    return not any(term in broad or term in fine for term in bad_terms)


def main() -> None:
    summary = pd.read_csv(SC_DIR / "single_cell_target_expression_summary_all.csv")
    delta = pd.read_csv(SC_DIR / "single_cell_ipf_control_delta_by_celltype.csv")
    targets = pd.read_csv(SC_DIR / "single_cell_target_genes.csv")
    qc = pd.read_csv(SC_DIR / "single_cell_validation_triple_qc.csv")

    clean_delta = delta[delta.apply(is_clean_celltype, axis=1)].copy()
    clean_summary = summary[summary.apply(is_clean_celltype, axis=1)].copy()

    clean_delta = clean_delta[
        (clean_delta["ipf_cells"] >= 20)
        & (clean_delta["control_cells"] >= 20)
    ].copy()

    clean_delta["abs_delta"] = clean_delta["ipf_minus_control_log1p_mean_norm"].abs()
    clean_delta = clean_delta.sort_values(["abs_delta", "source_count"], ascending=[False, False])
    clean_delta.to_csv(SC_DIR / "single_cell_clean_ipf_control_delta_by_celltype.csv", index=False, encoding="utf-8-sig")
    clean_delta.head(200).to_csv(SC_DIR / "single_cell_clean_top_celltype_gene_changes.csv", index=False, encoding="utf-8-sig")

    localization = (
        clean_summary.groupby(["gene_symbol", "series_id", "broad_celltype"], as_index=False)
        .agg(
            mean_log1p_expression=("log1p_mean_norm", "mean"),
            mean_detection_fraction=("detection_fraction", "mean"),
            total_cells=("cells", "sum"),
        )
        .sort_values(["gene_symbol", "mean_log1p_expression"], ascending=[True, False])
    )
    top_localization = localization.sort_values(["gene_symbol", "mean_log1p_expression"], ascending=[True, False]).groupby("gene_symbol").head(3)
    top_localization.to_csv(SC_DIR / "single_cell_gene_top_celltype_localization.csv", index=False, encoding="utf-8-sig")

    manuscript_qc = pd.DataFrame(
        [
            {
                "available_expression_datasets": int(qc["triple_qc_pass"].sum()),
                "expression_datasets_triple_qc_pass": ";".join(qc.loc[qc["triple_qc_pass"] == True, "series_id"].astype(str)),
                "excluded_dataset": "GSE122960",
                "excluded_dataset_reason": "cell-level expression matrix not available locally",
                "target_genes": targets.shape[0],
                "clean_celltype_delta_rows": clean_delta.shape[0],
                "clean_celltype_delta_requires_min_ipf_control_cells": 20,
                "cleaning_rule": "exclude Multiplet, Outlier, MT-tRNAs, and CellCycle-labeled fine cell types",
                "manuscript_summary_qc_pass": clean_delta.shape[0] > 0 and int(qc["triple_qc_pass"].sum()) >= 2,
            }
        ]
    )
    manuscript_qc.to_csv(SC_DIR / "single_cell_manuscript_summary_qc.csv", index=False, encoding="utf-8-sig")

    top_genes = targets.sort_values(["source_count", "best_source_rank"], ascending=[False, True]).head(25)["gene_symbol"].tolist()
    for series_id, sdf in clean_delta[clean_delta["gene_symbol"].isin(top_genes)].groupby("series_id"):
        broad = (
            sdf.groupby(["gene_symbol", "broad_celltype"], as_index=False)["ipf_minus_control_log1p_mean_norm"]
            .mean()
        )
        if broad.empty:
            continue
        matrix = broad.pivot(index="gene_symbol", columns="broad_celltype", values="ipf_minus_control_log1p_mean_norm")
        matrix = matrix.reindex([g for g in top_genes if g in matrix.index])
        vmax = np.nanmax(np.abs(matrix.to_numpy()))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1
        fig, ax = plt.subplots(figsize=(max(7, matrix.shape[1] * 1.1), max(6, matrix.shape[0] * 0.3)))
        im = ax.imshow(matrix.fillna(0).to_numpy(), cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(matrix.shape[1]))
        ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(matrix.shape[0]))
        ax.set_yticklabels(matrix.index, fontsize=8)
        ax.set_title(f"{series_id} clean cell types: IPF-Control delta")
        fig.colorbar(im, ax=ax, label="IPF - Control log1p mean normalized expression")
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"{series_id}_clean_single_cell_delta_heatmap.png", dpi=180)
        plt.close(fig)

    with pd.ExcelWriter(SC_DIR / "single_cell_validation_manuscript_summary.xlsx", engine="openpyxl") as writer:
        manuscript_qc.to_excel(writer, sheet_name="manuscript_qc", index=False)
        clean_delta.head(500).to_excel(writer, sheet_name="clean_top_changes", index=False)
        top_localization.to_excel(writer, sheet_name="gene_localization", index=False)
        targets.to_excel(writer, sheet_name="target_genes", index=False)

    print(manuscript_qc.to_string(index=False))
    print("\nClean top changes:")
    print(clean_delta[["series_id", "gene_symbol", "broad_celltype", "fine_celltype", "ipf_minus_control_log1p_mean_norm", "ipf_detection_fraction", "control_detection_fraction"]].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
