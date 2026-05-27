#!/usr/bin/env python
"""
Build a broad-celltype mapping table for single-cell localization outputs.

The pseudobulk validation uses broad cell types. This table records how fine
cell labels appearing in the expression summaries map into those broad groups.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
IN = PROJECT_DIR / "results" / "single_cell_validation" / "single_cell_target_expression_summary_all.csv"
OUT = PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_broad_celltype_mapping.csv"
QC_OUT = PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_broad_celltype_mapping_qc.csv"


def main() -> None:
    df = pd.read_csv(IN)
    mapping = (
        df.groupby(["series_id", "broad_celltype", "fine_celltype"], dropna=False)
        .agg(
            groups_observed=("group", lambda x: ";".join(sorted(set(map(str, x))))),
            max_cells_in_group=("cells", "max"),
            total_rows=("gene_symbol", "size"),
            unique_genes=("gene_symbol", "nunique"),
        )
        .reset_index()
        .sort_values(["series_id", "broad_celltype", "fine_celltype"])
    )
    mapping["aggregation_use"] = "used for donor/sample x broad-celltype pseudobulk summaries when donor metadata and matrix coverage were available"
    mapping.to_csv(OUT, index=False, encoding="utf-8-sig")

    qc = (
        mapping.groupby(["series_id", "broad_celltype"], dropna=False)
        .agg(fine_celltypes=("fine_celltype", "nunique"), rows=("fine_celltype", "size"))
        .reset_index()
    )
    qc.to_csv(QC_OUT, index=False, encoding="utf-8-sig")
    print(OUT)
    print(QC_OUT)


if __name__ == "__main__":
    main()
