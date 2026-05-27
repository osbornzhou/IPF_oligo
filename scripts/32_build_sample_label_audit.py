#!/usr/bin/env python
"""
Build a sample-level disease-label audit for all bulk mRNA and miRNA datasets.

The table is intended for Additional file 1 and records the original curated
label fields, inclusion status, final analysis label, exclusion reason, and
annotation source used for each sample.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
ANNOTATION = PROJECT_DIR / "metadata" / "all_bulk_mirna_annotation.csv"
OUT = PROJECT_DIR / "metadata" / "all_bulk_mirna_sample_label_audit.csv"
QC_OUT = PROJECT_DIR / "metadata" / "all_bulk_mirna_sample_label_audit_qc.csv"


def normalize_text(value: object) -> str:
    if pd.isna(value):
        return "NA"
    text = str(value).strip()
    return text if text else "NA"


def final_label(row: pd.Series) -> str:
    include = normalize_text(row.get("include")).lower()
    group = normalize_text(row.get("group"))
    if include != "yes":
        return "excluded"
    if group.lower() in {"ipf", "control"}:
        return group
    return group if group != "NA" else "included_unlabeled"


def exclusion_reason(row: pd.Series) -> str:
    include = normalize_text(row.get("include")).lower()
    group = normalize_text(row.get("group"))
    subgroup = normalize_text(row.get("subgroup"))
    notes = normalize_text(row.get("notes"))
    if include == "yes":
        return "included in analysis"
    parts = []
    if group != "NA":
        parts.append(f"group={group}")
    if subgroup != "NA" and subgroup != group:
        parts.append(f"subgroup={subgroup}")
    if notes != "NA":
        parts.append(f"notes={notes}")
    return "; ".join(parts) if parts else "excluded by curated annotation"


def main() -> None:
    df = pd.read_csv(ANNOTATION)
    rows = []
    for _, row in df.iterrows():
        original_label = " | ".join(
            [
                f"group={normalize_text(row.get('group'))}",
                f"subgroup={normalize_text(row.get('subgroup'))}",
                f"tissue={normalize_text(row.get('tissue'))}",
            ]
        )
        source_fields = "group; subgroup; include; tissue; notes; annotation_file; annotation_sheet"
        rows.append(
            {
                "series_id": normalize_text(row.get("series_id")),
                "sample_accession": normalize_text(row.get("geo_accession") or row.get("sample_id")),
                "sample_id": normalize_text(row.get("sample_id")),
                "data_type": normalize_text(row.get("data_type")),
                "original_disease_label": original_label,
                "included_excluded": "included" if normalize_text(row.get("include")).lower() == "yes" else "excluded",
                "final_analysis_label": final_label(row),
                "reason_for_exclusion_or_inclusion": exclusion_reason(row),
                "source_field_used_for_annotation": source_fields,
                "annotation_file": normalize_text(row.get("annotation_file")),
                "annotation_sheet": normalize_text(row.get("annotation_sheet")),
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT, index=False, encoding="utf-8-sig")

    qc = (
        audit.groupby(["series_id", "data_type", "included_excluded", "final_analysis_label"], dropna=False)
        .size()
        .reset_index(name="sample_count")
    )
    qc.to_csv(QC_OUT, index=False, encoding="utf-8-sig")
    print(OUT)
    print(QC_OUT)


if __name__ == "__main__":
    main()
