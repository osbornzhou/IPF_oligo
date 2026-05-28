#!/usr/bin/env python
"""Prepare unified single-cell sample and cell annotations.

Outputs:
- metadata/single_cell_sample_annotation.csv
- metadata/single_cell_cell_annotation.csv.gz
- metadata/single_cell/by_dataset/*_sample_annotation.csv
- metadata/single_cell/by_dataset/*_cell_annotation.csv.gz
"""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GEO_ROOT = PROJECT_ROOT / "data_raw" / "GEO"
OUT_ROOT = PROJECT_ROOT / "metadata" / "single_cell"
BY_DATASET = OUT_ROOT / "by_dataset"


SAMPLE_COLUMNS = [
    "sample_id",
    "series_id",
    "geo_accession",
    "sample_title",
    "data_type",
    "dataset_role",
    "tissue",
    "organism",
    "group",
    "subgroup",
    "include",
    "platform",
    "library_technology",
    "library_id",
    "subject_id",
    "batch",
    "source",
    "n_cells",
    "cell_type_major_count",
    "cell_type_minor_count",
    "notes",
]

CELL_COLUMNS = [
    "cell_id",
    "series_id",
    "sample_id",
    "subject_id",
    "library_id",
    "geo_accession",
    "group",
    "subgroup",
    "include",
    "platform",
    "library_technology",
    "cell_type_major",
    "cell_type_minor",
    "cell_type_subclass",
    "nUMI",
    "nGene",
    "percent_mt",
    "batch",
    "source",
    "notes",
]


def include_for_group(group: str) -> str:
    return "Yes" if group in {"IPF", "Control"} else "No"


def map_diagnosis_to_group(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"ipf", "idiopathic pulmonary fibrosis"}:
        return "IPF"
    if text in {"control", "donor", "normal"}:
        return "Control"
    return "Exclude"


def clean_subgroup(value: object) -> str:
    text = str(value).strip()
    mapping = {
        "ipf": "IPF",
        "control": "Control",
        "donor": "Control",
        "chp": "cHP",
        "copd": "COPD",
        "idiopathic pulmonary fibrosis": "IPF",
        "hypersensitivity pneumonitis": "HP",
        "systemic slcerosis-associated interstitial lung disease": "SSc-ILD",
        "myositis-associated interstitial lng disease": "Myositis-ILD",
    }
    return mapping.get(text.lower(), text)


def join_unique(values: pd.Series) -> str:
    clean = [str(x) for x in values.dropna().unique() if str(x) != "nan"]
    return ";".join(sorted(clean)) if clean else pd.NA


def write_csv(df: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8", compression=compression)


def read_gzip_lines(path: Path, n: int | None = None) -> list[str]:
    lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for i, line in enumerate(handle):
            if n is not None and i >= n:
                break
            lines.append(line.rstrip("\n"))
    return lines


def parse_series_matrix_fields(path: Path) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                fields.setdefault(parts[0], parts[1:])
    return fields


def gse135893_platform_map() -> pd.DataFrame:
    rows = []
    matrix_root = GEO_ROOT / "GSE135893" / "matrix"
    for matrix_path in sorted(matrix_root.glob("GSE135893-GPL*_series_matrix.txt.gz")):
        fields = parse_series_matrix_fields(matrix_path)
        titles = fields.get("!Sample_title", [])
        geos = fields.get("!Sample_geo_accession", [pd.NA] * len(titles))
        platforms = fields.get("!Sample_platform_id", [pd.NA] * len(titles))
        for title, geo, platform in zip(titles, geos, platforms):
            sample_name = str(title).split(",")[-1].strip()
            rows.append(
                {
                    "sample_title_geo": title,
                    "geo_accession": geo,
                    "platform": platform,
                    "sample_name_from_title": sample_name,
                    "sample_name_base": sample_name.split("-")[0],
                }
            )
    return pd.DataFrame(rows)


def validate_cell_ids(cell_ids: pd.Series, barcode_path: Path) -> dict[str, object]:
    barcode_ids = set(read_gzip_lines(barcode_path))
    ids = set(cell_ids.astype(str))
    return {
        "barcode_count": len(barcode_ids),
        "annotation_count": len(ids),
        "missing_in_annotation": len(barcode_ids - ids),
        "missing_in_barcodes": len(ids - barcode_ids),
    }


def prepare_gse135893() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    series_id = "GSE135893"
    meta_path = GEO_ROOT / series_id / "supplementary" / "GSE135893_IPF_metadata.csv.gz"
    barcode_path = GEO_ROOT / series_id / "supplementary" / "GSE135893_barcodes.tsv.gz"
    raw = pd.read_csv(meta_path)
    platform_map = gse135893_platform_map()
    platform_by_library = platform_map.set_index("sample_name_from_title")["platform"].to_dict()
    platform_by_sample = platform_map.drop_duplicates("sample_name_base").set_index("sample_name_base")["platform"].to_dict()
    geo_by_sample = platform_map.drop_duplicates("sample_name_base").set_index("sample_name_base")["geo_accession"].to_dict()

    cell = pd.DataFrame(
        {
            "cell_id": raw["Unnamed: 0"].astype(str),
            "series_id": series_id,
            # Use orig.ident as sample_id because it is the library/barcode prefix
            # used by the sparse matrix barcodes, e.g. F00431_AAAC....
            "sample_id": raw["orig.ident"].astype(str),
            "subject_id": raw["Sample_Name"].astype(str),
            "library_id": raw["orig.ident"].astype(str),
            "geo_accession": raw["Sample_Name"].astype(str).map(geo_by_sample),
            "group": raw["Diagnosis"].map(map_diagnosis_to_group),
            "subgroup": raw["Diagnosis"].map(clean_subgroup),
            "include": raw["Diagnosis"].map(map_diagnosis_to_group).map(include_for_group),
            "platform": raw["orig.ident"].astype(str).map(platform_by_library).fillna(
                raw["Sample_Name"].astype(str).map(platform_by_sample)
            ),
            "library_technology": "10x Genomics 5' GEX",
            "cell_type_major": raw["population"].astype(str),
            "cell_type_minor": raw["celltype"].astype(str),
            "cell_type_subclass": raw["celltype"].astype(str),
            "nUMI": raw["nCount_RNA"],
            "nGene": raw["nFeature_RNA"],
            "percent_mt": raw["percent.mt"],
            "batch": raw["Sample_Source"].astype(str),
            "source": raw["Sample_Source"].astype(str),
            "notes": pd.NA,
        }
    )[CELL_COLUMNS]
    cell["platform"] = cell["platform"].fillna("NA")
    cell["geo_accession"] = cell["geo_accession"].fillna("NA")

    grouped = cell.groupby("sample_id", dropna=False)
    sample = grouped.agg(
        series_id=("series_id", "first"),
        geo_accession=("geo_accession", "first"),
        sample_title=("subject_id", join_unique),
        data_type=("series_id", lambda _: "scRNA-seq"),
        dataset_role=("series_id", lambda _: "single_cell_validation"),
        tissue=("series_id", lambda _: "Lung tissue"),
        organism=("series_id", lambda _: "Homo sapiens"),
        group=("group", join_unique),
        subgroup=("subgroup", join_unique),
        include=("include", lambda x: "Yes" if set(x) <= {"Yes"} else "No"),
        platform=("platform", join_unique),
        library_technology=("library_technology", join_unique),
        library_id=("library_id", join_unique),
        subject_id=("subject_id", join_unique),
        batch=("batch", join_unique),
        source=("source", join_unique),
        n_cells=("cell_id", "size"),
        cell_type_major_count=("cell_type_major", "nunique"),
        cell_type_minor_count=("cell_type_minor", "nunique"),
        notes=("notes", "first"),
    ).reset_index()
    sample = sample[SAMPLE_COLUMNS]

    validation = validate_cell_ids(cell["cell_id"], barcode_path)
    validation["series_id"] = series_id
    validation["n_cells"] = len(cell)
    validation["n_samples"] = sample.shape[0]
    return sample, cell, validation


def prepare_gse136831() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    series_id = "GSE136831"
    meta_path = GEO_ROOT / series_id / "supplementary" / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
    barcode_path = GEO_ROOT / series_id / "supplementary" / "GSE136831_AllCells.cellBarcodes.txt.gz"
    raw = pd.read_csv(meta_path, sep="\t")

    disease = raw["Disease_Identity"].map(clean_subgroup)
    group = disease.map(map_diagnosis_to_group)
    cell = pd.DataFrame(
        {
            "cell_id": raw["CellBarcode_Identity"].astype(str),
            "series_id": series_id,
            "sample_id": raw["Library_Identity"].astype(str),
            "subject_id": raw["Subject_Identity"].astype(str),
            "library_id": raw["Library_Identity"].astype(str),
            "geo_accession": pd.NA,
            "group": group,
            "subgroup": disease,
            "include": group.map(include_for_group),
            "platform": "GPL20301",
            "library_technology": "10x Genomics scRNA-seq",
            "cell_type_major": raw["CellType_Category"].astype(str),
            "cell_type_minor": raw["Manuscript_Identity"].astype(str),
            "cell_type_subclass": raw["Subclass_Cell_Identity"].astype(str),
            "nUMI": raw["nUMI"],
            "nGene": raw["nGene"],
            "percent_mt": pd.NA,
            "batch": raw["Library_Identity"].astype(str),
            "source": pd.NA,
            "notes": pd.NA,
        }
    )[CELL_COLUMNS]

    grouped = cell.groupby("sample_id", dropna=False)
    sample = grouped.agg(
        series_id=("series_id", "first"),
        geo_accession=("geo_accession", "first"),
        sample_title=("sample_id", "first"),
        data_type=("series_id", lambda _: "scRNA-seq"),
        dataset_role=("series_id", lambda _: "single_cell_validation"),
        tissue=("series_id", lambda _: "Lung tissue"),
        organism=("series_id", lambda _: "Homo sapiens"),
        group=("group", join_unique),
        subgroup=("subgroup", join_unique),
        include=("include", lambda x: "Yes" if set(x) <= {"Yes"} else "No"),
        platform=("platform", join_unique),
        library_technology=("library_technology", join_unique),
        library_id=("library_id", join_unique),
        subject_id=("subject_id", join_unique),
        batch=("batch", join_unique),
        source=("source", "first"),
        n_cells=("cell_id", "size"),
        cell_type_major_count=("cell_type_major", "nunique"),
        cell_type_minor_count=("cell_type_minor", "nunique"),
        notes=("notes", "first"),
    ).reset_index()
    sample = sample[SAMPLE_COLUMNS]

    validation = validate_cell_ids(cell["cell_id"], barcode_path)
    validation["series_id"] = series_id
    validation["n_cells"] = len(cell)
    validation["n_samples"] = sample.shape[0]
    return sample, cell, validation


def parse_gse122960_sample_metadata() -> pd.DataFrame:
    series_id = "GSE122960"
    matrix_path = GEO_ROOT / series_id / "matrix" / "GSE122960_series_matrix.txt.gz"
    wanted = {
        "!Sample_title": None,
        "!Sample_geo_accession": None,
        "!Sample_source_name_ch1": None,
        "!Sample_characteristics_ch1": None,
        "!Sample_platform_id": None,
    }
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            key = line.split("\t", 1)[0]
            if key in wanted:
                wanted[key] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                break

    titles = wanted["!Sample_title"] or []
    geos = wanted["!Sample_geo_accession"] or [pd.NA] * len(titles)
    sources = wanted["!Sample_source_name_ch1"] or ["Lung"] * len(titles)
    conditions = wanted["!Sample_characteristics_ch1"] or [pd.NA] * len(titles)
    platforms = wanted["!Sample_platform_id"] or ["GPL20301"] * len(titles)

    rows = []
    for title, geo, source, condition, platform in zip(titles, geos, sources, conditions, platforms):
        diagnosis = re.sub(r"^disease condition:\s*", "", str(condition), flags=re.I).strip()
        subgroup = clean_subgroup(diagnosis)
        group = map_diagnosis_to_group(subgroup)
        rows.append(
            {
                # GSE122960 stores one HDF5 per GSM accession. Use GSM as the
                # stable sample_id; keep the human-readable title separately.
                "sample_id": geo,
                "series_id": series_id,
                "geo_accession": geo,
                "sample_title": title,
                "data_type": "scRNA-seq",
                "dataset_role": "single_cell_validation_backup",
                "tissue": source,
                "organism": "Homo sapiens",
                "group": group,
                "subgroup": subgroup,
                "include": include_for_group(group),
                "platform": platform,
                "library_technology": "10x Genomics scRNA-seq",
                "library_id": geo,
                "subject_id": title,
                "batch": pd.NA,
                "source": pd.NA,
                "n_cells": pd.NA,
                "cell_type_major_count": pd.NA,
                "cell_type_minor_count": pd.NA,
                "notes": "Local cell-level matrix/metadata not available; GEO points to processed HDF5 in RAW.tar.",
            }
        )
    return pd.DataFrame(rows, columns=SAMPLE_COLUMNS)


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    BY_DATASET.mkdir(parents=True, exist_ok=True)

    sample_frames = []
    cell_frames = []
    validations = []

    for series_id, prepare in [
        ("GSE135893", prepare_gse135893),
        ("GSE136831", prepare_gse136831),
    ]:
        sample, cell, validation = prepare()
        sample_frames.append(sample)
        cell_frames.append(cell)
        validations.append(validation)
        write_csv(sample, BY_DATASET / f"{series_id}_sample_annotation.csv")
        write_csv(cell, BY_DATASET / f"{series_id}_cell_annotation.csv.gz", compression="gzip")

    gse122960_sample = parse_gse122960_sample_metadata()
    sample_frames.append(gse122960_sample)
    write_csv(gse122960_sample, BY_DATASET / "GSE122960_sample_annotation.csv")

    all_samples = pd.concat(sample_frames, ignore_index=True)
    all_cells = pd.concat(cell_frames, ignore_index=True)
    validation_df = pd.DataFrame(validations)

    write_csv(all_samples, OUT_ROOT / "single_cell_sample_annotation.csv")
    write_csv(all_cells, OUT_ROOT / "single_cell_cell_annotation.csv.gz", compression="gzip")
    write_csv(validation_df, OUT_ROOT / "single_cell_annotation_validation.csv")

    print("Wrote", OUT_ROOT / "single_cell_sample_annotation.csv", all_samples.shape)
    print("Wrote", OUT_ROOT / "single_cell_cell_annotation.csv.gz", all_cells.shape)
    print("Wrote", OUT_ROOT / "single_cell_annotation_validation.csv", validation_df.shape)
    print("\nSample summary:")
    print(pd.crosstab(all_samples["series_id"], all_samples["group"], dropna=False).to_string())
    print("\nCell summary:")
    print(pd.crosstab(all_cells["series_id"], all_cells["group"], dropna=False).to_string())


if __name__ == "__main__":
    main()
