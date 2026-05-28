#!/usr/bin/env python
"""
Build reproducibility assets for the manuscript package.

This script creates compact tables that make the custom parts of the study
auditable: validation-support rules by dataset, machine-learning hyperparameter
specifications, random seeds, script-to-output mapping, software notes, and a
manifest for key generated files.
"""

from __future__ import annotations

import hashlib
import importlib.metadata as importlib_metadata
import re
import sys
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
REPRO_DIR = PROJECT_DIR / "manuscript" / "reproducibility_package"
ROBUST_DIR = PROJECT_DIR / "results" / "robust_candidates"
ML_DIR = PROJECT_DIR / "results" / "models" / "ml_outputs_discovery_only_mrna"
ML_SENS_DIR = PROJECT_DIR / "results" / "ml_sensitivity"
ADDITIONAL_DIR = PROJECT_DIR / "manuscript" / "additional_files"

REPRO_DIR.mkdir(parents=True, exist_ok=True)
ROBUST_DIR.mkdir(parents=True, exist_ok=True)
ML_DIR.mkdir(parents=True, exist_ok=True)
ML_SENS_DIR.mkdir(parents=True, exist_ok=True)


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT_DIR)).replace("\\", "/")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def package_version(package: str) -> str:
    try:
        return importlib_metadata.version(package)
    except importlib_metadata.PackageNotFoundError:
        return "not installed in bundled Python; used in original analysis environment where applicable"


def sklearn_version_from_artifact() -> str:
    artifact = ML_DIR / "artifacts" / "elastic_net_final_model.joblib"
    if not artifact.exists():
        return "not detected; final model artifact unavailable"
    text = artifact.read_bytes()
    match = re.search(rb"_sklearn_version\x94\x8c\x05([^\x94]+)\x94", text)
    if match:
        return match.group(1).decode("ascii", errors="replace")
    return "not detected in final model artifact"


def build_validation_support_rule() -> None:
    matrix = pd.read_csv(PROJECT_DIR / "metadata" / "matrix_processing_decisions.csv")
    validation_roles = {
        "GSE110147": ("bulk mRNA", "mRNA validation", "FDR+abs_logFC"),
        "GSE150910": ("bulk mRNA", "mRNA validation", "FDR+abs_logFC"),
        "GSE53845": ("bulk mRNA", "mRNA validation", "direction_supporting_only"),
        "GSE92592": ("bulk mRNA", "mRNA validation", "FDR+abs_logFC"),
        "GSE21394": ("miRNA", "miRNA validation", "FDR+abs_logFC"),
        "GSE27430": ("miRNA", "miRNA validation", "direction_supporting_only"),
        "GSE32537": ("bulk mRNA", "mRNA discovery", "discovery_FDR+abs_logFC"),
        "GSE32538": ("miRNA", "miRNA discovery", "discovery_FDR+abs_logFC"),
    }
    rows = []
    for series_id, (data_type, role, support_rule) in validation_roles.items():
        row = matrix[matrix["series_id"].eq(series_id)]
        if row.empty:
            input_type = "not found in matrix_processing_decisions"
            method = "NA"
            note = "NA"
        else:
            record = row.iloc[0]
            input_type = record.get("input_matrix_type", "NA")
            method = record.get("differential_expression_method", "NA")
            note = record.get("validation_threshold_note", "NA")
        rows.append(
            {
                "series_id": series_id,
                "data_type": data_type,
                "dataset_role": role,
                "input_matrix_type": input_type,
                "differential_expression_method": method,
                "validation_support_rule": support_rule,
                "rule_interpretation": (
                    "FDR and absolute logFC were used for log-scale/count-derived datasets"
                    if support_rule == "FDR+abs_logFC"
                    else "Centered/transformed datasets were used as direction-supporting validation layers with cautious effect-size interpretation"
                    if support_rule == "direction_supporting_only"
                    else "Discovery dataset threshold used to define the discovery candidate universe"
                ),
                "source_threshold_note": note,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(ROBUST_DIR / "validation_support_rule_by_dataset.csv", index=False, encoding="utf-8-sig")
    qc = pd.DataFrame(
        [
            {
                "qc_item": "all_expected_series_present",
                "value": out["series_id"].nunique(),
                "pass": out["series_id"].nunique() == len(validation_roles),
            },
            {
                "qc_item": "centered_transformed_flagged_direction_supporting",
                "value": ";".join(out.loc[out["validation_support_rule"].eq("direction_supporting_only"), "series_id"]),
                "pass": set(out.loc[out["validation_support_rule"].eq("direction_supporting_only"), "series_id"])
                == {"GSE53845", "GSE27430"},
            },
        ]
    )
    qc.to_csv(ROBUST_DIR / "validation_support_rule_qc.csv", index=False, encoding="utf-8-sig")


def build_ml_reproducibility_tables() -> None:
    grid_rows = [
        {"model": "lasso_logistic", "search_parameter": "select__k", "values": "15;30;50;100 where permitted"},
        {"model": "lasso_logistic", "search_parameter": "clf__C", "values": "0.03;0.1;0.3;1.0"},
        {"model": "elastic_net", "search_parameter": "select__k", "values": "15;30;50;100 where permitted"},
        {"model": "elastic_net", "search_parameter": "clf__C", "values": "0.03;0.1;0.3;1.0"},
        {"model": "elastic_net", "search_parameter": "clf__l1_ratio", "values": "0.2;0.5;0.8"},
        {"model": "linear_svm", "search_parameter": "select__k", "values": "15;30;50;100 where permitted"},
        {"model": "linear_svm", "search_parameter": "clf__C", "values": "0.03;0.1;0.3;1.0"},
        {"model": "rbf_svm", "search_parameter": "select__k", "values": "30;50 where permitted"},
        {"model": "rbf_svm", "search_parameter": "clf__C", "values": "0.1;1.0;3.0"},
        {"model": "rbf_svm", "search_parameter": "clf__gamma", "values": "scale"},
        {"model": "random_forest", "search_parameter": "select__k", "values": "50;100 where permitted"},
        {"model": "random_forest", "search_parameter": "clf__max_depth", "values": "3;5;None"},
        {"model": "random_forest", "search_parameter": "clf__min_samples_leaf", "values": "1;3"},
        {"model": "gradient_boosting", "search_parameter": "select__k", "values": "30;50 where permitted"},
        {"model": "gradient_boosting", "search_parameter": "clf__n_estimators", "values": "80;150"},
        {"model": "gradient_boosting", "search_parameter": "clf__learning_rate", "values": "0.03;0.08"},
        {"model": "gradient_boosting", "search_parameter": "clf__max_depth", "values": "1;2"},
        {"model": "mlp_small", "search_parameter": "select__k", "values": "30;50 where permitted"},
        {"model": "mlp_small", "search_parameter": "clf__hidden_layer_sizes", "values": "(16);(16,8);(32,16)"},
        {"model": "mlp_small", "search_parameter": "clf__alpha", "values": "0.001;0.01;0.1"},
    ]
    pd.DataFrame(grid_rows).to_csv(ML_DIR / "ml_hyperparameter_grid_summary.csv", index=False, encoding="utf-8-sig")

    final_spec = pd.DataFrame(
        [
            {"field": "final_model", "value": "Elastic Net logistic regression"},
            {"field": "training_dataset", "value": "GSE32537"},
            {"field": "external_validation_datasets", "value": "GSE110147;GSE150910;GSE53845;GSE92592"},
            {"field": "feature_universe", "value": "329 common discovery-only mRNA features"},
            {"field": "final_selected_panel_size", "value": "25 genes"},
            {"field": "imputation", "value": "SimpleImputer(strategy='median')"},
            {"field": "scaling", "value": "StandardScaler"},
            {"field": "feature_selection", "value": "SelectKBest(f_classif), final k=50 in final fitted artifact"},
            {"field": "classifier", "value": "LogisticRegression(penalty='elasticnet', solver='saga', class_weight='balanced', max_iter=8000)"},
            {"field": "final_C", "value": "0.03"},
            {"field": "final_l1_ratio", "value": "0.2"},
            {"field": "random_state", "value": "20260524 for main ML; 20260525 for random-panel sensitivity"},
            {"field": "outer_cv", "value": "RepeatedStratifiedKFold, 5 folds x 10 repeats"},
            {"field": "inner_cv", "value": "StratifiedKFold, 3 folds for nested tuning; 5 folds for final refit tuning"},
            {"field": "final_model_artifact", "value": "results/models/ml_outputs_discovery_only_mrna/artifacts/elastic_net_final_model.joblib"},
        ]
    )
    final_spec.to_csv(ML_DIR / "ml_final_elastic_net_specification.csv", index=False, encoding="utf-8-sig")

    notes = pd.DataFrame(
        [
            {
                "item": "leakage_control",
                "note": "External validation cohorts were not used for feature selection, hyperparameter tuning, or final refitting.",
            },
            {
                "item": "reported_gene_set_benchmarking",
                "note": "Published comparators were retrained reported gene sets under the same framework, not original model-object reproductions.",
            },
            {
                "item": "decision_curve_outputs",
                "note": "Calibration and decision-curve outputs were retained as exploratory model-review checks, not as clinical threshold claims.",
            },
        ]
    )
    notes.to_csv(ML_DIR / "ml_reproducibility_notes.csv", index=False, encoding="utf-8-sig")


def build_clean_ml_sensitivity_aliases() -> None:
    expected_files = [
        "elastic_net_external_per_cohort_performance.csv",
        "leave_one_validation_cohort_sensitivity.csv",
        "random_robust_25_gene_panel_baseline_summary.csv",
        "random_robust_25_gene_panel_baseline.csv",
        "random_robust_25_gene_panel_cohort_performance.csv",
        "machine_learning_sensitivity_qc.csv",
    ]
    rows = []
    for clean_name in expected_files:
        dst = ML_SENS_DIR / clean_name
        if dst.exists():
            df = pd.read_csv(dst)
            rows.append(
                {
                    "clean_file": rel(dst),
                    "status": "included",
                    "rows": len(df),
                    "source_note": "generated from locked Elastic Net predictions or random robust-panel sensitivity outputs",
                }
            )
        else:
            rows.append({"clean_file": rel(dst), "status": "missing", "rows": 0, "source_note": "expected machine-learning sensitivity output was not found"})
    pd.DataFrame(rows).to_csv(ML_SENS_DIR / "machine_learning_sensitivity_manifest.csv", index=False, encoding="utf-8-sig")


def build_reproducibility_package() -> None:
    script_rows = [
        ("01_download_geo.py", "Download and organize GEO datasets", "data_raw/GEO", "dataset folders"),
        ("04_merge_annotations_and_check.py", "Merge curated annotation and sample matching", "data_raw/GEO annotations", "metadata/all_bulk_mirna_annotation.csv"),
        ("05_extract_expression_matrices.py", "Extract expression matrices and matrix QC", "GEO matrix files", "metadata/expression_matrix_qc.csv"),
        ("07_annotate_de_features.py", "Map probes/transcripts/miRNA probes", "DE feature IDs", "results/feature_annotation"),
        ("08_select_robust_candidates.py", "Discovery-validation robust candidate selection", "DE results", "results/robust_candidates"),
        ("09_integrate_mirna_mrna_targets.py", "miRTarBase opposite-direction axis integration", "robust mRNA/miRNA; miRTarBase", "results/mirna_mrna_axes"),
        ("13_train_ml_models.py", "Leakage-controlled model training and external validation", "feature matrices", "results/models/ml_outputs_discovery_only_mrna"),
        ("15_single_cell_validation.py", "Single-cell localization summaries", "single-cell matrices/annotations", "results/single_cell_validation"),
        ("21_submission_enhancements.py", "Evidence grading and integrated priority table", "robust/ML/single-cell outputs", "results/submission_enhancements"),
        ("22_mechanistic_extension_tnik.py", "Module, ligand-receptor proxy, perturbation-priority, TNIK bridge", "bulk and single-cell outputs", "results/mechanistic_extension"),
        ("24_machine_learning_sensitivity.py", "ML leave-one-cohort and random robust-panel sensitivity", "locked Elastic Net outputs", "results/ml_sensitivity"),
        ("26_oligonucleotide_actionability_index.py", "Perturbation-triage validation-planning score and triage-map figure", "integrated priority outputs", "results/oligonucleotide_actionability"),
        ("27_single_cell_pseudobulk_validation.py", "Donor-aware pseudobulk validation for core candidates", "single-cell matrices/annotations", "results/single_cell_pseudobulk"),
        ("29_actionability_weight_sensitivity.py", "Perturbation-triage score weight-sensitivity analysis", "perturbation-triage score", "results/oligonucleotide_actionability"),
        ("30_mirna_exact_axis_support_table.py", "Exact hsa-miR-375 axis cross-check table", "miRNA axes", "results/submission_enhancements"),
        ("31_build_pseudobulk_figure.py", "High-resolution donor-aware pseudobulk figure panel", "pseudobulk summaries", "results/single_cell_pseudobulk"),
        ("32_build_sample_label_audit.py", "All-dataset sample-label audit", "curated annotations", "metadata/all_bulk_mirna_sample_label_audit.csv"),
        ("33_build_single_cell_broad_celltype_mapping.py", "Broad-celltype mapping for pseudobulk", "single-cell annotations", "results/single_cell_pseudobulk"),
        ("34_build_reproducibility_assets.py", "Reproducibility package", "project outputs", "manuscript/reproducibility_package"),
        ("35_strengthen_submission_package.py", "AUC confidence intervals, pseudobulk donor dotplots, and final submission-support notes", "locked predictions and pseudobulk outputs", "results/models/ml_outputs_discovery_only_mrna; results/single_cell_pseudobulk"),
        ("36_boundary_and_mirna_program_support.py", "External disease-state boundary tests and miRNA target-program stress tests", "feature matrices, locked predictions, miRTarBase targets", "results/ml_stress_tests; results/mirna_program_support"),
        ("37_build_editable_main_figures.py", "Editable SVG source figures for main Figures 1-6", "manuscript figure data tables", "manuscript/editable_figures"),
    ]
    pd.DataFrame(script_rows, columns=["script", "purpose", "primary_input", "primary_output"]).to_csv(
        REPRO_DIR / "script_input_output_map.csv", index=False, encoding="utf-8-sig"
    )

    seed_rows = [
        ("13_train_ml_models.py", "main machine-learning CV and model refit", 20260524),
        ("24_machine_learning_sensitivity.py", "random robust-panel baseline", 20260525),
        ("29_actionability_weight_sensitivity.py", "seeded +/-20% perturbation-triage score weight perturbations", 20260526),
        ("35_strengthen_submission_package.py", "stratified bootstrap confidence intervals for external validation metrics", 20260526),
        ("36_boundary_and_mirna_program_support.py", "matched random discovery-feature panels and permutation target-release score", 20260527),
    ]
    pd.DataFrame(seed_rows, columns=["script", "random_process", "seed"]).to_csv(REPRO_DIR / "random_seeds.csv", index=False)

    packages = ["pandas", "numpy", "openpyxl", "matplotlib", "Pillow", "scipy", "pypdf"]
    software_rows = [{"software": "Python", "version_or_note": sys.version.split()[0]}]
    software_rows.extend({"software": p, "version_or_note": package_version(p)} for p in packages)
    software_rows.extend(
        [
            {"software": "scikit-learn", "version_or_note": sklearn_version_from_artifact()},
            {"software": "joblib", "version_or_note": "used to serialize final model artifacts; exact package version not encoded separately in artifacts"},
            {"software": "R", "version_or_note": "4.6.0"},
            {"software": "limma", "version_or_note": "3.68.3"},
            {"software": "edgeR", "version_or_note": "4.10.0"},
            {"software": "clusterProfiler", "version_or_note": "4.20.0"},
            {"software": "ReactomePA", "version_or_note": "1.56.0"},
            {"software": "org.Hs.eg.db", "version_or_note": "3.23.1"},
        ]
    )
    pd.DataFrame(software_rows).to_csv(REPRO_DIR / "software_environment.csv", index=False, encoding="utf-8-sig")

    manifest_targets = [
        PROJECT_DIR / "manuscript" / "ipf_oligo_ml_bmc_genomics_manuscript_draft.md",
        ML_DIR / "ml_final_elastic_net_specification.csv",
        ML_DIR / "ml_hyperparameter_grid_summary.csv",
        ML_DIR / "ml_external_validation_bootstrap_ci.csv",
        ML_DIR / "ml_metric_definitions_and_interpretation_notes.csv",
        ROBUST_DIR / "validation_support_rule_by_dataset.csv",
        PROJECT_DIR / "results" / "submission_enhancements" / "mirna_exact_axis_support_table.csv",
        PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_pseudobulk_differential_summary.csv",
        PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_broad_celltype_mapping.csv",
        PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_pseudobulk_broad_celltype_cell_count_summary.csv",
        PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_pseudobulk_method_notes.csv",
        PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_pseudobulk_core_candidate_donor_dotplots.pdf",
        PROJECT_DIR / "results" / "ml_stress_tests" / "gse110147_excluded_ild_disease_control_summary.csv",
        PROJECT_DIR / "results" / "ml_stress_tests" / "matched_random_discovery_feature_panel_summary.csv",
        PROJECT_DIR / "results" / "ml_stress_tests" / "cohort_adjusted_external_logistic_association.csv",
        PROJECT_DIR / "results" / "ml_stress_tests" / "additional_figure_s13_auc_stress_tests.pdf",
        PROJECT_DIR / "results" / "mirna_program_support" / "robust_mirna_target_set_enrichment.csv",
        PROJECT_DIR / "results" / "mirna_program_support" / "hsa_mir_375_target_repression_release_score.csv",
        PROJECT_DIR / "results" / "mirna_program_support" / "additional_figure_s14_mirna_program_stress_test.pdf",
        PROJECT_DIR / "results" / "oligonucleotide_actionability" / "actionability_weight_sensitivity_rank_summary.csv",
        PROJECT_DIR / "manuscript" / "editable_figures" / "README.md",
        PROJECT_DIR / "manuscript" / "editable_figures" / "editable_figure_manifest.csv",
    ]
    manifest_targets.extend(sorted((PROJECT_DIR / "manuscript" / "editable_figures" / "svg").glob("Figure_*_editable.svg")))
    manifest_targets.extend(sorted(ADDITIONAL_DIR.glob("Additional_file_*.xlsx")))
    manifest_targets.extend(sorted(ADDITIONAL_DIR.glob("Additional_file_*.zip")))
    manifest = []
    for path in manifest_targets:
        if path.exists() and path.is_file():
            manifest.append(
                {
                    "path": rel(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    pd.DataFrame(manifest).to_csv(REPRO_DIR / "file_manifest.csv", index=False, encoding="utf-8-sig")

    readme = """# Reproducibility package

This folder indexes the machine-readable materials used to regenerate and audit
the IPF oligonucleotide-focused perturbation-triage Research article.

## Contents

- `script_input_output_map.csv`: numbered script order with the primary input and output of each step.
- `software_environment.csv`: software and package notes for reproducing the analyses.
- `random_seeds.csv`: random seeds used in model training and sensitivity analyses.
- `file_manifest.csv`: checksums for key manuscript, result, and additional-file artifacts available at package-build time.

The main manuscript data tables are provided as Additional files 1-10, the
numbered scripts and reproducibility tables are bundled in Additional file 11,
and supplementary PDF plots are provided as Additional files 12-14. Editable
SVG sources for main Figures 1-6 are provided in `manuscript/editable_figures`
for post-generation alignment and label edits. Public GEO accessions are listed
in the manuscript and in Additional file 1. This repository is publicly
available at https://github.com/osbornzhou/IPF_oligo. A Zenodo DOI can be added
if the repository is archived as a release.
"""
    (REPRO_DIR / "README.md").write_text(readme, encoding="utf-8")


def build_pseudobulk_cell_count_summary() -> None:
    source = PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_pseudobulk_donor_celltype_expression.csv"
    if not source.exists():
        return
    df = pd.read_csv(source)
    donor_celltype = df[["series_id", "donor_id", "group", "broad_celltype", "cells"]].drop_duplicates()
    summary = (
        donor_celltype.groupby(["series_id", "broad_celltype", "group"], as_index=False)
        .agg(
            donors=("donor_id", "nunique"),
            total_cells=("cells", "sum"),
            median_cells_per_donor=("cells", "median"),
            min_cells_per_donor=("cells", "min"),
            max_cells_per_donor=("cells", "max"),
        )
        .sort_values(["series_id", "broad_celltype", "group"])
    )
    summary.to_csv(
        PROJECT_DIR / "results" / "single_cell_pseudobulk" / "single_cell_pseudobulk_broad_celltype_cell_count_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )


def main() -> None:
    build_validation_support_rule()
    build_ml_reproducibility_tables()
    build_clean_ml_sensitivity_aliases()
    build_pseudobulk_cell_count_summary()
    build_reproducibility_package()
    print(REPRO_DIR)


if __name__ == "__main__":
    main()
