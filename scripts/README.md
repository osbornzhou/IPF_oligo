# Numbered analysis scripts

This folder contains the numbered scripts used for the submitted IPF multi-cohort transcriptomic actionability framework. The authoritative current script order is also recorded in:

`../manuscript/reproducibility_package/script_input_output_map.csv`

## Submitted script order

| Script | Purpose |
| --- | --- |
| `01_download_geo.py` | Download and organize GEO datasets. |
| `04_merge_annotations_and_check.py` | Merge curated annotations and perform sample-matching checks. |
| `05_extract_expression_matrices.py` | Extract expression matrices and generate expression-matrix QC. |
| `07_annotate_de_features.py` | Map probe, transcript, and miRNA feature identifiers. |
| `08_select_robust_candidates.py` | Select discovery-validation robust mRNA and miRNA candidates. |
| `09_integrate_mirna_mrna_targets.py` | Integrate robust candidates with miRTarBase opposite-direction axes. |
| `13_train_ml_models.py` | Train leakage-controlled machine-learning models and run external validation. |
| `15_single_cell_validation.py` | Generate single-cell candidate-localization summaries. |
| `21_submission_enhancements.py` | Build evidence-graded axes and integrated priority tables. |
| `22_mechanistic_extension_tnik.py` | Generate module, ligand-receptor proxy, perturbation-priority, and TNIK bridge outputs. |
| `24_machine_learning_sensitivity.py` | Run leave-one-validation-cohort and random robust-panel sensitivity analyses. |
| `26_oligonucleotide_actionability_index.py` | Build the oligonucleotide validation-planning score and Figure 7. |
| `27_single_cell_pseudobulk_validation.py` | Run donor-aware pseudobulk validation for core candidates. |
| `29_actionability_weight_sensitivity.py` | Run actionability weight-sensitivity analyses. |
| `30_mirna_exact_axis_support_table.py` | Build exact hsa-miR-375 axis support and cross-check tables. |
| `31_build_pseudobulk_figure.py` | Build high-resolution Figure 6C. |
| `32_build_sample_label_audit.py` | Build all-dataset sample-label audit tables. |
| `33_build_single_cell_broad_celltype_mapping.py` | Build broad-celltype mapping tables for pseudobulk audit. |
| `34_build_reproducibility_assets.py` | Build reproducibility metadata and checksums. |
| `35_strengthen_submission_package.py` | Add AUC confidence intervals and donor-level pseudobulk dotplots. |
| `36_boundary_and_mirna_program_support.py` | Add external disease-state boundary tests and miRNA target-program stress tests. |

## Manuscript and supplement builders

| Script | Purpose |
| --- | --- |
| `23_build_bmc_submission_ready_docx.py` | Build the editable BMC Genomics submission DOCX from the manuscript markdown and figure assets. |
| `25_build_bmc_additional_files.py` | Build BMC-style machine-readable Additional files 1-14. |

## Notes for reviewers

- Script numbering is historical, so not every integer is used.
- Scripts are designed to consume locally available public GEO files and generated intermediate outputs.
- Random seeds are recorded in `../manuscript/reproducibility_package/random_seeds.csv`.
- Software notes are recorded in `../manuscript/reproducibility_package/software_environment.csv`.
- Input/output mapping is recorded in `../manuscript/reproducibility_package/script_input_output_map.csv`.
