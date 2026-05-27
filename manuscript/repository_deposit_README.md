# Repository deposit package

This archive is prepared for public deposition of the IPF oligonucleotide actionability framework manuscript files.

Public repository: https://github.com/osbornzhou/IPF_oligo

## Contents

- `scripts/`: numbered analysis scripts used to generate processed results, figures, additional files, and manuscript artifacts.
- `metadata/`: curated dataset annotation and sample metadata used by the analysis workflow.
- `additional_files/`: machine-readable additional files and supplementary plot PDFs referenced by the manuscript, including AUC disease-state stress tests and miRNA target-program stress tests.
- `reproducibility_package/`: script manifest, input/output map, software environment notes, random seeds, and file manifest.
- `figures/`: high-resolution manuscript figure assets used for final submission checks.
- `ipf_oligo_ml_bmc_genomics_submission_ready.docx`: editable BMC Genomics submission manuscript.
- `ipf_oligo_ml_bmc_genomics_submission_ready.pdf`: reviewer preview PDF generated from the editable manuscript.
- `ipf_oligo_ml_bmc_genomics_manuscript_draft.md`: source manuscript text used by the document builder.
- `repository_deposit_README.md`: this deposit note.

## Reproducibility note

The analysis uses public GEO datasets and processed non-sensitive outputs only. Raw GEO matrices should be retrieved from the accessions listed in the manuscript and Additional file 1. The `reproducibility_package/README.md`, `script_input_output_map.csv`, `software_environment.csv`, and `random_seeds.csv` files map numbered scripts to manuscript figures, tables, and additional files, and record the software versions and random seeds used for reproducible regeneration of processed outputs.

For journal submission, this archive can be deposited as a public repository release and, if required, archived through Zenodo to obtain a persistent DOI.
