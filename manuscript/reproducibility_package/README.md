# Reproducibility package

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
