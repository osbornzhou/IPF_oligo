# IPF boundary-tested transcriptomic perturbation-triage study

This repository contains reproducibility materials for the BMC Genomics Research article:

**A boundary-tested transcriptomic perturbation-triage framework identifies oligonucleotide validation candidates and fibrotic disease-state markers in idiopathic pulmonary fibrosis**

The project analyzes public human idiopathic pulmonary fibrosis (IPF) bulk mRNA, miRNA, and single-cell transcriptomic datasets to triage reproducible IPF/fibrotic transcriptomic abnormalities for oligonucleotide-focused validation planning. The workflow is a multi-cohort transcriptomic study, not a clinical diagnostic deployment workflow and not a therapeutic-candidate validation study.

## Current submitted workflow

The submitted workflow includes:

- public GEO dataset curation and sample-label audit;
- expression-matrix QC and matrix-type adjudication;
- differential expression with limma or edgeR-limma voom as appropriate to matrix type;
- cross-cohort robust mRNA and miRNA candidate selection;
- evidence-graded miRTarBase miRNA-mRNA candidate axes;
- leakage-controlled machine learning with external validation;
- machine-learning boundary stress tests, including disease-control scoring, matched random discovery-feature panels, non-perfect-cohort summaries, and cohort-adjusted score tests;
- single-cell localization and donor-aware pseudobulk validation;
- exploratory coexpression, curated ligand-receptor, perturbation-priority, and TNIK bridge analyses;
- oligonucleotide-focused perturbation-triage scoring and score-sensitivity checks;
- machine-readable additional files, manuscript figures, and editable SVG main-figure sources.

## Repository layout

```text
ipf_oligo_ml/
  data_raw/                  Public GEO source files where available locally
  data_processed/            Intermediate processed matrices
  data_external/             External reference resources
  metadata/                  Curated dataset annotations, sample-label audits, QC tables
  scripts/                   Numbered analysis and manuscript-build scripts
  results/                   Generated analysis outputs
  manuscript/                Submission DOCX/PDF, figures, editable figure sources, additional files, and local deposit package
  references/                Reference notes and source-tracking materials
  logs/                      Local run logs
```

## Reproducibility entry points

For journal review, the most important audit files are:

- Public repository: https://github.com/osbornzhou/IPF_oligo.
- Additional files 1-14 are provided through the journal submission system; the repository provides the scripts and reproducibility metadata needed to audit or regenerate those outputs from local source data.
- `manuscript/reproducibility_package/script_input_output_map.csv`: current numbered script order and input/output map.
- `manuscript/reproducibility_package/software_environment.csv`: software and package notes.
- `manuscript/reproducibility_package/random_seeds.csv`: random seeds used in model and sensitivity analyses.
- `manuscript/reproducibility_package/file_manifest.csv`: checksums for key generated artifacts.
- `manuscript/editable_figures/`: editable SVG sources for main Figures 1-6; these are vector/text source files for post-generation figure alignment edits.
- The local submission workspace includes a repository deposit archive for upload or release preparation; the public GitHub repository itself provides the current scripts and reproducibility metadata.

The numbered scripts are intended to be read together with `script_input_output_map.csv`. Some scripts assume that public GEO source files and processed intermediates already exist locally under the paths recorded in the metadata and result manifests.

## Interpretation boundary

The perturbation-triage score is a structured validation-planning heuristic. It is not a learned model of oligonucleotide efficacy, not evidence of therapeutic readiness, and not a substitute for qPCR, protein-level, siRNA/ASO, miRNA mimic, or animal-model validation.
