# IPF multi-cohort transcriptomic actionability framework

This repository contains the reproducibility package for the manuscript:

**A multi-cohort transcriptomic actionability framework prioritizes IPF-associated mRNA and miRNA-mRNA axes for oligonucleotide-focused validation**

The project analyzes public human idiopathic pulmonary fibrosis (IPF) bulk mRNA, miRNA, and single-cell transcriptomic datasets to build a QC-traceable computational resource for candidate prioritization. The workflow is designed for validation planning, not for therapeutic-candidate validation or clinical diagnostic deployment.

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
- oligonucleotide validation-planning actionability scoring and score-sensitivity checks;
- machine-readable additional files and manuscript figures.

## Repository layout

```text
ipf_oligo_ml/
  data_raw/                  Public GEO source files where available locally
  data_processed/            Intermediate processed matrices
  data_external/             External reference resources
  metadata/                  Curated dataset annotations, sample-label audits, QC tables
  scripts/                   Numbered analysis and manuscript-build scripts
  results/                   Generated analysis outputs
  manuscript/                Submission DOCX/PDF, figures, additional files, and deposit package
  references/                Reference notes and source-tracking materials
  logs/                      Local run logs
```

## Reproducibility entry points

For journal review, the most important audit files are:

- `manuscript/additional_files/`: Additional files 1-14 referenced by the manuscript.
- `manuscript/reproducibility_package/script_input_output_map.csv`: current numbered script order and input/output map.
- `manuscript/reproducibility_package/software_environment.csv`: software and package notes.
- `manuscript/reproducibility_package/random_seeds.csv`: random seeds used in model and sensitivity analyses.
- `manuscript/reproducibility_package/file_manifest.csv`: checksums for key generated artifacts.
- `manuscript/repository_deposit_package.zip`: ready-to-upload repository archive containing scripts, additional files, reproducibility metadata, figures, and manuscript artifacts.

The numbered scripts are intended to be read together with `script_input_output_map.csv`. Some scripts assume that public GEO source files and processed intermediates already exist locally under the paths recorded in the metadata and result manifests.

## Interpretation boundary

The actionability score is a structured validation-planning heuristic. It is not a learned model of oligonucleotide efficacy, not evidence of therapeutic readiness, and not a substitute for qPCR, protein-level, siRNA/ASO, miRNA mimic, or animal-model validation.
