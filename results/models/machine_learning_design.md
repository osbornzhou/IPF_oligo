# Machine Learning Module Design for IPF Oligonucleotide-Target Project

## Overall Goal

Build a rigorous, multi-layer machine learning module that supports the manuscript narrative:

1. Identify a compact and stable IPF diagnostic mRNA biomarker panel.
2. Compare classical machine learning models with a small neural-network model.
3. Prioritize oligonucleotide-targetable miRNA-mRNA regulatory modules using model evidence plus biological evidence.
4. Validate model performance across independent GEO cohorts and distinct platforms.

This module should not be a simple single-dataset classification exercise. It should be a discovery-validation framework with leakage control, cross-platform robustness, model stability analysis, and biological interpretability.

## Literature-Informed Rationale

The design follows three ideas from BMC Genomics-style studies:

- Biomarker studies should reduce high-dimensional omics data into compact and biologically meaningful panels, rather than reporting hundreds of genes.
- Feature selection should be stable across resampling, because high-dimensional small-sample omics data can produce unstable signatures.
- Neural networks should be compared with classical models rather than presented alone; omics sample sizes are usually too small for deep architectures to be trusted without strong validation.

The provided BMC Genomics PDF by Lin and Qiu also suggests a useful principle for this project: use multi-evidence consistency. Their study used correlation-pattern features, curated miRNA-target labels, XGBoost, held-out validation, and perturbation-direction checks. Here, the analogous strategy is:

- use robust discovery-validation mRNA and miRNA candidates as model inputs;
- use external GEO datasets as held-out validation;
- use miRNA-mRNA negative regulation, PPI hubness, and enrichment evidence as biological consistency checks.

## Data Partitioning

### Training / Discovery

- Main training cohort: `GSE32537`
- Task: IPF vs Control
- Excluded samples: all samples marked `Exclude`
- Candidate feature universe:
  - Primary: `robust_mrna_candidates_strict.csv`
  - Secondary sensitivity sets:
    - PPI top hub genes
    - miRNA-mRNA axis target genes
    - robust mRNA up/down subsets

### External Validation

External validation cohorts:

- `GSE110147`
- `GSE150910`
- `GSE53845`
- `GSE92592`

The validation datasets must never be used for model fitting, feature selection, hyperparameter tuning, or threshold optimization.

## Expression Matrix Harmonization

For each mRNA dataset:

1. Convert probe/transcript IDs to `standard_feature_id`.
2. Collapse multiple probes per gene by keeping the probe with highest variance within that dataset.
3. Keep only genes available in the model feature panel.
4. Apply dataset-wise transformation:
   - count-like datasets: `log2(CPM + 1)` or already voom-compatible log-like matrix;
   - array/normalized datasets: use existing normalized scale.
5. Standardize each gene within each dataset using z-score.

Important leakage rule:

- For internal cross-validation inside `GSE32537`, scaling parameters must be estimated inside each training fold and applied to the held-out fold.
- For external validation, scaling can be done within each external dataset using all external samples because no labels are used. However, label-informed operations are forbidden.

## Feature Selection Strategy

Use a nested feature-selection design.

### Stage 1: Biological Pre-filter

Allowed input feature spaces:

1. `robust_mrna_strict`: 280 robust mRNA candidates.
2. `hub_augmented`: robust mRNA plus top PPI hub evidence.
3. `axis_target_focused`: genes from miRNA-mRNA negative axes.
4. `combined_evidence`: genes ranked by combined evidence:
   - discovery FDR
   - number of validation datasets with same-direction support
   - STRING hub score
   - presence in miRNA-mRNA negative axis
   - enrichment/pathway relevance

Primary model should use `robust_mrna_strict`. The other sets are sensitivity analyses.

### Stage 2: Statistical / ML Feature Selection

Inside training folds only:

- LASSO logistic regression
- Elastic Net logistic regression
- SVM-RFE
- Random Forest importance
- XGBoost importance
- Mutual information or univariate logistic filter as a lightweight baseline

Final candidate panel should be selected by stability:

- gene selected in at least 50% of repeated CV runs; or
- top-ranked by at least 3 independent selection methods; or
- included in a compact model with strong external validation AUC.

## Models to Compare

### Baseline Models

1. Logistic regression using all selected features.
2. Penalized logistic regression:
   - LASSO
   - Elastic Net

### Classical ML Models

3. Support Vector Machine:
   - linear SVM
   - RBF SVM as sensitivity
4. Random Forest
5. XGBoost / gradient boosting
6. k-nearest neighbors as a simple non-linear baseline

### Neural Network Model

7. Small MLP neural network:
   - input: selected gene panel only, not thousands of genes
   - hidden layers: 1-2 small layers, e.g. 16 and 8 nodes
   - dropout: 0.2-0.4
   - L2 regularization
   - early stopping

Do not use a large deep model. The sample size is too small and would be easy to overfit.

## Cross-Validation Design

Inside `GSE32537`:

- repeated stratified 5-fold CV;
- repeat 50-100 times for stability;
- nested hyperparameter tuning inside training folds;
- metrics averaged across repeats.

External validation:

- train final model on full `GSE32537`;
- test once on each validation dataset;
- report per-dataset AUC and pooled external AUC;
- do not retune on validation datasets.

## Metrics

Classification metrics:

- ROC AUC
- PR AUC
- accuracy
- sensitivity
- specificity
- balanced accuracy
- F1 score
- calibration slope/intercept if feasible

Feature stability metrics:

- selection frequency
- mean rank across repeats
- Jaccard similarity of selected panels
- coefficient/importance direction consistency

Biological evidence metrics:

- same-direction validation count
- PPI hub score
- miRNA-mRNA axis membership
- enrichment/pathway membership
- oligonucleotide targetability class

## Final Model Selection

Choose the final model by a hierarchy:

1. External validation AUC is strong and consistent.
2. Feature panel is small enough for biological interpretation.
3. Selected genes are stable across resampling.
4. Model is interpretable enough for manuscript presentation.
5. Key genes connect to miRNA-mRNA axes, PPI hubs, or IPF-relevant pathways.

Preferred final model if performance is similar:

- Elastic Net logistic regression or linear SVM.

Reason:

- easier to explain than black-box models;
- more suitable for BMC Genomics and biomedical reviewers;
- provides coefficients and direction.

Use Random Forest, XGBoost, and MLP as comparative/confirmatory models.

## Oligonucleotide-Target Prioritization Layer

After selecting the diagnostic panel, assign each candidate gene or miRNA-mRNA axis a target priority score.

Potential score components:

- ML selection frequency
- model coefficient or SHAP-like importance
- robust differential expression evidence
- external validation support count
- miRNA-mRNA negative axis evidence
- STRING hub score
- pathway relevance
- direction of therapeutic intervention:
  - upregulated mRNA: ASO/siRNA knockdown candidate
  - downregulated miRNA: miRNA mimic candidate
  - upregulated miRNA: antagomir/ASO candidate

High-priority examples to evaluate:

- `hsa-miR-375 -> CLDN1`
- `hsa-miR-375 -> MNS1`
- `hsa-miR-92a -> TP63`
- `hsa-miR-30a -> CDH2`
- `hsa-miR-92a -> GOLM1`
- `hsa-miR-92a -> LTBP1`

## Triple QC for Machine Learning

### QC1: Data and Feature Integrity QC

Purpose: prevent wrong input and sample/feature mismatch.

Checks:

- sample labels are only IPF or Control;
- no `Exclude` samples enter model training;
- sample IDs match expression matrix columns;
- no duplicated sample IDs;
- feature IDs are standardized gene symbols;
- no missing or infinite expression values after harmonization;
- train and validation cohorts share required feature columns;
- class balance is recorded for every cohort;
- no validation dataset is used during feature selection.

Pass criteria:

- all model matrices have valid labels and numeric features;
- feature overlap is sufficient for each planned panel;
- no sample leakage between train and validation.

### QC2: Resampling, Tuning, and Leakage QC

Purpose: ensure model estimates are honest.

Checks:

- stratified folds preserve IPF/Control ratio;
- scaling/imputation is fit only on training folds;
- feature selection is nested inside training folds;
- hyperparameter tuning is nested inside training folds;
- random seeds are fixed and recorded;
- each model has a complete parameter log;
- external validation is untouched until final evaluation.

Pass criteria:

- all folds complete successfully;
- no fold has one class only;
- no validation metrics are used for tuning;
- final model can be exactly reproduced.

### QC3: Performance, Stability, and Biological Plausibility QC

Purpose: avoid overfitted or biologically implausible signatures.

Checks:

- internal CV AUC and external AUC are both reported;
- performance drop from internal to external validation is quantified;
- sensitivity and specificity are both acceptable;
- selected features have stable selection frequency;
- selected genes show direction consistency across validation cohorts;
- final genes overlap with robust DE, PPI hub, pathway, or miRNA-mRNA evidence;
- permutation-label control performs near chance;
- optional: random-gene-set control does not outperform final model.

Pass criteria:

- final model has acceptable external AUC in multiple validation cohorts;
- no severe internal-external performance collapse;
- top features are stable and biologically interpretable;
- permutation control confirms the model is not learning noise.

## Planned Outputs

Tables:

- `ml_dataset_qc.csv`
- `ml_feature_matrix_qc.csv`
- `ml_model_performance_internal_cv.csv`
- `ml_model_performance_external_validation.csv`
- `ml_feature_selection_stability.csv`
- `ml_final_biomarker_panel.csv`
- `ml_oligonucleotide_target_priority.csv`
- `ml_triple_qc.csv`

Figures:

- ROC curves for internal CV and each external validation dataset
- PR curves
- model comparison barplot
- feature selection stability plot
- final biomarker heatmap
- final model coefficient/importance plot
- oligonucleotide target prioritization plot

## Manuscript Position

This module should appear after differential expression, miRNA-mRNA axis, enrichment, and PPI analyses.

Suggested Results subsection title:

Machine learning identifies a compact and externally validated IPF biomarker panel

Suggested Methods subsection title:

Machine learning model development, external validation, and feature stability analysis

