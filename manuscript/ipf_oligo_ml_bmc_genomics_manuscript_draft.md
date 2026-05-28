# A boundary-tested transcriptomic perturbation-triage framework identifies oligonucleotide validation candidates and fibrotic disease-state markers in idiopathic pulmonary fibrosis

Authors: Yunyi Zhou¹; Yanli Zhang¹*

Affiliation: ¹State Key Laboratory of Common Mechanism Research for Major Diseases, Department of Biochemistry and Molecular Biology, Institute of Basic Medical Sciences, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, China.

Corresponding author: Yanli Zhang, zhangyanli@ibms.pumc.edu.cn

## Abstract
Background: Idiopathic pulmonary fibrosis (IPF) is a progressive fibrotic lung disease with limited therapeutic options. Public transcriptomic datasets identify many reproducible IPF-associated signals, but only a subset are suitable for oligonucleotide-focused perturbation validation. We asked which reproducible IPF/fibrotic transcriptomic abnormalities should be advanced as perturbation-screening candidates and which should instead be triaged as disease-state markers, restoration hypotheses, or low-confidence regulatory hypotheses.

Results: Eight bulk mRNA or miRNA GEO datasets passed three-layer quality control, comprising 691 profiled samples. Cross-cohort screening retained 280 robust mRNA candidates and 10 robust miRNA candidates, defining a reproducible IPF/fibrotic transcriptomic evidence base rather than a stand-alone biomarker list. A leakage-controlled Elastic Net score trained with discovery-only mRNA features showed strong external disease-state separation, but boundary stress tests using excluded GSE110147 NSIP and mixed IPF-NSIP samples, matched random-panel baselines, non-perfect-cohort summaries, and cohort-adjusted modeling indicated that the score captured a fibrotic interstitial-lung-disease state rather than IPF-specific diagnostic specificity. miRTarBase integration yielded 22 opposite-direction miRNA-mRNA candidate axes; strict evidence gates retained 3 exact mature hsa-miR-375 axes for main-text interpretation and downgraded 19 arm-agnostic axes to exploratory status. Target-program stress testing did not support expansion to a broader miRNA derepression program. Donor-aware single-cell pseudobulk validation localized key candidate signals, including SPP1 in myeloid cells, COL1A1/COL3A1/PTGFRN in stromal cells, CD24 in epithelial cells, and reduced GPX3 in stromal cells. The final perturbation-triage framework separated knockdown-screening candidates, context-dependent fibrotic disease-state markers, restoration/pathway markers, exact miRNA-axis hypotheses, and an externally motivated TNIK bridge.

Conclusions: This multi-cohort transcriptomic study provides a boundary-tested perturbation-triage framework for organizing IPF/fibrotic transcriptomic abnormalities before oligonucleotide-focused validation. The findings support epithelial, stromal, immune, and ciliary remodeling programs in IPF while emphasizing perturbation triage rather than diagnostic deployment or therapeutic-candidate validation.

## Keywords
Idiopathic pulmonary fibrosis; oligonucleotide therapeutics; machine learning; miRNA; transcriptomics; single-cell RNA sequencing; candidate prioritization

## Background
Idiopathic pulmonary fibrosis (IPF) is characterized by progressive distortion of lung architecture, aberrant epithelial repair, fibroblast activation, extracellular matrix accumulation, and immune remodeling [23,24]. Although antifibrotic therapies can slow functional decline, they do not reverse established disease, and there remains a need for mechanistically grounded molecular targets [23,24]. Oligonucleotide therapeutics provide a direct route to modulate transcripts or miRNA activity, but the translational value of a candidate depends on reproducibility across cohorts, compatibility with target directionality, and evidence that the candidate is active in relevant lung cell populations [27,28].

Public transcriptomic resources make it possible to evaluate these criteria systematically [1,2]. However, many IPF biomarker studies rely on single discovery cohorts, limited validation, or feature selection strategies that risk information leakage when validation data are used before model testing [29-34]. In addition, mRNA and miRNA evidence is often analyzed separately, making it difficult to nominate coordinated regulatory axes for oligonucleotide intervention [20,21].

Here, we designed a multi-cohort transcriptomic perturbation-triage study that combines GEO mRNA and miRNA datasets, external validation, miRTarBase target evidence, pathway enrichment, STRING protein interaction analysis, machine learning, disease-boundary stress testing, and single-cell localization [1,2,17-22,25,26]. The primary objective was not to propose another IPF biomarker classifier, but to determine which reproducible IPF/fibrotic transcriptomic signals are suitable for oligonucleotide-focused perturbation validation and which should be downgraded as disease-state markers, restoration hypotheses, or low-confidence regulatory hypotheses.

## Results
### QC-controlled multi-cohort dataset and sample-label audit
The workflow integrated bulk mRNA, miRNA, and single-cell transcriptomic data from public GEO datasets (Figure 1) [1,2,25,26]. Eight bulk or miRNA datasets passed three independent quality-control checks: annotation completeness, expression-sample cross-matching, and matrix integrity. These datasets contained 495 bulk mRNA samples (304 IPF and 191 controls) and 196 miRNA samples (128 IPF and 68 controls). Two single-cell datasets were available for expression-level validation, while GSE122960 was excluded because a cell-level expression matrix was not available locally.

Table 1 summarizes the quality-controlled expression datasets.

| series_id | data_type | dataset_role | samples | IPF/control | matrix_type | DE_method | QC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE110147 | bulk mRNA | validation | 33 | 22/11 | microarray normalized log-intensity | limma | pass |
| GSE150910 | bulk mRNA | validation | 206 | 103/103 | supplementary gene-level raw counts | edgeR-voom | pass |
| GSE21394 | miRNA | validation | 15 | 9/6 | miRNA microarray normalized log-intensity | limma | pass |
| GSE27430 | miRNA | validation | 25 | 13/12 | centered/transformed miRNA expression | limma | pass |
| GSE32537 | bulk mRNA | discovery | 169 | 119/50 | microarray normalized log-intensity | limma | pass |
| GSE32538 | miRNA | discovery | 156 | 106/50 | centered/transformed miRNA expression | limma | pass |
| GSE53845 | bulk mRNA | validation | 48 | 40/8 | centered/transformed microarray expression | limma | pass |
| GSE92592 | bulk mRNA | validation | 39 | 20/19 | supplementary mapped gene-count matrix | edgeR-voom | pass |

Dataset source publications were cited where available for the major bulk, miRNA, and single-cell cohorts [6,25,26,36-41]. GEO accession records were retained in the reference list to provide persistent dataset-level identifiers.

GSE110147 was treated as an IPF-vs-normal validation subset from a study that also included NSIP and mixed IPF-NSIP samples. The main analysis included 22 IPF samples and 11 normal controls; 10 NSIP and 5 mixed IPF-NSIP samples were excluded before differential expression and machine-learning validation. Disease labels were assigned from GEO sample metadata and manually checked against sample titles and source descriptions; all included and excluded GSE110147 samples are listed in the sample-level label audit in Additional file 1.

Sample-level inclusion, exclusion, and disease-label audit tables for all bulk and miRNA datasets are provided in Additional file 1, including sample accession, original curated disease label, final analysis label, inclusion status, exclusion reason, and source annotation fields.

### Cross-cohort screening defines a reproducible IPF/fibrotic transcriptomic evidence base
Differential expression was performed using limma for log-scale or normalized matrices and edgeR-limma voom for count-like matrices. In the mRNA discovery dataset GSE32537, 483 features passed FDR < 0.05 and absolute log fold change >= 1; after feature annotation and gene-level harmonization, 449 mRNA genes were used for discovery-validation screening. In the miRNA discovery dataset GSE32538, 47 significant miRNAs were detected.

Direction-consistency screening across independent validation datasets retained 280 strict robust mRNA candidates and 10 strict robust miRNA candidates. The highest-ranked robust mRNAs included ASPN, COL14A1, NECAB1, CDH3, LRRC17, CD24, CFH, DIO2, DCLK1, CXCL14, while the robust miRNA set included hsa-miR-375, hsa-miR-30d, hsa-miR-30a, hsa-miR-423-5p, hsa-miR-31, hsa-miR-205, hsa-miR-34c-5p, hsa-miR-92a, hsa-miR-203, hsa-miR-141. These candidates were treated as a reproducible IPF/fibrotic transcriptomic evidence base for downstream perturbation triage rather than as a final biomarker list. They were carried forward for candidate-axis construction, enrichment, network analysis, disease-state scoring, single-cell localization, and perturbation-triage classification.

Complete dataset-level differential-expression summaries, annotated feature-level results, and gene- or miRNA-level collapsed outputs are provided in Additional file 2.

**Supporting enrichment and protein-interaction analyses.**
GO, KEGG, and Reactome enrichment analyses used the gene-level feature universe tested in the GSE32537 discovery analysis as background [17,18]. For the 280 strict robust mRNAs, 251 symbols mapped successfully and 81 significant enrichment terms were detected. The most significant terms included cilium movement, cilium movement involved in cell motility, microtubule-based movement, and axoneme assembly. Terms containing sperm or flagellar labels were interpreted as shared axonemal/ciliary structural programs rather than reproductive biology.

STRING analysis mapped 267 of 280 robust mRNAs at medium confidence, yielding 950 edges and a largest connected component of 219 nodes [19]. The top-ranked hubs were DNAI1, TEKT1, CFAP52, COL1A1, COL3A1, POSTN, RSPH4A, SPAG17, TTC25, COL6A3, COL1A2, VCAM1.

Complete GO, KEGG, Reactome, STRING mapping, edge, and hub metric outputs are provided in Additional files 5 and 6.

### Disease-boundary testing supports a fibrotic ILD state score rather than IPF-specific diagnostic specificity
To avoid validation leakage, machine learning used only mRNA features selected from the GSE32537 discovery analysis. After harmonizing gene-level matrices across the discovery and external validation cohorts, 329 common discovery-only mRNA features were available. Seven model families were evaluated with imputation, scaling, and SelectKBest feature selection inside scikit-learn pipelines [22]; feature selection and tuning were nested within resampling folds. Elastic Net achieved the highest mean external ROC AUC. Across GSE110147, GSE150910, GSE53845, and GSE92592, the best model had mean external ROC AUC 0.971, minimum external ROC AUC 0.910, mean external PR AUC 0.964, and mean balanced accuracy 0.942. These headline values are arithmetic means of per-cohort metrics unless otherwise stated. A label-permutation control produced a mean AUC of 0.531, supporting non-random predictive signal under the modeling workflow while not excluding cohort-composition or platform effects.

The final stable panel contained 25 genes: ASPN, PTGFRN, COL14A1, GPX3, FNDC1, NECAB1, CDH3, CFH, HMCN1, CD24, LRRC17, KLHL13, NPR1, CSF3R, SYTL2, PSD3, SLCO4A1, PDE7B, GRASP, CBS, MGAM, IL1R2, GLT1D1, TNFRSF19, TPPP3. Many top-ranked panel genes were also robustly direction-consistent across validation datasets and were annotated with cautious follow-up labels, such as conceptual siRNA or antisense knockdown compatibility for upregulated mRNAs. These labels were used to organize validation priorities and should not be read as evidence of therapeutic readiness. Pooled external interpretability analysis gave a ROC AUC of 0.924 across 326 validation samples, with GPX3 as the largest permutation-importance feature followed by GRASP, CSF3R, CFH, ID1, FASN. Calibration and decision-curve outputs were generated to support model-review discussions rather than to claim immediate clinical deployability.

Two additional sensitivity checks addressed whether the high external performance was dominated by a single validation cohort or by generic robust-gene signal. Per-cohort Elastic Net ROC AUC values were 1.000 in GSE110147, 0.910 in GSE150910, 1.000 in GSE53845, and 0.974 in GSE92592. When each validation cohort was omitted from the summary, the mean cohort ROC AUC across the remaining cohorts ranged from 0.961 to 0.991 and the minimum retained cohort ROC AUC remained at least 0.910. In 500 random 25-gene panels sampled from the common robust mRNA universe and refit with the same fixed Elastic Net classifier, the mean external ROC AUC was 0.889 +/- 0.066, with a 95th percentile of 0.971. The observed 25-gene panel refit achieved mean external ROC AUC 0.964 and minimum external ROC AUC 0.907. These results support enrichment of predictive signal in the selected panel. At the same time, the random-panel analysis indicates that robust IPF-associated genes broadly carry disease-state information; therefore, the model is best interpreted as a cross-cohort disease-state prioritization tool rather than a clinically deployable classifier.

Because the external ROC AUC values were high, we added disease-state boundary stress tests rather than additional model families (Additional file 8; Figure 3). Excluded GSE110147 NSIP and mixed IPF-NSIP samples were scored without using them for model training. Their scores were high and close to IPF samples (median score 0.930 for NSIP, 0.996 for mixed IPF-NSIP, 0.971 for IPF, and 0.095 for normal controls), indicating that the model captured a fibrotic interstitial-lung-disease state rather than IPF-specific diagnostic specificity. When validation cohorts with perfect ROC AUC were removed, the retained GSE150910 and GSE92592 cohorts still showed mean ROC AUC 0.942, mean PR AUC 0.929, and mean balanced accuracy 0.932. A cohort-adjusted external logistic model retained an association between the Elastic Net logit score and disease status after validation-cohort fixed effects (coefficient 0.218; Wald P = 6.39e-24). In a matched random discovery-feature baseline, the observed final-panel refit had mean ROC AUC 0.946, close to the matched random-panel mean of 0.944 and below its 95th percentile of 0.968, while balanced accuracy remained at the high end of the matched distribution. These stress tests support use of the score as a disease-state prioritization layer, not as a uniquely optimized IPF diagnostic signature.

Table 3 summarizes the external performance and sensitivity checks. Per-cohort bootstrap 95% confidence intervals for ROC AUC and PR AUC, the full 25-gene panel, feature stability, predictions, hyperparameter grids, final Elastic Net specification, random-panel outputs, and disease-state boundary stress tests are provided in Additional files 7 and 8.

| analysis | metric | value |
| --- | --- | --- |
| Final Elastic Net | mean external ROC AUC | 0.971 |
| Final Elastic Net | minimum external ROC AUC | 0.910 |
| Final Elastic Net | mean external PR AUC | 0.964 |
| Final Elastic Net | mean balanced accuracy | 0.942 |
| Label permutation | mean ROC AUC | 0.531 |
| Leave-one-validation-cohort | retained mean ROC AUC range | 0.961-0.991 |
| Random robust 25-gene panels | mean external ROC AUC +/- SD | 0.889 +/- 0.066 |
| Observed 25-gene refit | mean external ROC AUC | 0.964 |

**Supporting reported-gene-set benchmarking.**
Because most published IPF machine-learning studies do not provide deployable model objects or full coefficients, we evaluated reported gene signatures as fixed external comparator feature sets [29-34]. Each comparator was trained only in GSE32537 and externally validated under the same framework used for the proposed model. The best published comparator, an explainable machine-learning reported-gene set, achieved mean external ROC AUC 0.946 [32]. Under this harmonized reported-gene-set benchmarking design, the proposed discovery-only Elastic Net panel showed higher mean external ROC AUC (0.971) and a higher minimum external ROC AUC (0.910), suggesting improved cross-cohort stability within this specific benchmarking framework. This comparison should be interpreted as reported gene-set benchmarking rather than direct comparison with fully deployable published models.

Complete comparator metrics are provided in Additional file 7.

### Strict miRNA evidence gates retain exact hsa-miR-375 axes and downgrade weak regulatory hypotheses
We integrated the robust miRNA and mRNA candidate sets with the human miRTarBase 2025 v10 interaction table [20,21]. Of 10 robust miRNAs and 280 robust mRNAs, 7936 miRTarBase records matched candidate miRNAs. Applying an opposite-direction rule between miRNA and mRNA discovery log fold changes identified 22 negative candidate axes involving 5 miRNAs and 16 target genes. Evidence grading separated 3 exact mature-miRNA matches from 19 arm-agnostic matches. hsa-miR-375 was downregulated in discovery (logFC -1.52) and showed same-direction decreases in both miRNA validation datasets (GSE21394 logFC -0.72; GSE27430 logFC -1.10). The exact hsa-miR-375 -> CLDN1, hsa-miR-375 -> MNS1, and hsa-miR-375 -> RPGRIP1L axes were retained as main-text prioritized hypotheses. The focused support table for these axes lists miRTarBase IDs, support type, validation method, PMID, miRNA/target species, and targeted cross-check status against accessible TargetScan/literature, miRDB, miRWalk, and CLIP-supported evidence. CLDN1 had additional literature support as a TargetScan-nominated and experimentally tested miR-375 target in lung cancer cells [35], whereas MNS1 and RPGRIP1L remained miRTarBase-derived exact-axis hypotheses without independently verified prediction or CLIP support in the current version. Arm-agnostic matches, including hsa-miR-92a with TP63 and hsa-miR-30a with CDH2, were retained only as exploratory supplementary hypotheses because they require mature-arm-specific validation before therapeutic interpretation.

The exact mature-miRNA axes prioritized in the main text are shown in Table 2. Arm-agnostic axes are retained as exploratory hypotheses in Additional file 4.

| axis | match_type | axis_score |
| --- | --- | --- |
| hsa-miR-375 -> CLDN1 | exact | 49.000 |
| hsa-miR-375 -> MNS1 | exact | 37.000 |
| hsa-miR-375 -> RPGRIP1L | exact | 37.000 |

Because exact mature-miRNA axes were intentionally stringent, we also tested whether robust downregulated miRNAs showed broader target-program evidence among IPF-upregulated robust mRNAs (Additional file 4; Figure 4). miRTarBase target sets for robust downregulated miRNAs did not show FDR-supported enrichment in the IPF-upregulated robust mRNA foreground. For hsa-miR-375, 250 miRTarBase exact or arm-recoverable targets were present in the discovery mRNA background, but only 3 overlapped the upregulated robust mRNA foreground (CLDN1, MNS1, and RPGRIP1L; odds ratio 0.63; FDR = 1.0). hsa-miR-375 targets also did not show a release-like increase in discovery mRNA log fold change relative to non-target genes (one-sided permutation P = 0.932). Therefore, the miRNA layer was not expanded beyond the three exact hsa-miR-375 axes in the main interpretation; arm-agnostic and relaxed axes remain exploratory.

### Donor-aware single-cell pseudobulk localizes candidate signals
Single-cell validation was performed in GSE135893 and GSE136831, which together provided 332798 IPF/control metadata-matched cells [25,26]. The pipeline streamed sparse matrices directly, validated metadata alignment, and filtered non-informative cell labels including Multiplet, Outlier, MT-tRNAs, and CellCycle categories. Of 59 requested candidate genes, 59 were found in GSE135893 and 58 were found in GSE136831.

The most prominent cell-type signals localized SPP1 to macrophage or myeloid compartments, COL1A1, COL1A2, COL3A1, and POSTN to stromal myofibroblast compartments, and GPX3 downregulation to fibroblast or myofibroblast populations. These cell-level summaries provide localization context rather than donor-level causal inference. They support a disease model in which epithelial/ciliary changes, stromal matrix remodeling, and macrophage activation converge in IPF and provide cellular context for oligonucleotide-focused follow-up.

To reduce pseudoreplication concerns, we added a donor-aware pseudobulk sensitivity analysis for core candidates (Figure 5C). Counts were aggregated at the donor/sample x broad-celltype level before IPF-control comparison. Both single-cell datasets passed donor metadata and matrix-integrity QC, with 22 donors in GSE135893 and 60 donors in GSE136831. The pseudobulk analysis supported SPP1 upregulation in myeloid cells, COL1A1/COL3A1 and PTGFRN upregulation in stromal cells, POSTN upregulation in endothelial/stromal contexts, CD24 upregulation in epithelial cells, and GPX3 downregulation in stromal cells in GSE136831; COL14A1 was supported in mesenchymal cells in GSE135893. POSTN pseudobulk support in the endothelial-context comparison was interpreted as dataset-specific localization support, while POSTN was retained overall as a context-dependent matrix/pathway candidate rather than a simple knockdown-screening candidate. All displayed core pseudobulk comparisons were direction-consistent with the bulk evidence and passed BH-adjusted FDR < 0.01, with the largest donor-aware effects observed for SPP1 in myeloid cells, POSTN in endothelial cells, COL1A1/COL3A1 in stromal cells, and reduced GPX3 in stromal cells. These donor-level summaries strengthen the localization layer while remaining validation-oriented rather than a full single-cell differential-expression study.

Complete single-cell localization summaries, donor-aware pseudobulk outputs, cell-count summaries, and candidate-gene context tables are provided in Additional file 9.

**Supporting immune-stromal module and single-cell communication analyses.**
To avoid treating the final candidates as isolated biomarkers, we added three exploratory mechanistic layers after the main discovery and validation analyses. First, coexpression modules were used to test whether robust and machine-learning candidates clustered into coordinated IPF-associated transcriptional programs. Second, curated ligand-receptor scoring was used to place prioritized genes in an intercellular signaling context. Third, a coexpression-neighborhood perturbation-priority proxy was used to organize knockdown-oriented follow-up candidates. These analyses were designed to support biological interpretation and experimental prioritization rather than to establish causal mechanisms.

The module analysis addressed whether the 25-gene panel and robust candidates were scattered classifier features or components of broader disease programs. In the GSE32537 discovery cohort, ten WGCNA-like coexpression modules were detected. M02 contained 211 robust mRNA candidates and 15 machine-learning panel genes, correlated with IPF status (r = 0.66), and aligned with epithelial/ciliary and matrix-remodeling marker scores. M09 contained 68 robust mRNA candidates, 10 machine-learning panel genes, and TNIK; it showed the strongest IPF correlation (r = 0.74). Thus, the module analysis provided a bridge between feature prioritization and pathway-level disease biology, while remaining an exploratory correlation analysis rather than a causal network model.

The curated ligand-receptor score was calculated as the product of mean normalized ligand expression in a sender cell type and mean normalized receptor expression in a receiver cell type, followed by the IPF-control difference for the corresponding sender-receiver pair. A positive IPF-control score therefore indicates higher expression-based communication potential in IPF than in controls, not a CellChat probability, statistical interaction test, or experimentally validated signaling flux. After excluding non-informative cell labels such as multiplets and outliers, 934 interaction changes were scored. The strongest IPF-increased expression proxies involved MIF-CD74 across endothelial, epithelial, and immune compartments; COL14A1-ITGB1 within mesenchymal cells; CXCL12-CXCR4 from mesenchymal to immune cells; POSTN-ITGB1; and SPP1-integrin/CD44 signaling. These results place the candidates in an immune-stromal communication context but should be interpreted as hypothesis-generating expression proxies.

Finally, the coexpression-neighborhood perturbation-priority proxy was not intended to simulate molecular knockdown. Instead, it ranked whether suppressing an upregulated candidate would be expected to align with reversal of its local IPF-associated coexpression neighborhood in the discovery cohort. This proxy prioritized CD24, COL14A1, POSTN, and PTGFRN as internally supported perturbation-screening candidates, but POSTN was retained as a context-dependent matrix/pathway candidate because it is a broad fibrotic matrix marker rather than a simple reversal target. SPP1 remained biologically important through macrophage/myeloid localization and ligand-receptor evidence, but its negative proxy score indicates that it was not as well supported as a simple bulk coexpression-neighborhood reversal target. TNIK had validation-cohort differential-expression and Wnt/TNIK bridge evidence, but it was not significant in the GSE32537 discovery analysis and was not selected by the final machine-learning panel. Therefore, TNIK is best positioned as an externally supplied translational hypothesis to test alongside the primary IPF target modules rather than as a target discovered de novo by this study.

Complete coexpression module, curated ligand-receptor, coexpression-neighborhood perturbation-priority, and TNIK evidence tables are provided in Additional file 10. The main-text interpretation is intentionally limited to whether these outputs support the primary candidate modules or motivate cautious external hypotheses.

### Perturbation-triage separates reversal-screening candidates from disease-state markers and restoration hypotheses
To convert the multi-omic results into a practical shortlist for oligonucleotide-focused follow-up, we built an evidence-weighted prioritization table that combined robust differential expression, machine-learning panel membership, PPI hub status, pathway membership, miRNA-axis evidence, single-cell localization, and direction-compatible oligonucleotide strategy labels. The priority score was designed to organize follow-up experiments, not to rank therapeutic promise directly. This distinction is important because upregulated candidates such as COL14A1, PTGFRN, CD24, ASPN, and CDH3 can be considered for knockdown screening, whereas broad matrix or inflammatory markers such as POSTN, COL1A1, COL3A1, and SPP1 require context-dependent pathway screening and downregulated candidates such as GPX3 and NECAB1 are better interpreted as pathway-state markers or restoration hypotheses. In total, 288 robust genes were scored, of which 10 were assigned Tier 1 priority. The highest-ranked genes were COL14A1, GPX3, ASPN, COL1A1, NECAB1, SPP1, PTGFRN, CD24, POSTN, CDH3.

To further separate this study from biomarker-only discovery, we converted the integrated evidence into an oligonucleotide-focused perturbation-triage map (Figure 6). This triage score did not add a new disease-association test; instead, it translated the existing evidence layers into experimental validation classes. Knockdown-screening candidates included upregulated genes with cross-cohort support and a plausible knockdown-screening route, such as CD24, COL14A1, PTGFRN, ASPN, and CDH3. Context-dependent candidates, including POSTN, SPP1, COL1A1, and COL3A1, were retained as biologically important but not necessarily simple reversal candidates. GPX3 and NECAB1 were assigned to restoration or pathway-marker classes because they were downregulated in IPF. The exact hsa-miR-375 target genes CLDN1, MNS1, and RPGRIP1L were grouped as miRNA-axis hypotheses, and TNIK was retained as an external Wnt/TNIK bridge rather than a primary discovery candidate. The perturbation-triage score should be interpreted as a structured validation-planning heuristic, not a validated predictive model of oligonucleotide efficacy.

Integrated priority score and perturbation-triage score serve different purposes. The integrated priority score ranks disease-association evidence across robust expression, machine-learning, PPI, enrichment, miRNA-axis, and single-cell layers. The perturbation-triage score reorganizes a subset of candidates by experimental route and direction-compatible validation logic. Therefore, the two rankings are not expected to be identical; for example, CD24 can rank highly by perturbation triage because it is upregulated, machine-learning-selected, cell-localized, and perturbation-priority supported, whereas COL14A1 ranks highly by integrated disease-association evidence.

Because the perturbation-triage score uses manually specified validation-planning weights, we added a weight-sensitivity analysis. Candidate ranks were recalculated with equal weights, leave-one-evidence-layer-out scores, and 1000 seeded +/-20% perturbations of the base weights. The five knockdown-screening candidates remained present across all sensitivity scenarios and each retained top-10 placement in at least 95% of scenarios, supporting the stability of the experimental class assignment while preserving the interpretation that the score is a prioritization aid rather than an efficacy estimate.

The final biological model organized the findings into five interpretable modules: macrophage/myeloid activation, stromal matrix remodeling, stromal antioxidant/metabolic loss, epithelial/ciliary remodeling, and an exact hsa-miR-375 target-axis hypothesis module. This organization separates biomarker classifier genes, disease-state markers, and putative perturbation candidates, which should not be treated as interchangeable categories.

Table 4 lists the top integrated candidate-prioritization results. Full scoring columns and additional candidates are provided in Additional file 7.

| gene_symbol | direction | main evidence layers | cellular context | follow-up class |
| --- | --- | --- | --- | --- |
| COL14A1 | upregulated | robust; ML; PPI; single-cell | fibroblast/myofibroblast | knockdown screening |
| GPX3 | downregulated | robust; ML; single-cell | fibroblast/myofibroblast | restoration or pathway marker |
| ASPN | upregulated | robust; ML; stromal module | pericyte/myofibroblast | knockdown screening |
| COL1A1 | upregulated | robust; PPI; single-cell | myofibroblast/stromal | context-dependent matrix/pathway screening |
| NECAB1 | downregulated | robust; ML | mesothelial/epithelial | pathway marker |
| SPP1 | upregulated | robust; PPI; single-cell; LR | macrophage/myeloid | context-dependent inflammatory/pathway screening |
| PTGFRN | upregulated | robust; ML; module | myofibroblast | knockdown screening |
| CD24 | upregulated | robust; ML; perturbation-priority proxy | AT2/transitional epithelial | knockdown screening |
| POSTN | upregulated | robust; PPI; LR; single-cell | myofibroblast/stromal | context-dependent matrix/pathway screening |
| CDH3 | upregulated | robust; ML | goblet/epithelial | knockdown screening |

The final biological model organized these findings into macrophage/myeloid activation, stromal matrix remodeling, stromal antioxidant/metabolic loss, epithelial/ciliary remodeling, and exact hsa-miR-375 target-axis hypothesis modules.

## Discussion
This multi-cohort transcriptomic study integrates mRNA, miRNA, disease-boundary testing, interaction-network, machine-learning, and single-cell evidence to triage IPF/fibrotic transcriptomic abnormalities for oligonucleotide-focused follow-up. Several points strengthen the resulting candidate set. First, mRNA and miRNA candidates were not selected from a single cohort alone; instead, discovery signals were required to show directionally consistent validation across independent datasets. Second, the machine-learning workflow used discovery-only feature prefiltering and nested resampling to reduce leakage, and the high-scoring output was explicitly bounded as a fibrotic ILD disease-state score rather than a diagnostic deployment model. Third, strict miRNA evidence gates and target-program stress tests prevented expansion beyond exact mature-miRNA hypotheses when broader target-set support was absent. Finally, donor-aware single-cell pseudobulk localized key candidates to biologically plausible lung compartments.

The intended novelty of this work is not de novo discovery of individual IPF genes. Several highly ranked candidates, including COL1A1, COL3A1, POSTN, SPP1, ASPN, and COL14A1, are already biologically plausible in fibrotic lung disease. The distinction from a biomarker-only analysis is the conversion of multi-cohort, direction-aware, cell-contextual evidence into a perturbation-triage framework. This framing separates knockdown-screening upregulated candidates, context-dependent matrix or inflammatory disease-state markers, restoration/pathway-marker candidates, exact miRNA-axis hypotheses, and the externally motivated TNIK bridge. The framework is designed to organize validation experiments, not to imply that any candidate is already therapeutically validated.

The supporting mechanistic extension adds a second layer of interpretation. Rather than treating the disease-state panel as a black-box biomarker list, marker-score and WGCNA-like module analysis showed that candidate genes concentrate in IPF-associated epithelial/ciliary and matrix-remodeling modules. Curated ligand-receptor scoring then connected these modules to altered intercellular signaling proxies, particularly MIF-CD74, COL14A1-ITGB1, CXCL12-CXCR4, POSTN-ITGB1, and SPP1-integrin/CD44 interactions. These results provide a coherent biological argument: the nominated candidates are embedded in disease-associated epithelial, endothelial, immune, and mesenchymal contexts. They remain hypothesis-generating and do not establish causal cell-cell signaling.

The robust mRNA set captured several expected IPF-associated programs. Upregulated matrix and fibroblast-related genes, including ASPN, COL14A1, COL1A1, COL3A1, POSTN, and THY1, were supported by PPI and single-cell evidence. The enrichment of cilium movement and axoneme assembly terms is also notable. Although several ontology labels reference sperm motility, the shared genes encode axonemal and microtubule-associated machinery relevant to motile cilia and epithelial biology. This supports a broader interpretation involving epithelial remodeling and mucociliary dysfunction rather than a reproductive process.

The miRNA analysis suggested several downregulated miRNAs with potential replacement or modulation relevance. hsa-miR-375, hsa-miR-30a, hsa-miR-30d, and hsa-miR-92a were directionally consistent across miRNA validation datasets, and miRTarBase integration connected these miRNAs to upregulated mRNA candidates. We explicitly separated exact mature-miRNA evidence from arm-agnostic evidence because these categories do not carry the same evidentiary weight. The hsa-miR-375 -> CLDN1/MNS1/RPGRIP1L axes were exact mature-miRNA matches and therefore suitable for main-text prioritization, whereas hsa-miR-30a and hsa-miR-92a axes remain exploratory until mature-arm expression and target repression are validated in IPF-relevant cells. The limited number of main-text miRNA axes therefore reflects strict evidence grading rather than an attempt to claim a broad regulatory miRNA network. These relationships should be considered prioritization hypotheses, because miRNA effects are context-dependent and can involve many targets.

From an oligonucleotide-focused validation perspective, the final panel contains two broad candidate classes. Upregulated mRNAs with stable validation and cell-type localization may be considered for knockdown screening with siRNA or antisense oligonucleotides. Downregulated miRNAs with inverse mRNA target relationships may motivate miRNA mimic or pathway-restoration experiments. The integrated priority table helps distinguish putative perturbation candidates such as COL14A1, ASPN, PTGFRN, CD24, and CDH3 from pathway-state markers such as GPX3, where restoration or pathway protection may be more biologically plausible than direct knockdown. However, delivery to fibrotic lung compartments, target-cell uptake, off-target effects, and disease-stage specificity remain major translational barriers.

The high external AUC values require cautious interpretation. Although validation cohorts were not used for feature selection or model tuning, public transcriptomic studies can still be influenced by platform, tissue-source, processing, clinical-annotation, and cohort-composition differences. Boundary stress tests supported interpretation as a fibrotic ILD disease-state score and a feature-prioritization layer, not as IPF-specific diagnostic deployment. Although leave-one-validation-cohort summaries, disease-boundary tests, and random robust-gene panel baselines were included, further prospective validation and independent clinical-cohort testing would be needed before clinical claims.

TNIK requires cautious interpretation. Public-data evidence placed TNIK in an IPF-associated coexpression module and showed differential expression in several validation cohorts, but it was not significant in the primary discovery cohort and was not selected into the final machine-learning panel. Therefore, a TNIK-targeting small nucleic acid can be connected to this manuscript only as an externally motivated translational bridge, preferably through Wnt/TNIK pathway readouts and co-culture rescue experiments. If the TNIK reagent is a knockdown oligonucleotide, the public IPF transcriptomic direction does not by itself support a simple disease-reversal claim; its value should be tested experimentally in the specific profibrotic cell state where TNIK activity is hypothesized to be pathogenic.

Experimental follow-up should first confirm candidate expression and perturbation efficiency in human lung fibroblasts, alveolar epithelial cells, and macrophage or monocyte-derived macrophage co-culture models. Perturbation experiments under TGF-beta1, epithelial-injury, or macrophage-conditioned contexts could then test candidate siRNA/ASO or miRNA mimic effects on matrix deposition, epithelial injury, inflammatory signaling, and Wnt/TNIK pathway readouts. These studies would be required to move from the present prioritization map to experimental candidate validation.

This study has limitations. It is retrospective and relies on public datasets generated across different platforms, tissues, protocols, and clinical annotation schemes. Although we applied strict sample matching and matrix QC, residual cohort effects cannot be eliminated fully. The machine-learning model was externally validated across public cohorts but not prospectively validated, and very high external AUC values should be interpreted with attention to possible cohort structure. Published model comparison was based on reported gene signatures rather than original trained model objects, because most prior studies did not provide deployable model files or coefficients. The miRNA layer used miRTarBase evidence grading and a targeted cross-check table. CLDN1 had additional TargetScan/literature support outside miRTarBase, but MNS1 and RPGRIP1L lacked independently verified prediction or CLIP support in the current version; therefore, all exact mature-miRNA axes remain prioritization hypotheses rather than confirmed IPF regulatory mechanisms. Finally, single-cell validation supports cellular localization and donor-aware pseudobulk sensitivity for selected candidates, but it does not prove causal function. Experimental perturbation in relevant human lung cell systems and in vivo models is required before therapeutic interpretation.

## Conclusions
The integrated analysis identified 280 robust mRNA candidates, 10 robust miRNA candidates, 22 opposite-direction miRNA-mRNA axes, and a 25-gene fibrotic disease-state panel with strong external separation but bounded diagnostic interpretation. Evidence grading retained three exact hsa-miR-375 target axes as miRNA-axis hypotheses and kept arm-agnostic axes as exploratory signals; target-program stress testing did not support expanding the miRNA interpretation beyond these strict axes. Perturbation triage highlighted validation classes that include knockdown-screening candidates such as COL14A1, ASPN, PTGFRN, CD24, and CDH3; context-dependent fibrotic disease-state markers such as COL1A1, SPP1, and POSTN; restoration/pathway-marker candidates such as GPX3 and NECAB1; exact hsa-miR-375 axis hypotheses; and an external TNIK bridge. The evidence converges on ciliary/epithelial remodeling, extracellular matrix activation, macrophage-associated SPP1 signaling, and stromal GPX3 loss as important IPF-associated patterns. These results provide a boundary-tested perturbation-triage map for downstream oligonucleotide-focused validation.

## Methods
### Data collection and annotation
Public GEO datasets relevant to human IPF mRNA, miRNA, and single-cell transcriptomics were downloaded and organized locally. Non-human data were excluded. For each bulk or miRNA dataset, sample annotation files were curated with consistent identifiers, disease labels, data type, dataset role, and inclusion status. The expression matrix column identifier was treated as the primary sample_id, while GEO accession identifiers were retained separately when available.

### Expression extraction and three-layer QC
Expression matrices were extracted from GEO series matrix files or supplementary count files [1,2]. Each dataset underwent three independent QC checks: annotation completeness, sample cross-match between annotation and expression columns, and expression-matrix integrity including numeric coercion, missingness, duplicate feature checks, and scale assessment. Only datasets passing all three checks were used in downstream bulk or miRNA analyses.

Matrix type was adjudicated from the GEO source, supplementary filename, expression range, zero fraction, and value distribution. GSE150910 and GSE92592 were treated as supplementary gene-count matrices and analyzed with edgeR-limma voom. GSE110147, GSE21394, and GSE32537 were treated as processed log-scale or normalized intensity matrices and analyzed with limma. GSE27430, GSE32538, and GSE53845 contained centered or transformed expression values, including negative values; these datasets were analyzed with limma, and cross-cohort validation emphasized FDR and direction consistency because absolute logFC values in centered matrices are not directly interchangeable with raw count-derived log fold changes. Dataset-level matrix-processing decisions are provided in Additional file 1.

### Differential expression
Differential expression was performed in R using limma for processed log-scale or normalized matrices and edgeR with limma-voom for supplementary count matrices [14-16]. The primary contrast was IPF versus control. Discovery significance used FDR < 0.05 and absolute log fold change >= 1. Probe IDs, transcript IDs, and miRNA probe IDs were harmonized to gene symbols or miRNA names before cross-cohort comparison. Mature-arm identifiers were retained where available; miRNA-family or arm-unspecified records were flagged as arm-agnostic.

### Robust candidate selection
GSE32537 served as the mRNA discovery cohort and GSE32538 as the miRNA discovery cohort. mRNA validation used GSE110147, GSE150910, GSE53845, and GSE92592. miRNA validation used GSE21394 and GSE27430. Strict robust mRNAs required same-direction validation in at least two datasets. For log-scale or count-derived datasets, validation support used FDR < 0.05 and absolute log fold change >= 1; centered/transformed datasets were treated as direction-supporting validation layers with cautious fold-change interpretation because their absolute logFC values are not directly interchangeable with raw count-derived log fold changes. Strict robust miRNAs required at least one same-direction validation dataset meeting the miRNA validation rule.

### miRNA-mRNA target-axis integration
Robust miRNAs and mRNAs were integrated with the human miRTarBase 2025 v10 interaction table. Candidate axes were retained when the miRNA and mRNA discovery log fold changes were in opposite directions and the target gene was present in the robust mRNA set. Exact mature-miRNA matches and arm-agnostic matches were annotated separately.

### Evidence grading and integrated candidate prioritization
miRNA-mRNA axes were graded by match specificity and experimental-support category. Exact mature-miRNA matches were eligible for main-text prioritization, while arm-agnostic matches were retained as exploratory hypotheses requiring mature-arm-specific validation. Robust mRNA candidates were then scored using a transparent additive framework incorporating robust-expression support, machine-learning panel evidence, PPI hub status, enrichment membership, miRNA-axis evidence, single-cell disease-control signal, and oligonucleotide strategy compatibility. The score was used to rank follow-up candidates, not to estimate causal effect size.

### Perturbation-triage validation-planning score
The perturbation-triage validation-planning score was built from the integrated priority table, miRNA-axis evidence, single-cell localization summaries, and coexpression-neighborhood perturbation-priority proxy. It was defined as an additive validation-planning score rather than a therapeutic-readiness score. The components were disease direction (1.0 for upregulated candidates, 0.35 for downregulated candidates), cross-cohort support (same-direction significant validation count divided by four and capped at 1.0), single-cell localization (number of supporting single-cell datasets divided by two and capped at 1.0), positive perturbation-priority support normalized to 0-1, machine-learning selection frequency, and miRNA-axis support (exact mature-miRNA evidence plus 0.25-weighted arm-agnostic evidence, capped at 1.0). Directionality and cross-cohort support were assigned the highest weights because they determine whether an oligonucleotide perturbation route is conceptually compatible and reproducible across independent cohorts; cell localization and perturbation-priority evidence were treated as intermediate validation layers, while machine-learning and miRNA-axis evidence were used as supportive layers. The final score used weights of 1.30 for disease direction, 1.20 for cross-cohort support, 1.00 for cell localization, 1.00 for perturbation-priority support, 0.90 for machine-learning evidence, and 0.80 for miRNA-axis evidence, with a 0.15 context penalty for broad matrix, inflammatory, or downregulated markers where simple knockdown interpretation is less direct. In formula form: score = 1.30 x direction + 1.20 x cross-cohort + 1.00 x cell localization + 1.00 x perturbation-priority + 0.90 x machine-learning + 0.80 x miRNA-axis - 0.15 x context penalty. The 0.35 downregulated-candidate direction value and 0.15 context penalty were chosen to encode validation-priority logic rather than learned from outcome data.

Candidate classes were assigned using the same evidence layers and biological caveats. Upregulated candidates with robust validation and a plausible knockdown-screening route were assigned to a knockdown-screening class. Broad matrix or macrophage-associated genes with important disease biology but less straightforward perturbation interpretation were assigned to a context-dependent class. Downregulated genes were assigned to restoration or pathway-marker classes. Exact hsa-miR-375 target genes were assigned to the miRNA-axis hypothesis class, and TNIK was assigned to an external-bridge class rather than scored as a primary discovery candidate. The context penalty was binary and non-accumulating; a candidate received one 0.15 penalty if it fell into any broad matrix, inflammatory, or downregulated marker category with less direct simple-knockdown interpretation. The index was intended to organize experimental follow-up, not to infer therapeutic readiness or causal efficacy.

Perturbation-triage score sensitivity was evaluated using three perturbation designs: equal-weight scoring, leave-one-evidence-layer-out scoring, and 1000 seeded random perturbations in which each base weight was multiplied by a factor sampled uniformly from 0.8 to 1.2. Candidate ranks and top-rank frequencies were summarized across all scenarios. This analysis tested whether the proposed validation classes were stable to modest changes in the manually specified score weights.

### Marker-score modules, curated ligand-receptor scoring, and perturbation-priority proxy
Bulk immune, stromal, epithelial/ciliary, endothelial, and Wnt/TNIK marker scores were computed as within-dataset mean z scores of curated marker genes. In the GSE32537 discovery cohort, highly variable genes plus forced inclusion of robust candidates, machine-learning panel genes, ligand-receptor genes, and TNIK were clustered by correlation distance to generate WGCNA-like coexpression modules. Module eigengenes were calculated by singular-value decomposition and correlated with disease status and marker scores. Single-cell communication potential was approximated by curated ligand-receptor scoring: mean normalized ligand expression in sender cell types was multiplied by mean normalized receptor expression in receiver cell types, and IPF-control score differences were calculated after excluding multiplet, outlier, cell-cycle, and mitochondrial tRNA labels. A coexpression-neighborhood perturbation-priority proxy ranked candidate perturbations by discovery-cohort coexpression-neighborhood structure and target directionality. These analyses are mechanism-generating extensions and should not be interpreted as outputs from full immune-deconvolution, formal cell-communication, or network-knockout packages.

### Enrichment and PPI analysis
GO, KEGG, and Reactome enrichment analyses were performed with clusterProfiler and ReactomePA using the gene-level background tested in the GSE32537 discovery analysis. STRING protein interaction networks were queried for Homo sapiens at medium and high confidence thresholds. Hubs were ranked by combined degree, weighted degree, betweenness, and closeness metrics.

### Machine learning and published-signature comparison
Gene-level matrices were generated by collapsing multiple probes or transcripts to gene symbols using within-dataset variance. The main machine-learning analysis used only discovery-significant mRNA features from GSE32537 that were common across external validation datasets. Seven model families were evaluated: lasso logistic regression, Elastic Net, linear support vector machine, radial-basis support vector machine, random forest, gradient boosting, and a small multilayer perceptron. Imputation used median replacement, scaling used StandardScaler, and feature selection used SelectKBest with the f_classif score. Candidate feature counts were 15, 30, 50, and 100 where permitted by the available feature universe. Hyperparameter tuning was performed inside GridSearchCV with ROC AUC scoring. Logistic and linear-SVM C values were 0.03, 0.1, 0.3, and 1.0; Elastic Net additionally tested l1_ratio values of 0.2, 0.5, and 0.8; RBF SVM tested C values of 0.1, 1.0, and 3.0 with gamma set to scale; random forest tested max_depth values of 3, 5, and unrestricted with min_samples_leaf values of 1 and 3; gradient boosting tested n_estimators of 80 and 150, learning rates of 0.03 and 0.08, and max_depth values of 1 and 2; the neural-network classifier tested hidden-layer sizes of 16, 16/8, and 32/16 with alpha values of 0.001, 0.01, and 0.1. The internal design used repeated stratified five-fold outer cross-validation repeated 10 times, with three-fold inner cross-validation for tuning. Final models were refit in GSE32537 using five-fold internal tuning, and external validation was performed only after model selection. Label permutation was used as a negative control.

Published IPF signatures were evaluated as fixed gene-set comparators because deployable trained model objects or full coefficients were generally unavailable [29-34]. Each comparator was trained in GSE32537 and externally validated on the same four mRNA validation cohorts.

Model interpretability and utility checks were performed on pooled external validation samples using permutation importance, calibration curves, and decision-curve analysis. These outputs were used to assess whether model predictions were driven by biologically plausible features and whether predicted probabilities had reviewable calibration and threshold-utility behavior.

Sensitivity analyses used the locked Elastic Net external predictions for leave-one-validation-cohort summaries. For the random robust-gene baseline, 500 random panels with the same size as the observed 25-gene panel were sampled without replacement from the common robust mRNA feature universe. Each random panel was trained in GSE32537 using the final Elastic Net hyperparameters and externally evaluated in the four validation cohorts. The observed 25-gene panel was refit with the same fixed Elastic Net specification to provide an algorithm-matched comparison against the random panels.

External disease-state boundary tests were added to define the interpretation of high external AUC values. First, excluded NSIP and mixed IPF-NSIP samples from GSE110147 were scored using a final-panel regularized refit without including these samples in training or model selection. Second, 500 matched random discovery-feature panels were generated by sampling 25 features from the common discovery-feature universe while matching discovery direction, validation-support count, discovery absolute logFC stratum, and discovery FDR stratum as closely as possible to the observed panel. Third, external metrics were summarized after excluding validation cohorts with perfect ROC AUC. Fourth, a pooled external logistic regression tested disease status as a function of the Elastic Net logit score with validation cohort included as a fixed effect. These analyses were intended to bound classifier interpretation rather than to support clinical deployment.

### miRNA target-program stress tests
Because exact mature-miRNA axis criteria were intentionally conservative, robust downregulated miRNAs were also evaluated at the target-set level. Human miRTarBase functional interactions were harmonized to exact or arm-recoverable mature-miRNA names and compared against IPF-upregulated robust mRNAs using the discovery mRNA background. For each robust downregulated miRNA, Fisher exact tests evaluated whether its miRTarBase target set was enriched among IPF-upregulated robust mRNAs, followed by Benjamini-Hochberg correction across tested miRNAs. For hsa-miR-375, a target release-like score compared discovery mRNA logFC values for hsa-miR-375 targets versus non-target genes and used 5000 seeded label permutations to estimate a one-sided enrichment P value. A relaxed exploratory axis set was also generated for supplementary audit by allowing exact or arm-recoverable miRNA matches with opposite-direction robust mRNA candidates. Paired miRNA-mRNA inverse correlation was not attempted because GSE32537 and GSE32538 did not provide explicit one-to-one shared donor identifiers in the available metadata.

### Single-cell validation
Single-cell sparse expression matrices were processed by streaming matrix entries to avoid excessive memory use [25,26]. Candidate genes were evaluated by cell type and disease group. Cell labels corresponding to Multiplet, Outlier, MT-tRNAs, or CellCycle categories were excluded from manuscript-level summaries, and cell-type comparisons required at least 20 IPF and 20 control cells. These single-cell analyses were used as localization and prioritization evidence. They were not treated as causal differential-expression tests.

For donor-aware sensitivity analysis, core candidates were aggregated at the donor/sample x broad-celltype level in GSE135893 and GSE136831 after the same non-informative label filtering used for manuscript-level single-cell summaries. GSE135893 donor IDs were taken from Sample_Name and GSE136831 donor IDs from Subject_Identity. For each donor-celltype-gene combination, raw target-gene counts were summed and normalized by the summed cell-level library size, then log1p transformed. IPF-control differences were estimated across donors using two-sided Welch tests within each dataset, broad cell type, and gene, with Benjamini-Hochberg correction across all tested core candidate x dataset x broad-celltype pseudobulk comparisons. This pseudobulk layer was restricted to core candidates and used as validation support rather than a full single-cell differential-expression discovery analysis. Figure 5C displayed direction-consistent comparisons with at least 10 IPF and 10 control donors and BH-adjusted FDR < 0.01. No additional per-donor minimum cell-count threshold was imposed after broad-celltype aggregation; donor counts and per-donor cell-count summaries are provided for audit. The n values shown in Figure 5C denote the number of IPF/control donors represented in the indicated broad-celltype category after filtering. The original-label to broad-celltype mapping and broad-celltype cell-count summaries are provided in Additional file 9.

### Software and reproducibility
Analyses were implemented with R, Python, Bioconductor packages including limma, edgeR, clusterProfiler, ReactomePA, org.Hs.eg.db, and Python packages including pandas, numpy, scikit-learn, scipy, matplotlib, and openpyxl. Software versions, random seeds, script order, and input/output mapping are provided in Additional file 11. All generated outputs are stored in the project results directory and are traceable through the numbered scripts in the scripts directory.

## Declarations
### Ethics approval and consent to participate
Not applicable. This study analyzed publicly available de-identified datasets.

### Consent for publication
Not applicable.

### Availability of data and materials
All raw datasets analyzed in this study are publicly available from the Gene Expression Omnibus under the accession numbers reported in Table 1 and the single-cell validation section. Processed non-sensitive outputs supporting the conclusions are included as Additional files 1-10 and Additional files 12-14. The numbered analysis scripts, README file, software environment notes, and reproducibility instructions are provided as Additional file 11 and are publicly available at https://github.com/osbornzhou/IPF_oligo.

### Competing interests
The authors declare that they have no competing interests.

### Funding
No specific funding was received for this study.

### Authors' contributions
Yunyi Zhou curated the datasets, performed the computational analyses, generated the figures and tables, and drafted the manuscript. Yanli Zhang supervised the study design, interpretation, and manuscript revision. Both authors read and approved the final manuscript.

### Acknowledgements
The authors thank the investigators who generated and shared the public GEO datasets used in this study.

### Use of AI-assisted tools
AI-assisted drafting and coding support were used to organize analysis scripts and prepare an initial manuscript draft. The authors are responsible for verifying all analyses, interpretations, citations, and the final submitted text.

### Additional files
Additional file 1. XLSX. Dataset metadata and QC. Bulk, miRNA, and single-cell annotation and expression-matrix QC tables.

Additional file 2. XLSX. Differential expression results. Annotated differential-expression summaries and QC outputs.

Additional file 3. XLSX. Robust mRNA and miRNA candidates. Discovery-validation robust candidate tables, direction matrices, and dataset-specific validation-support rules.

Additional file 4. XLSX. miRNA-mRNA candidate axes. miRTarBase-derived opposite-direction candidate axes, evidence grading, exact-axis support table, target-set stress tests, hsa-miR-375 target release-like score, relaxed exploratory axes, and paired-correlation audit.

Additional file 5. XLSX. GO/KEGG/Reactome enrichment results. Enrichment results, gene-set inputs, and enrichment QC tables.

Additional file 6. XLSX. STRING PPI and hub results. STRING mapping, network nodes, edges, hub metrics, and QC tables.

Additional file 7. XLSX. Machine-learning performance, feature stability, and perturbation triage. Model performance, feature-selection stability, predictions, hyperparameter grids, final Elastic Net specification, candidate-prioritization outputs, perturbation-triage validation-planning score, and score weight-sensitivity analysis.

Additional file 8. XLSX. Machine-learning sensitivity analyses. Leave-one-validation-cohort sensitivity, random robust-gene panel baseline outputs, disease-control score stress test, matched random discovery-feature panel baseline, non-perfect-cohort summary, and cohort-adjusted score test.

Additional file 9. XLSX. Single-cell localization summaries. Single-cell target-expression summaries, donor-aware pseudobulk validation, broad-celltype mapping, localization tables, and validation QC.

Additional file 10. XLSX. Mechanistic extension outputs. Marker scores, coexpression modules, curated ligand-receptor scoring, coexpression-neighborhood perturbation-priority proxy, and TNIK bridge tables.

Additional file 11. ZIP. Analysis code archive. Numbered analysis scripts, README file, software environment notes, and reproducibility instructions.

Additional file 12. PDF. Donor-level pseudobulk plots. Supplementary donor-level pseudobulk dot/median plots for the core single-cell validation comparisons shown in Figure 5C.

Additional file 13. PDF. AUC disease-state stress-test plots. Supplementary disease-control score distribution and matched random discovery-feature panel baseline for external ML boundary assessment.

Additional file 14. PDF. miRNA target-program stress-test plots. Supplementary miRTarBase target-set stress test and hsa-miR-375 target release-like score used to bound miRNA interpretation.

## Figure and table plan
Figure 1. Perturbation-triage study design from GEO data collection, sample-label audit, matrix QC, robust transcriptomic screening, disease-boundary testing, miRNA evidence gating, donor-aware single-cell localization, and perturbation-triage classification.

Figure 2. Cross-cohort transcriptomic evidence base. Suggested panels: discovery volcano plot; cross-cohort logFC/effect estimates and direction matrix for top robust mRNAs and robust miRNAs.

Figure 3. Disease-boundary stress tests. Excluded GSE110147 NSIP and mixed IPF-NSIP samples are scored without training use, and matched random discovery-feature panels contextualize the final-panel refit.

Figure 4. miRNA evidence-gate stress tests. Exact hsa-miR-375 axes are retained under strict criteria, while target-set enrichment and hsa-miR-375 release-like score define why the miRNA interpretation is not expanded to a broad regulatory program.

Figure 5. Donor-aware single-cell localization of candidate signals. Suggested panels: cell-level disease-control delta heatmaps, largest candidate shifts, and donor-aware pseudobulk validation with IPF/control donor counts and FDR labels.

Figure 6. Perturbation-triage map. Candidate genes and miRNA-axis targets are grouped by validation route rather than therapeutic readiness.

Additional Figure S1. Integrated candidate-prioritization score for the top robust candidates.

Additional Figure S2. Machine-learning permutation importance, calibration, and decision-curve checks in pooled external validation samples.

Additional Figure S3. Biological interpretation model linking robust transcriptomic signals, cell context, and oligonucleotide development hypotheses.

Additional Figure S4. WGCNA-like module-trait heatmap showing IPF-associated epithelial/ciliary and matrix-remodeling modules.

Additional Figure S5. Curated ligand-receptor score changes in IPF single-cell datasets.

Additional Figure S6. Coexpression-neighborhood perturbation-priority ranking for oligonucleotide follow-up candidates including TNIK.

Additional Figure S7. TNIK evidence summary across bulk cohorts.

Additional Figure S8. Machine-learning sensitivity checks showing random robust 25-gene panel baseline performance.

Additional Figure S9. miRNA-mRNA evidence grading workflow showing strict mature-miRNA evidence prioritization and exploratory arm-agnostic axes.

Additional Figure S10. Donor-aware pseudobulk validation of core single-cell candidates. Red bars indicate IPF-increased pseudobulk expression and blue bars indicate IPF-decreased pseudobulk expression. Labels indicate IPF/control donor counts represented in each broad-celltype category and Benjamini-Hochberg adjusted FDR values.

Additional Figure S11. Perturbation-triage score weight-sensitivity analysis. Lines indicate the rank range across equal-weight, leave-one-layer-out, and +/-20% weight-perturbation scenarios; lower rank indicates higher prioritization.

Additional Figure S12. Donor-level pseudobulk dot/median plots for the core comparisons displayed in Figure 5C. Each dot represents one donor-level pseudobulk value; vertical bars mark group medians.

Additional file 13 contains the supporting PDF version and machine-readable source outputs for the disease-boundary stress tests shown in Figure 3.

Additional file 14 contains the supporting PDF version and machine-readable source outputs for the miRNA evidence-gate stress tests shown in Figure 4.

Table 1. Quality-controlled datasets.
Table 2. Exact mature-miRNA candidate axes.
Table 3. Machine-learning external validation and sensitivity summary.
Table 4. Integrated candidate-prioritization shortlist.
Supplementary Table 1. Full merged annotation and matrix QC.
Supplementary Table 2. Full differential expression results.
Supplementary Table 3. Robust mRNA and miRNA candidates.
Supplementary Table 4. miRNA-mRNA axes.
Supplementary Table 5. Enrichment and PPI outputs.
Supplementary Table 6. Machine-learning predictions, performance, and feature stability.
Supplementary Table 7. Single-cell target localization.
Supplementary Table 8. Evidence-graded miRNA-mRNA axes and integrated candidate-prioritization outputs.
Supplementary Table 9. Mechanistic extension, curated ligand-receptor scoring, perturbation-priority proxy, and TNIK bridge outputs.
Supplementary Table 10. Machine-learning sensitivity analyses including leave-one-validation-cohort summaries and random robust-gene panel baselines.

## References
1. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res. 2002;30:207-210. doi:10.1093/nar/30.1.207.
2. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets-update. Nucleic Acids Res. 2013;41:D991-D995. doi:10.1093/nar/gks1193.
3. GEO Series GSE32537. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32537. Accessed 25 May 2026.
4. GEO Series GSE32538. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32538. Accessed 25 May 2026.
5. GEO Series GSE110147. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110147. Accessed 25 May 2026.
6. Cecchini MJ, Hosein K, Howlett CJ, Joseph M, Mura M. Comprehensive gene expression profiling identifies distinct and overlapping transcriptional profiles in non-specific interstitial pneumonia and idiopathic pulmonary fibrosis. Respir Res. 2018;19:153. doi:10.1186/s12931-018-0857-1.
7. GEO Series GSE150910. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150910. Accessed 25 May 2026.
8. GEO Series GSE53845. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53845. Accessed 25 May 2026.
9. GEO Series GSE92592. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92592. Accessed 25 May 2026.
10. GEO Series GSE21394. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21394. Accessed 25 May 2026.
11. GEO Series GSE27430. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27430. Accessed 25 May 2026.
12. GEO Series GSE135893. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893. Accessed 25 May 2026.
13. GEO Series GSE136831. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831. Accessed 25 May 2026.
14. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47. doi:10.1093/nar/gkv007.
15. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26:139-140. doi:10.1093/bioinformatics/btp616.
16. Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014;15:R29. doi:10.1186/gb-2014-15-2-r29.
17. Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. Innovation (Camb). 2021;2:100141. doi:10.1016/j.xinn.2021.100141.
18. Yu G, He QY. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Mol Biosyst. 2016;12:477-479. doi:10.1039/C5MB00663E.
19. Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Res. 2023;51:D638-D646. doi:10.1093/nar/gkac1000.
20. Huang HY, Lin YC, Li J, et al. miRTarBase 2020: updates to the experimentally validated microRNA-target interaction database. Nucleic Acids Res. 2020;48:D148-D154. doi:10.1093/nar/gkz896.
21. miRTarBase. https://mirtarbase.cuhk.edu.cn/. Accessed 25 May 2026.
22. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. J Mach Learn Res. 2011;12:2825-2830.
23. Raghu G, Remy-Jardin M, Richeldi L, et al. Idiopathic pulmonary fibrosis: an update of the ATS/ERS/JRS/ALAT clinical practice guideline. Am J Respir Crit Care Med. 2022;205:e18-e47. doi:10.1164/rccm.202202-0399ST.
24. Lederer DJ, Martinez FJ. Idiopathic pulmonary fibrosis. N Engl J Med. 2018;378:1811-1823. doi:10.1056/NEJMra1705751.
25. Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis. Sci Adv. 2020;6:eaba1972. doi:10.1126/sciadv.aba1972.
26. Adams TS, Schupp JC, Poli S, et al. Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in idiopathic pulmonary fibrosis. Sci Adv. 2020;6:eaba1983. doi:10.1126/sciadv.aba1983.
27. Roberts TC, Langer R, Wood MJA. Advances in oligonucleotide drug delivery. Nat Rev Drug Discov. 2020;19:673-694. doi:10.1038/s41573-020-0075-7.
28. Crooke ST, Witztum JL, Bennett CF, Baker BF. RNA-targeted therapeutics. Cell Metab. 2018;27:714-739. doi:10.1016/j.cmet.2018.03.004.
29. Zeng Y, Huang J, Guo R, Cao S, Yang H, Ouyang W. Identification and validation of metabolism-related hub genes in idiopathic pulmonary fibrosis. Front Genet. 2023;14:1058582. doi:10.3389/fgene.2023.1058582.
30. Li Z, Wang S, Zhao H, et al. Artificial neural network identified the significant genes to distinguish idiopathic pulmonary fibrosis. Sci Rep. 2023;13:1225. doi:10.1038/s41598-023-28536-w.
31. Wan H, Huang X, Cong P, et al. Identification of hub genes and pathways associated with idiopathic pulmonary fibrosis via bioinformatics analysis. Front Mol Biosci. 2021;8:711239. doi:10.3389/fmolb.2021.711239.
32. Fanidis D, Pezoulas VC, Fotiadis DI, Aidinis V. An explainable machine learning-driven proposal of pulmonary fibrosis biomarkers. Comput Struct Biotechnol J. 2023;21:2305-2315. doi:10.1016/j.csbj.2023.03.043.
33. Rosas IO, Richards TJ, Konishi K, et al. MMP1 and MMP7 as potential peripheral blood biomarkers in idiopathic pulmonary fibrosis. PLoS Med. 2008;5:e93. doi:10.1371/journal.pmed.0050093.
34. Chen B, Lu H, Lu J, Yuan J. Shapley additive explanations based feature selection reveals CXCL14 as a key immune-related gene in predicting idiopathic pulmonary fibrosis. Front Med. 2025;12:1608078. doi:10.3389/fmed.2025.1608078.
35. Yoda S, Soejima K, Hamamoto J, et al. Claudin-1 is a novel target of miR-375 in non-small-cell lung cancer. Lung Cancer. 2014;85:366-372. doi:10.1016/j.lungcan.2014.06.009.
36. Yang IV, Coldren CD, Leach SM, et al. Expression of cilium-associated genes defines novel molecular subtypes of idiopathic pulmonary fibrosis. Thorax. 2013;68:1114-1121. doi:10.1136/thoraxjnl-2012-202943.
37. Furusawa H, Cardwell JH, Okamoto T, et al. Chronic hypersensitivity pneumonitis, an interstitial lung disease with distinct molecular signatures. Am J Respir Crit Care Med. 2020;202:1430-1444. doi:10.1164/rccm.202001-0134OC.
38. DePianto DJ, Chandriani S, Abbas AR, et al. Heterogeneous gene expression signatures correspond to distinct lung pathologies and biomarkers of disease severity in idiopathic pulmonary fibrosis. Thorax. 2015;70:48-56. doi:10.1136/thoraxjnl-2013-204596.
39. Schafer MJ, White TA, Iijima K, et al. Cellular senescence mediates fibrotic pulmonary disease. Nat Commun. 2017;8:14532. doi:10.1038/ncomms14532.
40. Cho JH, Gelinas R, Wang K, et al. Systems biology of interstitial lung diseases: integration of mRNA and microRNA expression changes. BMC Med Genomics. 2011;4:8. doi:10.1186/1755-8794-4-8.
41. Milosevic J, Pandit K, Magister M, et al. Profibrotic role of miR-154 in pulmonary fibrosis. Am J Respir Cell Mol Biol. 2012;47:879-887. doi:10.1165/rcmb.2011-0377OC.
