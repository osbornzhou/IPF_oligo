#!/usr/bin/env Rscript

user_library <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", paste0(R.version$major, ".", R.version$minor))
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
  library(ggplot2)
})

project_dir <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
annotation_path <- file.path(project_dir, "metadata", "all_bulk_mirna_annotation.csv")
expression_qc_path <- file.path(project_dir, "metadata", "expression_matrix_qc.csv")
expression_dir <- file.path(project_dir, "data_processed", "expression")
output_dir <- file.path(project_dir, "results", "differential_expression")
plot_dir <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
}

stop_if_missing(annotation_path)
stop_if_missing(expression_qc_path)

annotation <- fread(annotation_path, na.strings = c("", "NA", "N/A", "na"))
expression_qc <- fread(expression_qc_path, na.strings = c("", "NA", "N/A", "na"))

expression_qc <- expression_qc[triple_qc_pass %in% c(TRUE, "True", "TRUE", "true", 1)]
if (nrow(expression_qc) == 0) {
  stop("No expression matrices passed triple QC.")
}

clean_feature_id <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- paste0("feature_", seq_len(sum(is.na(x) | x == "")))
  make.unique(x, sep = "_dup")
}

prepare_annotation <- function(series_id) {
  sid <- series_id
  pheno <- annotation[
    series_id == sid &
      include %in% c("Yes", "YES", "yes", TRUE, "TRUE", "True") &
      group %in% c("IPF", "Control")
  ]
  pheno[, group := factor(group, levels = c("Control", "IPF"))]
  pheno <- pheno[!is.na(sample_id) & !duplicated(sample_id)]
  pheno
}

read_expression <- function(path) {
  expr_dt <- fread(path)
  if (ncol(expr_dt) < 3) {
    stop("Expression matrix has too few columns: ", path)
  }
  feature_col <- names(expr_dt)[1]
  feature_id <- clean_feature_id(expr_dt[[feature_col]])
  expr_dt[, (feature_col) := NULL]
  mat <- as.matrix(expr_dt)
  mode(mat) <- "numeric"
  rownames(mat) <- feature_id
  mat
}

select_batch <- function(pheno) {
  if (!"batch" %in% names(pheno)) {
    return(NULL)
  }
  batch <- as.character(pheno$batch)
  batch[is.na(batch) | batch == ""] <- NA_character_
  if (sum(!is.na(batch)) != length(batch)) {
    return(NULL)
  }
  if (length(unique(batch)) < 2) {
    return(NULL)
  }
  batch_tab <- table(batch)
  if (any(batch_tab < 2)) {
    return(NULL)
  }
  cross_tab <- table(batch, pheno$group)
  if (any(rowSums(cross_tab > 0) < 2)) {
    return(NULL)
  }
  factor(batch)
}

build_design <- function(pheno) {
  batch <- select_batch(pheno)
  if (is.null(batch)) {
    design <- model.matrix(~ group, data = pheno)
    contrast <- "groupIPF"
    formula_label <- "~ group"
    batch_adjusted <- FALSE
  } else {
    design <- model.matrix(~ group + batch, data = cbind(pheno, batch = batch))
    contrast <- "groupIPF"
    formula_label <- "~ group + batch"
    batch_adjusted <- TRUE
  }
  if (qr(design)$rank < ncol(design)) {
    design <- model.matrix(~ group, data = pheno)
    formula_label <- "~ group"
    batch_adjusted <- FALSE
  }
  list(design = design, coef = contrast, formula = formula_label, batch_adjusted = batch_adjusted)
}

filter_normalized_matrix <- function(mat) {
  finite_rows <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[finite_rows, , drop = FALSE]
  variances <- matrixStats::rowVars(mat)
  mat[is.finite(variances) & variances > 0, , drop = FALSE]
}

run_limma_normalized <- function(mat, pheno, design_info) {
  mat <- filter_normalized_matrix(mat)
  fit <- lmFit(mat, design_info$design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = design_info$coef, number = Inf, sort.by = "P")
  res$feature_id <- rownames(res)
  setDT(res)
  setcolorder(res, c("feature_id", setdiff(names(res), "feature_id")))
  list(result = res, tested_features = nrow(mat), method = "limma")
}

run_limma_voom <- function(mat, pheno, design_info) {
  mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
  mat[mat < 0] <- 0
  mat <- round(mat)
  dge <- DGEList(counts = mat, group = pheno$group)
  keep <- filterByExpr(dge, design = design_info$design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design_info$design, plot = FALSE)
  fit <- lmFit(v, design_info$design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = design_info$coef, number = Inf, sort.by = "P")
  res$feature_id <- rownames(res)
  setDT(res)
  setcolorder(res, c("feature_id", setdiff(names(res), "feature_id")))
  list(result = res, tested_features = nrow(v$E), method = "edgeR_limma_voom")
}

plot_volcano <- function(res, series_id, method, path) {
  plot_dt <- copy(res)
  plot_dt[, neg_log10_fdr := -log10(pmax(adj.P.Val, .Machine$double.xmin))]
  plot_dt[, significant := adj.P.Val < 0.05 & abs(logFC) >= 1]
  p <- ggplot(plot_dt, aes(x = logFC, y = neg_log10_fdr, color = significant)) +
    geom_point(alpha = 0.7, size = 0.9) +
    scale_color_manual(values = c("FALSE" = "#6f7785", "TRUE" = "#b42318")) +
    labs(
      title = paste0(series_id, " IPF vs Control"),
      subtitle = method,
      x = "log2 fold change",
      y = "-log10 adjusted P value",
      color = "FDR < 0.05\n|logFC| >= 1"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  ggsave(path, p, width = 6.4, height = 4.8, dpi = 160)
}

qc_rows <- list()

for (i in seq_len(nrow(expression_qc))) {
  q <- expression_qc[i]
  series_id <- q$series_id
  message("Running differential expression for ", series_id)

  expr_path <- file.path(expression_dir, paste0(series_id, "_expression_in_annotation_order.csv.gz"))
  if (!file.exists(expr_path)) {
    expr_path <- q$annotation_order_output_path
  }
  stop_if_missing(expr_path)

  pheno <- prepare_annotation(series_id)
  mat <- read_expression(expr_path)

  matched_samples <- intersect(colnames(mat), pheno$sample_id)
  pheno <- pheno[match(matched_samples, sample_id)]
  mat <- mat[, matched_samples, drop = FALSE]

  if (nrow(pheno) < 4 || length(unique(pheno$group)) < 2) {
    warning("Skipping ", series_id, ": insufficient IPF/Control samples.")
    next
  }

  design_info <- build_design(pheno)
  scale_guess <- as.character(q$value_scale_guess)
  is_count_like <- scale_guess %in% c("count_or_abundance")

  result <- if (is_count_like) {
    run_limma_voom(mat, pheno, design_info)
  } else {
    run_limma_normalized(mat, pheno, design_info)
  }

  de <- result$result
  de[, series_id := series_id]
  de[, data_type := as.character(q$data_type)]
  de[, dataset_role := as.character(q$dataset_role)]
  de[, method := result$method]
  de[, comparison := "IPF_vs_Control"]
  setcolorder(
    de,
    c("series_id", "data_type", "dataset_role", "method", "comparison", "feature_id",
      setdiff(names(de), c("series_id", "data_type", "dataset_role", "method", "comparison", "feature_id")))
  )

  result_path <- file.path(output_dir, paste0(series_id, "_IPF_vs_Control_de_results.csv"))
  significant_path <- file.path(output_dir, paste0(series_id, "_IPF_vs_Control_de_significant_fdr0.05_logfc1.csv"))
  plot_path <- file.path(plot_dir, paste0(series_id, "_volcano.png"))

  fwrite(de, result_path)
  fwrite(de[adj.P.Val < 0.05 & abs(logFC) >= 1], significant_path)
  plot_volcano(de, series_id, result$method, plot_path)

  qc_rows[[series_id]] <- data.table(
    series_id = series_id,
    data_type = as.character(q$data_type),
    dataset_role = as.character(q$dataset_role),
    method = result$method,
    value_scale_guess = scale_guess,
    design_formula = design_info$formula,
    batch_adjusted = design_info$batch_adjusted,
    input_features = nrow(mat),
    tested_features = result$tested_features,
    samples_total = nrow(pheno),
    ipf_samples = sum(pheno$group == "IPF"),
    control_samples = sum(pheno$group == "Control"),
    significant_fdr_0_05 = sum(de$adj.P.Val < 0.05),
    significant_fdr_0_05_logfc_1 = sum(de$adj.P.Val < 0.05 & abs(de$logFC) >= 1),
    top_feature = de$feature_id[1],
    top_logFC = de$logFC[1],
    top_adj_p = de$adj.P.Val[1],
    result_path = normalizePath(result_path, winslash = "/", mustWork = FALSE),
    significant_path = normalizePath(significant_path, winslash = "/", mustWork = FALSE),
    volcano_path = normalizePath(plot_path, winslash = "/", mustWork = FALSE)
  )
}

de_qc <- rbindlist(qc_rows, fill = TRUE)
fwrite(de_qc, file.path(output_dir, "differential_expression_qc.csv"))
fwrite(
  rbindlist(lapply(
    list.files(output_dir, pattern = "_de_significant_fdr0.05_logfc1\\.csv$", full.names = TRUE),
    fread
  ), fill = TRUE),
  file.path(output_dir, "all_significant_de_features_fdr0.05_logfc1.csv")
)

print(de_qc)
