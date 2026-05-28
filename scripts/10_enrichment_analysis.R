#!/usr/bin/env Rscript

user_library <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", paste0(R.version$major, ".", R.version$minor))
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(AnnotationDbi)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
robust_dir <- file.path(project_dir, "results", "robust_candidates")
axes_dir <- file.path(project_dir, "results", "mirna_mrna_axes")
de_annotated_dir <- file.path(project_dir, "results", "differential_expression_annotated")
output_dir <- file.path(project_dir, "results", "enrichment")
plot_dir <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
min_gene_set_size <- 5

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
}

mrna_path <- file.path(robust_dir, "robust_mrna_candidates_strict.csv")
axes_path <- file.path(axes_dir, "robust_mirna_mrna_negative_axes_mirtarbase.csv")
background_path <- file.path(de_annotated_dir, "GSE32537_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv")

stop_if_missing(mrna_path)
stop_if_missing(axes_path)
stop_if_missing(background_path)

robust_mrna <- fread(mrna_path, na.strings = c("", "NA", "N/A"))
axes <- fread(axes_path, na.strings = c("", "NA", "N/A"))
background_de <- fread(background_path, na.strings = c("", "NA", "N/A"))

normalize_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | toupper(x) == "NA"] <- NA_character_
  unique(na.omit(toupper(x)))
}

background_symbols <- normalize_symbol(background_de[is_annotated %in% c(TRUE, "True", "TRUE", "true", 1), standard_feature_id])
robust_symbols <- normalize_symbol(robust_mrna$standard_feature_id)
robust_up <- normalize_symbol(robust_mrna[discovery_logFC > 0, standard_feature_id])
robust_down <- normalize_symbol(robust_mrna[discovery_logFC < 0, standard_feature_id])
axis_targets <- normalize_symbol(axes$target_gene)
axis_targets_exact <- normalize_symbol(axes[match_type == "exact", target_gene])

gene_sets <- list(
  robust_mrna_strict = robust_symbols,
  robust_mrna_up = robust_up,
  robust_mrna_down = robust_down,
  mirna_mrna_axis_targets_all = axis_targets,
  mirna_mrna_axis_targets_exact = axis_targets_exact
)

all_gene_set_sizes <- data.table(
  gene_set = names(gene_sets),
  input_symbols = vapply(gene_sets, length, integer(1)),
  included_for_enrichment = vapply(gene_sets, length, integer(1)) >= min_gene_set_size,
  minimum_required_symbols = min_gene_set_size
)
fwrite(all_gene_set_sizes, file.path(output_dir, "gene_set_inclusion_qc.csv"))

gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_gene_set_size]

map_symbols_to_entrez <- function(symbols, set_name) {
  symbols <- normalize_symbol(symbols)
  if (length(symbols) == 0) {
    return(data.table())
  }
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID", "GENENAME")
  )
  setDT(mapped)
  mapped[, gene_set := set_name]
  mapped[, input_symbol := SYMBOL]
  mapped <- unique(mapped[, .(gene_set, input_symbol, SYMBOL, ENTREZID, GENENAME)])
  mapped
}

all_mapping <- rbindlist(
  c(
    list(map_symbols_to_entrez(background_symbols, "background_gse32537_all_tested_annotated_mrna")),
    lapply(names(gene_sets), function(nm) map_symbols_to_entrez(gene_sets[[nm]], nm))
  ),
  fill = TRUE
)

fwrite(all_mapping, file.path(output_dir, "gene_id_mapping.csv"))

background_entrez <- unique(na.omit(all_mapping[
  gene_set == "background_gse32537_all_tested_annotated_mrna",
  ENTREZID
]))

gene_set_table <- rbindlist(lapply(names(gene_sets), function(nm) {
  data.table(gene_set = nm, gene_symbol = gene_sets[[nm]])
}))
fwrite(gene_set_table, file.path(output_dir, "gene_sets_used.csv"))

gene_ratio_numeric <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 2) return(NA_real_)
    as.numeric(p[1]) / as.numeric(p[2])
  }, numeric(1))
}

save_plot_pair <- function(res_dt, gene_set_name, database) {
  if (nrow(res_dt) == 0) {
    return(invisible(FALSE))
  }
  plot_dt <- copy(res_dt[order(p.adjust)][seq_len(min(.N, 15))])
  plot_dt[, GeneRatioNumeric := gene_ratio_numeric(GeneRatio)]
  plot_dt[, DescriptionWrapped := stringr::str_wrap(Description, width = 48)]
  plot_dt[, DescriptionWrapped := factor(DescriptionWrapped, levels = rev(DescriptionWrapped))]

  dot <- ggplot(plot_dt, aes(x = GeneRatioNumeric, y = DescriptionWrapped)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
    scale_color_gradient(low = "#b42318", high = "#1f4e79", trans = "reverse") +
    labs(
      title = paste(gene_set_name, database, sep = " - "),
      x = "Gene ratio",
      y = NULL,
      size = "Count",
      color = "Adjusted P"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  bar <- ggplot(plot_dt, aes(x = Count, y = DescriptionWrapped, fill = p.adjust)) +
    geom_col(width = 0.75) +
    scale_fill_gradient(low = "#b42318", high = "#1f4e79", trans = "reverse") +
    labs(
      title = paste(gene_set_name, database, sep = " - "),
      x = "Gene count",
      y = NULL,
      fill = "Adjusted P"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  safe_name <- gsub("[^A-Za-z0-9_]+", "_", paste(gene_set_name, database, sep = "_"))
  ggsave(file.path(plot_dir, paste0(safe_name, "_dotplot.png")), dot, width = 8.2, height = 5.6, dpi = 180)
  ggsave(file.path(plot_dir, paste0(safe_name, "_barplot.png")), bar, width = 8.2, height = 5.6, dpi = 180)
  invisible(TRUE)
}

run_enrichment <- function(entrez_ids, gene_set_name, database, fun) {
  if (length(entrez_ids) < min_gene_set_size) {
    return(data.table())
  }
  result <- tryCatch(fun(entrez_ids), error = function(e) e)
  if (inherits(result, "error") || is.null(result)) {
    warning("Enrichment failed for ", gene_set_name, " / ", database, ": ", conditionMessage(result))
    return(data.table())
  }
  res <- as.data.table(as.data.frame(result))
  if (nrow(res) == 0) {
    return(data.table())
  }
  res[, gene_set := gene_set_name]
  res[, database := database]
  res[, input_entrez_count := length(entrez_ids)]
  setcolorder(res, c("gene_set", "database", "input_entrez_count", setdiff(names(res), c("gene_set", "database", "input_entrez_count"))))
  res
}

enrichment_results <- list()
qc_rows <- list()

for (gene_set_name in names(gene_sets)) {
  symbols <- gene_sets[[gene_set_name]]
  mapping <- all_mapping[gene_set == gene_set_name]
  mapped_entrez <- unique(na.omit(mapping$ENTREZID))
  mapped_symbols <- unique(na.omit(mapping[!is.na(ENTREZID), input_symbol]))

  qc_rows[[gene_set_name]] <- data.table(
    gene_set = gene_set_name,
    input_symbols = length(symbols),
    mapped_symbols = length(mapped_symbols),
    mapped_entrez_ids = length(mapped_entrez),
    unmapped_symbols = length(setdiff(symbols, mapped_symbols)),
    mapping_rate = round(length(mapped_symbols) / max(length(symbols), 1), 6),
    background_symbols = length(background_symbols),
    background_entrez_ids = length(background_entrez),
    qc1_input_pass = length(symbols) >= min_gene_set_size && length(intersect(symbols, background_symbols)) >= min_gene_set_size,
    qc2_mapping_pass = length(mapped_entrez) >= min_gene_set_size && length(mapped_symbols) / max(length(symbols), 1) >= 0.80
  )

  go_bp <- run_enrichment(
    mapped_entrez, gene_set_name, "GO_BP",
    function(ids) enrichGO(
      gene = ids, universe = background_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    )
  )
  go_cc <- run_enrichment(
    mapped_entrez, gene_set_name, "GO_CC",
    function(ids) enrichGO(
      gene = ids, universe = background_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
      ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    )
  )
  go_mf <- run_enrichment(
    mapped_entrez, gene_set_name, "GO_MF",
    function(ids) enrichGO(
      gene = ids, universe = background_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
      ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    )
  )
  kegg <- run_enrichment(
    mapped_entrez, gene_set_name, "KEGG",
    function(ids) enrichKEGG(
      gene = ids, universe = background_entrez, organism = "hsa",
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
    )
  )
  reactome <- run_enrichment(
    mapped_entrez, gene_set_name, "Reactome",
    function(ids) enrichPathway(
      gene = ids, universe = background_entrez, organism = "human",
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    )
  )

  set_results <- rbindlist(list(go_bp, go_cc, go_mf, kegg, reactome), fill = TRUE)
  if (nrow(set_results) > 0) {
    set_results[, significant := p.adjust < fdr_cutoff]
    fwrite(set_results, file.path(output_dir, paste0(gene_set_name, "_enrichment_all_databases.csv")))
    fwrite(set_results[p.adjust < fdr_cutoff], file.path(output_dir, paste0(gene_set_name, "_enrichment_significant_fdr0.05.csv")))
    for (db in unique(set_results$database)) {
      db_dt <- set_results[database == db & p.adjust < fdr_cutoff]
      save_plot_pair(db_dt, gene_set_name, db)
    }
  } else {
    fwrite(data.table(), file.path(output_dir, paste0(gene_set_name, "_enrichment_all_databases.csv")))
    fwrite(data.table(), file.path(output_dir, paste0(gene_set_name, "_enrichment_significant_fdr0.05.csv")))
  }
  enrichment_results[[gene_set_name]] <- set_results
}

all_results <- rbindlist(enrichment_results, fill = TRUE)
if (nrow(all_results) > 0) {
  fwrite(all_results, file.path(output_dir, "all_enrichment_results.csv"))
  fwrite(all_results[p.adjust < fdr_cutoff], file.path(output_dir, "all_enrichment_significant_fdr0.05.csv"))
} else {
  fwrite(data.table(), file.path(output_dir, "all_enrichment_results.csv"))
  fwrite(data.table(), file.path(output_dir, "all_enrichment_significant_fdr0.05.csv"))
}

qc <- rbindlist(qc_rows, fill = TRUE)
if (nrow(all_results) > 0) {
  db_qc <- all_results[, .(
    tested_terms = .N,
    significant_terms = sum(p.adjust < fdr_cutoff, na.rm = TRUE),
    min_adjusted_p = min(p.adjust, na.rm = TRUE)
  ), by = .(gene_set, database)]
} else {
  db_qc <- data.table(gene_set = character(), database = character(), tested_terms = integer(), significant_terms = integer(), min_adjusted_p = numeric())
}

qc <- merge(qc, db_qc[, .(
  total_tested_terms = sum(tested_terms),
  total_significant_terms = sum(significant_terms)
), by = gene_set], by = "gene_set", all.x = TRUE)
qc[is.na(total_tested_terms), total_tested_terms := 0]
qc[is.na(total_significant_terms), total_significant_terms := 0]
qc[, qc3_enrichment_pass := total_tested_terms > 0 & !is.na(total_significant_terms)]
qc[, triple_qc_pass := qc1_input_pass & qc2_mapping_pass & qc3_enrichment_pass]
qc[, organism := "Homo sapiens"]
qc[, background_definition := "All annotated gene-level features tested in GSE32537 discovery mRNA differential expression"]
qc[, ontology_sources := paste(
  paste0("GO.db ", as.character(packageVersion("GO.db"))),
  paste0("org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db"))),
  paste0("clusterProfiler ", as.character(packageVersion("clusterProfiler"))),
  paste0("ReactomePA ", as.character(packageVersion("ReactomePA"))),
  paste0("reactome.db ", as.character(packageVersion("reactome.db"))),
  paste0("KEGGREST ", as.character(packageVersion("KEGGREST")), " via KEGG REST"),
  sep = "; "
)]

fwrite(qc, file.path(output_dir, "enrichment_triple_qc.csv"))
fwrite(db_qc, file.path(output_dir, "enrichment_database_qc.csv"))

message("CSV and plot outputs are complete. Run scripts/10b_export_enrichment_excel.py to build the verified Excel workbook.")

print(qc)
print(db_qc)
