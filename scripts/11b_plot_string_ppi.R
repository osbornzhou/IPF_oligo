#!/usr/bin/env Rscript

user_library <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", paste0(R.version$major, ".", R.version$minor))
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ppi_dir <- file.path(project_dir, "results", "ppi_network")
plot_dir <- file.path(ppi_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

plot_hub_bar <- function(nodes_path, out_path, title) {
  nodes <- fread(nodes_path)
  top <- nodes[order(hub_rank)][1:min(.N, 20)]
  top[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p <- ggplot(top, aes(x = hub_score, y = gene_symbol, fill = degree)) +
    geom_col(width = 0.75) +
    scale_fill_gradient(low = "#87a9c4", high = "#b42318") +
    labs(title = title, x = "Hub score", y = NULL, fill = "Degree") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(out_path, p, width = 7.2, height = 5.8, dpi = 180)
}

plot_network <- function(edges_path, nodes_path, out_path, title, top_n = 50) {
  edges <- fread(edges_path)
  nodes <- fread(nodes_path)
  top_nodes <- nodes[order(hub_rank)][1:min(.N, top_n), gene_symbol]
  sub_edges <- edges[preferredName_A %in% top_nodes & preferredName_B %in% top_nodes]
  if (nrow(sub_edges) == 0) {
    return(invisible(FALSE))
  }
  vertices <- unique(nodes[gene_symbol %in% top_nodes], by = "gene_symbol")
  vertices[, name := gene_symbol]
  setcolorder(vertices, "name")
  graph <- graph_from_data_frame(
    sub_edges[, .(from = preferredName_A, to = preferredName_B, weight = combined_score)],
    directed = FALSE,
    vertices = vertices
  )
  set.seed(20260524)
  layout <- layout_with_fr(graph, weights = E(graph)$weight)
  png(out_path, width = 1800, height = 1500, res = 180)
  plot(
    graph,
    layout = layout,
    vertex.size = 4 + 10 * (V(graph)$degree / max(V(graph)$degree, 1)),
    vertex.label = V(graph)$name,
    vertex.label.cex = 0.62,
    vertex.label.color = "#172033",
    vertex.color = ifelse(V(graph)$hub_top20, "#b42318", "#5b84a4"),
    vertex.frame.color = "white",
    edge.width = 0.5 + 2.5 * E(graph)$weight,
    edge.color = grDevices::adjustcolor("#6f7785", alpha.f = 0.45),
    main = title
  )
  dev.off()
  invisible(TRUE)
}

plot_hub_bar(
  file.path(ppi_dir, "string_ppi_nodes_robust_mrna_strict_medium_confidence.csv"),
  file.path(plot_dir, "robust_mrna_string_medium_top20_hub_barplot.png"),
  "STRING PPI hub genes - robust mRNA, medium confidence"
)

plot_hub_bar(
  file.path(ppi_dir, "string_ppi_nodes_robust_mrna_strict_high_confidence.csv"),
  file.path(plot_dir, "robust_mrna_string_high_top20_hub_barplot.png"),
  "STRING PPI hub genes - robust mRNA, high confidence"
)

plot_network(
  file.path(ppi_dir, "string_ppi_edges_robust_mrna_strict_medium_confidence.csv"),
  file.path(ppi_dir, "string_ppi_nodes_robust_mrna_strict_medium_confidence.csv"),
  file.path(plot_dir, "robust_mrna_string_medium_top50_network.png"),
  "STRING PPI top hub subnetwork - medium confidence",
  top_n = 50
)

plot_network(
  file.path(ppi_dir, "string_ppi_edges_robust_mrna_strict_high_confidence.csv"),
  file.path(ppi_dir, "string_ppi_nodes_robust_mrna_strict_high_confidence.csv"),
  file.path(plot_dir, "robust_mrna_string_high_top50_network.png"),
  "STRING PPI top hub subnetwork - high confidence",
  top_n = 50
)

message("PPI plots written to: ", plot_dir)
