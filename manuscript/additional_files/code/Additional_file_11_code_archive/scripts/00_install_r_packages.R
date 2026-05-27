options(repos = c(CRAN = "https://cloud.r-project.org"))

user_library <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", paste0(R.version$major, ".", R.version$minor))
dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_library, .libPaths()))

message("Using R library: ", user_library)

cran_packages <- c(
  "data.table",
  "readr",
  "dplyr",
  "tibble",
  "ggplot2",
  "pheatmap",
  "matrixStats",
  "openxlsx"
)

bioc_packages <- c(
  "limma",
  "edgeR",
  "sva",
  "Biobase",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "enrichplot",
  "DOSE"
)

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, dependencies = TRUE)
  }
}

install_missing_cran("BiocManager")
install_missing_cran(cran_packages)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

all_packages <- c(cran_packages, bioc_packages)
status <- data.frame(
  package = all_packages,
  installed = vapply(all_packages, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)

print(status)

if (!all(status$installed)) {
  stop("Some required R packages failed to install.")
}
