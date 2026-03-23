#!/usr/bin/env Rscript

cran_packages <- c(
  "yaml", "dplyr", "ggplot2", "patchwork", "RColorBrewer", "Matrix",
  "readr", "tidyr", "viridis", "hdf5r", "Seurat", "harmony"
)

bioc_packages <- c(
  "scDblFinder", "SingleCellExperiment", "slingshot", "miloR", "edgeR", "zellkonverter"
)

optional_manual <- c("CHOIR", "CellChat", "monocle3", "SeuratWrappers", "CytoTRACE2")

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) installer(missing)
}

cat("Installing CRAN packages...\n")
install_if_missing(cran_packages, function(pkgs) install.packages(pkgs, repos = "https://cloud.r-project.org"))

cat("\nInstalling Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_if_missing(bioc_packages, function(pkgs) BiocManager::install(pkgs))

cat("\nOptional packages not installed automatically here:\n")
for (pkg in optional_manual) cat(" - ", pkg, "\n", sep = "")
cat("Install those separately in the environment where you plan to use the corresponding optional stages.\n")
