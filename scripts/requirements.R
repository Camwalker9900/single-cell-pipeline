#!/usr/bin/env Rscript

# =============================================================================
# envs/requirements.R
# Install all R packages required by the pipeline.
# Run once before first use: Rscript envs/requirements.R
# =============================================================================

cran_packages <- c("yaml", "dplyr", "ggplot2", "patchwork", "RColorBrewer", "Matrix")

bioc_packages <- c("Seurat", "harmony", "scDblFinder", "SingleCellExperiment")

# hdf5r for reading .h5 files (CellBender step)
hdf5r_packages <- c("hdf5r")

cat("Installing CRAN packages...\n")
install.packages(c(cran_packages, hdf5r_packages), repos = "https://cloud.r-project.org")

cat("\nInstalling Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(bioc_packages)

cat("\nAll packages installed.\n")
