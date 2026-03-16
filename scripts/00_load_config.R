#!/usr/bin/env Rscript

# =============================================================================
# 00_load_config.R
# Reads config/config.yaml and exposes all variables used by pipeline modules.
# Source this at the top of every module.
# =============================================================================

if (!requireNamespace("yaml", quietly = TRUE)) stop("Install 'yaml': install.packages('yaml')")
library(yaml)

CONFIG_PATH <- Sys.getenv("SCRNA_CONFIG", unset = "config.yaml")
if (!file.exists(CONFIG_PATH)) {
  stop("Config file not found: ", CONFIG_PATH,
       "\nRun from the project root, or set SCRNA_CONFIG env variable.")
}

cfg <- yaml::read_yaml(CONFIG_PATH)

# -----------------------------------------------------------------------------
# Project
# -----------------------------------------------------------------------------
project_name <- cfg$project$name

# -----------------------------------------------------------------------------
# Step toggles
# -----------------------------------------------------------------------------
run_cellbender  <- isTRUE(cfg$steps$run_cellbender)
run_cell_cycle  <- isTRUE(cfg$steps$run_cell_cycle)
run_integration <- isTRUE(cfg$steps$run_integration)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
cellbender_dir  <- cfg$paths$cellbender_output_dir
output_dir      <- cfg$paths$output_dir
checkpoints_dir <- cfg$paths$checkpoints_dir
plots_dir       <- file.path(output_dir, "plots")

# Create dirs if needed
for (d in c(output_dir, checkpoints_dir, plots_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Checkpoint file paths — each step reads and writes one of these
ckpt <- list(
  qc          = file.path(checkpoints_dir, "01_qc.rds"),
  cellbender  = file.path(checkpoints_dir, "02_cellbender.rds"),
  cell_cycle  = file.path(checkpoints_dir, "03_cell_cycle.rds"),
  norm_pca    = file.path(checkpoints_dir, "04_norm_pca.rds"),
  integration = file.path(checkpoints_dir, "05_integration.rds"),
  clustering  = file.path(checkpoints_dir, "06_clustering.rds")
)

# -----------------------------------------------------------------------------
# Samples
# -----------------------------------------------------------------------------
sample_names <- sapply(cfg$samples, `[[`, "name")

cellranger_dirs <- setNames(
  lapply(cfg$samples, `[[`, "cellranger_path"),
  sample_names
)

cb_files <- setNames(
  lapply(cfg$samples, function(s) file.path(cellbender_dir, s$cellbender_h5)),
  sample_names
)

barcode_files <- setNames(
  lapply(cfg$samples, function(s) file.path(cellbender_dir, s$cellbender_csv)),
  sample_names
)

sample_metadata <- setNames(
  lapply(cfg$samples, `[[`, "metadata"),
  sample_names
)

# All unique metadata field names across samples
all_meta_fields <- unique(unlist(lapply(sample_metadata, names)))

# -----------------------------------------------------------------------------
# Seurat creation
# -----------------------------------------------------------------------------
min_cells    <- cfg$seurat$min_cells
min_features <- cfg$seurat$min_features
mt_pattern   <- cfg$seurat$mt_pattern

# -----------------------------------------------------------------------------
# QC thresholds
# -----------------------------------------------------------------------------
qc_min_features        <- cfg$qc$min_features
qc_max_features        <- cfg$qc$max_features
qc_min_counts          <- cfg$qc$min_counts
qc_max_counts          <- cfg$qc$max_counts
qc_max_mt_percent      <- cfg$qc$max_mt_percent
run_doublet_detection  <- isTRUE(cfg$qc$run_doublet_detection)

# -----------------------------------------------------------------------------
# Normalization & PCA
# -----------------------------------------------------------------------------
n_variable_features <- cfg$analysis$n_variable_features
n_pca_components    <- cfg$analysis$n_pca_components

# -----------------------------------------------------------------------------
# Integration
# -----------------------------------------------------------------------------
harmony_group_vars <- cfg$integration$harmony_group_vars

# -----------------------------------------------------------------------------
# Clustering & UMAP
# -----------------------------------------------------------------------------
clustering_reduction   <- cfg$clustering$reduction
clustering_resolutions <- cfg$clustering$resolutions
selected_resolution    <- cfg$clustering$selected_resolution

cat("Config loaded:", CONFIG_PATH, "\n")
cat("Project:", project_name, "\n")
cat("Samples:", paste(sample_names, collapse = ", "), "\n")
cat("Steps  — CellBender:", run_cellbender,
    "| Cell Cycle:", run_cell_cycle,
    "| Integration:", run_integration, "\n\n")
