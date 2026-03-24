#!/usr/bin/env Rscript

# =============================================================================
# Step 10: Export to AnnData (.h5ad)
#
# Optional pipeline step that converts the latest Seurat checkpoint to .h5ad.
# Uses the intermediate directory approach for robust Seurat v5 conversion.
# =============================================================================

cat("Step 10: Export to AnnData\n\n")

source("scripts/conversion_utils.R")

# Load the best available upstream checkpoint
obj <- read_best_checkpoint(
  c(ckpt$downstream, ckpt$de, ckpt$annotation, ckpt$clustering,
    ckpt$integration, ckpt$norm_pca),
  missing_message = "No upstream checkpoint found for export."
)

cat("Loaded checkpoint:", ncol(obj), "cells,", nrow(obj), "features\n")
cat("Assays:", paste(names(obj@assays), collapse = ", "), "\n")
cat("Reductions:", paste(names(obj@reductions), collapse = ", "), "\n\n")

# Determine which assays to export
export_assays <- if (!is.null(export_cfg$assays)) {
  export_cfg$assays
} else {
  DefaultAssay(obj)
}

# Output path
h5ad_path <- file.path(exports_dir, paste0(project_name, ".h5ad"))
inter_dir <- file.path(exports_dir, ".seurat_intermediate")

# Phase 1: Seurat -> intermediate
cat("Exporting Seurat object to intermediate format...\n")
seurat_to_intermediate(
  obj, inter_dir,
  assays = export_assays,
  join_layers = isTRUE(export_cfg$join_layers %||% TRUE),
  include_graphs = isTRUE(export_cfg$include_graphs %||% TRUE)
)
cat("\n")

# Phase 2: intermediate -> AnnData via Python
cat("Building AnnData via Python...\n")
build_script <- "scripts/build_anndata.py"
python_bin <- export_cfg$python_bin %||% "python3"

status <- system2(
  python_bin,
  args = c(build_script, "--input", inter_dir, "--output", h5ad_path),
  stdout = "", stderr = ""
)

if (!identical(status, 0L)) {
  system2(python_bin,
          args = c(build_script, "--input", inter_dir, "--output", h5ad_path))
  stop("build_anndata.py failed with exit code ", status)
}

# Cleanup intermediate
unlink(inter_dir, recursive = TRUE)

cat("\nExport complete:", h5ad_path, "\n")
