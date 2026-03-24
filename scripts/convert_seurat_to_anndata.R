#!/usr/bin/env Rscript
#
# Standalone Seurat v5 -> AnnData (.h5ad) converter
#
# Usage:
#   Rscript scripts/convert_seurat_to_anndata.R \
#     --input  /path/to/seurat.rds \
#     --output /path/to/output.h5ad \
#     [--assays RNA,SCT] \
#     [--python python3] \
#     [--keep-intermediate] \
#     [--no-graphs]
#
# The conversion works in two phases:
#   1. R extracts the Seurat object into an intermediate directory
#   2. Python builds the AnnData object from the intermediate directory

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx + 1 > length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

has_flag <- function(args, flag) flag %in% args

input_path  <- parse_arg(args, "--input")
output_path <- parse_arg(args, "--output")
assay_str   <- parse_arg(args, "--assays")
python_bin  <- parse_arg(args, "--python", "python3")
keep_inter  <- has_flag(args, "--keep-intermediate")
no_graphs   <- has_flag(args, "--no-graphs")

if (is.null(input_path) || is.null(output_path)) {
  cat("Seurat v5 -> AnnData converter\n\n")
  cat("Usage:\n")
  cat("  Rscript scripts/convert_seurat_to_anndata.R --input <rds> --output <h5ad>\n\n")
  cat("Options:\n")
  cat("  --input   PATH        Input Seurat RDS file (required)\n")
  cat("  --output  PATH        Output .h5ad file (required)\n")
  cat("  --assays  RNA,SCT     Assays to export (default: default assay only)\n")
  cat("  --python  PATH        Python binary (default: python3)\n")
  cat("  --keep-intermediate   Keep the intermediate directory after conversion\n")
  cat("  --no-graphs           Skip exporting graph objects\n")
  quit(status = 1)
}

if (!file.exists(input_path)) stop("Input file not found: ", input_path)

# Resolve script paths relative to this script's location
get_script_dir <- function() {
  # Try multiple methods to find script location
  for (i in seq_len(sys.nframe())) {
    ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  # Fallback: check if we're in the repo root
  if (file.exists("scripts/conversion_utils.R")) return("scripts")
  stop("Cannot determine script directory. Run from the project root.")
}
script_dir <- get_script_dir()
source(file.path(script_dir, "conversion_utils.R"))

# Parse assay list
assays <- if (!is.null(assay_str)) {
  trimws(strsplit(assay_str, ",")[[1]])
} else {
  NULL
}

# Intermediate directory beside the output file
inter_dir <- paste0(tools::file_path_sans_ext(output_path), "_intermediate")

cat("=== Seurat -> AnnData Conversion ===\n")
cat("Input:  ", input_path, "\n")
cat("Output: ", output_path, "\n\n")

# Phase 1: Load Seurat object and export to intermediate
cat("Phase 1: Loading Seurat object...\n")
library(Seurat)
obj <- readRDS(input_path)
cat("  Loaded:", ncol(obj), "cells,", nrow(obj), "features\n")
cat("  Assays:", paste(names(obj@assays), collapse = ", "), "\n")
cat("  Reductions:", paste(names(obj@reductions), collapse = ", "), "\n\n")

cat("Phase 1: Exporting to intermediate format...\n")
seurat_to_intermediate(
  obj, inter_dir,
  assays = assays,
  join_layers = TRUE,
  include_graphs = !no_graphs
)
cat("\n")

# Phase 2: Build AnnData via Python
cat("Phase 2: Building AnnData via Python...\n")
build_script <- file.path(script_dir, "build_anndata.py")
if (!file.exists(build_script)) {
  stop("Python helper not found: ", build_script)
}

# Ensure output directory exists
out_parent <- dirname(output_path)
if (!dir.exists(out_parent)) dir.create(out_parent, recursive = TRUE)

status <- system2(
  python_bin,
  args = c(build_script, "--input", inter_dir, "--output", output_path),
  stdout = "", stderr = ""
)

if (!identical(status, 0L)) {
  # Re-run with output visible for debugging
  system2(python_bin,
          args = c(build_script, "--input", inter_dir, "--output", output_path))
  stop("Python build_anndata.py failed with exit code ", status)
}

cat("\n")

# Cleanup
if (!keep_inter && dir.exists(inter_dir)) {
  unlink(inter_dir, recursive = TRUE)
  cat("Cleaned up intermediate directory\n")
} else if (keep_inter) {
  cat("Intermediate directory kept at:", inter_dir, "\n")
}

cat("\n=== Conversion complete ===\n")
cat("Output: ", output_path, "\n")
