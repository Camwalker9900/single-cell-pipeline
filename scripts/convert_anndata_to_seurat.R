#!/usr/bin/env Rscript
#
# Standalone AnnData (.h5ad) -> Seurat v5 converter
#
# Usage:
#   Rscript scripts/convert_anndata_to_seurat.R \
#     --input  /path/to/input.h5ad \
#     --output /path/to/seurat.rds \
#     [--python python3] \
#     [--keep-intermediate]
#
# The conversion works in two phases:
#   1. Python exports the AnnData object to an intermediate directory
#   2. R rebuilds the Seurat v5 object from the intermediate directory

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
python_bin  <- parse_arg(args, "--python", "python3")
keep_inter  <- has_flag(args, "--keep-intermediate")

if (is.null(input_path) || is.null(output_path)) {
  cat("AnnData (.h5ad) -> Seurat v5 converter\n\n")
  cat("Usage:\n")
  cat("  Rscript scripts/convert_anndata_to_seurat.R --input <h5ad> --output <rds>\n\n")
  cat("Options:\n")
  cat("  --input   PATH        Input .h5ad file (required)\n")
  cat("  --output  PATH        Output Seurat RDS file (required)\n")
  cat("  --python  PATH        Python binary (default: python3)\n")
  cat("  --keep-intermediate   Keep the intermediate directory after conversion\n")
  quit(status = 1)
}

if (!file.exists(input_path)) stop("Input file not found: ", input_path)

# Resolve script paths
`%||%` <- function(x, y) if (is.null(x)) y else x
get_script_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  if (file.exists("scripts/conversion_utils.R")) return("scripts")
  stop("Cannot determine script directory. Run from the project root.")
}
script_dir <- get_script_dir()
source(file.path(script_dir, "conversion_utils.R"))

# Intermediate directory beside the output file
inter_dir <- paste0(tools::file_path_sans_ext(output_path), "_intermediate")

cat("=== AnnData -> Seurat Conversion ===\n")
cat("Input:  ", input_path, "\n")
cat("Output: ", output_path, "\n\n")

# Phase 1: Export AnnData to intermediate via Python
cat("Phase 1: Exporting AnnData via Python...\n")
export_script <- file.path(script_dir, "export_anndata.py")
if (!file.exists(export_script)) {
  stop("Python helper not found: ", export_script)
}

status <- system2(
  python_bin,
  args = c(export_script, "--input", input_path, "--output", inter_dir),
  stdout = "", stderr = ""
)

if (!identical(status, 0L)) {
  system2(python_bin,
          args = c(export_script, "--input", input_path, "--output", inter_dir))
  stop("Python export_anndata.py failed with exit code ", status)
}
cat("\n")

# Phase 2: Build Seurat v5 object from intermediate
cat("Phase 2: Building Seurat v5 object...\n")
library(Seurat)
obj <- intermediate_to_seurat(inter_dir)
cat("\n")

# Save
out_parent <- dirname(output_path)
if (!dir.exists(out_parent)) dir.create(out_parent, recursive = TRUE)

saveRDS(obj, output_path)
cat("Saved Seurat object:", output_path, "\n")

# Cleanup
if (!keep_inter && dir.exists(inter_dir)) {
  unlink(inter_dir, recursive = TRUE)
  cat("Cleaned up intermediate directory\n")
} else if (keep_inter) {
  cat("Intermediate directory kept at:", inter_dir, "\n")
}

cat("\n=== Conversion complete ===\n")
cat("Output: ", output_path, "\n")
