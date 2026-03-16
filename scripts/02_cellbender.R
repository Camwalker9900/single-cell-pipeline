#!/usr/bin/env Rscript

# =============================================================================
# 02_cellbender.R  —  CellBender Ambient RNA Filtering  (OPTIONAL)
# Input : checkpoints/01_qc.rds  — named list of QC-filtered Seurat objects
# Output: checkpoints/02_cellbender.rds  — named list, further filtered by CellBender
#         pipeline_outputs/cellbender_filtering_summary.csv
#
# If steps.run_cellbender is false in config, this script is skipped and the
# orchestrator passes checkpoints/01_qc.rds directly to the next step.
# =============================================================================

library(Seurat)
library(hdf5r)
library(dplyr)
library(Matrix)

source("scripts/00_load_config.R")

cat("=== STEP 2: CellBender Ambient RNA Filtering ===\n\n")

# =============================================================================
# LOAD PREVIOUS CHECKPOINT
# =============================================================================

if (!file.exists(ckpt$qc)) stop("QC checkpoint not found: ", ckpt$qc,
                                 "\nRun 01_qc.R first.")
seurat_objects <- readRDS(ckpt$qc)
cat("Loaded", length(seurat_objects), "QC-filtered objects from checkpoint\n\n")

# =============================================================================
# READ CELLBENDER FILTERED BARCODES
# =============================================================================

cat("--- 2A: Reading CellBender filtered barcodes ---\n")

read_cellbender_barcodes <- function(h5_file) {
  tryCatch({
    h5_data <- Read10X_h5(h5_file)
    if (class(h5_data) == "list") {
      barcodes <- colnames(h5_data[[1]])
    } else {
      barcodes <- colnames(h5_data)
    }
    return(barcodes)
  }, error = function(e) {
    cat("    Error reading", basename(h5_file), ":", e$message, "\n")
    return(NULL)
  })
}

read_barcode_csv <- function(csv_file) {
  if (file.exists(csv_file)) {
    barcodes <- read.csv(csv_file, header = FALSE, stringsAsFactors = FALSE)$V1
    return(as.character(barcodes))
  }
  return(NULL)
}

filtered_barcodes_list <- list()

for (sample in names(cb_files)) {
  cat("Reading CellBender results for", sample, "...\n")

  barcodes <- read_cellbender_barcodes(cb_files[[sample]])

  if (is.null(barcodes)) {
    cat("  .h5 failed, trying CSV...\n")
    barcodes <- read_barcode_csv(barcode_files[[sample]])
  }

  if (!is.null(barcodes)) {
    filtered_barcodes_list[[sample]] <- barcodes
    cat("  Found", length(barcodes), "filtered barcodes\n")
  } else {
    cat("  Warning: Could not read barcodes for", sample, "\n")
  }
}

cat("\n")

# =============================================================================
# FILTER SEURAT OBJECTS
# =============================================================================

cat("--- 2B: Filtering objects using CellBender barcodes ---\n")

filtered_seurat_objects <- list()
filtering_summary <- data.frame(
  Sample         = character(),
  Original_Cells = integer(),
  CellBender_Cells = integer(),
  Final_Cells    = integer(),
  Cells_Removed  = integer(),
  Retention_Rate = numeric(),
  stringsAsFactors = FALSE
)

for (sample_name in names(seurat_objects)) {
  cat("Filtering", sample_name, "...\n")

  seurat_obj     <- seurat_objects[[sample_name]]
  original_cells <- ncol(seurat_obj)

  if (!sample_name %in% names(filtered_barcodes_list)) {
    cat("  No CellBender results found — keeping all cells\n\n")
    filtered_seurat_objects[[sample_name]] <- seurat_obj
    next
  }

  cellbender_barcodes <- filtered_barcodes_list[[sample_name]]
  seurat_barcodes     <- colnames(seurat_obj)
  common_barcodes     <- intersect(seurat_barcodes, cellbender_barcodes)

  cat("  Original Seurat cells:", original_cells, "\n")
  cat("  CellBender filtered cells:", length(cellbender_barcodes), "\n")
  cat("  Common cells:", length(common_barcodes), "\n")

  if (length(common_barcodes) == 0) {
    cat("  ERROR: No common barcodes found!\n")
    cat("  Sample Seurat barcodes:", head(seurat_barcodes, 3), "\n")
    cat("  Sample CellBender barcodes:", head(cellbender_barcodes, 3), "\n")
    cat("  Skipping this sample...\n\n")
    next
  }

  filtered_obj <- subset(seurat_obj, cells = common_barcodes)

  filtered_obj@meta.data$cellbender_filtered <- TRUE
  filtered_obj@misc$filtering_info <- list(
    original_cells  = original_cells,
    cellbender_cells = length(cellbender_barcodes),
    final_cells     = ncol(filtered_obj),
    filtering_date  = Sys.Date()
  )

  filtered_seurat_objects[[sample_name]] <- filtered_obj

  final_cells    <- ncol(filtered_obj)
  cells_removed  <- original_cells - final_cells
  retention_rate <- round((final_cells / original_cells) * 100, 2)

  filtering_summary <- rbind(filtering_summary, data.frame(
    Sample           = sample_name,
    Original_Cells   = original_cells,
    CellBender_Cells = length(cellbender_barcodes),
    Final_Cells      = final_cells,
    Cells_Removed    = cells_removed,
    Retention_Rate   = retention_rate,
    stringsAsFactors = FALSE
  ))

  cat("  Final cells:", final_cells, "(", retention_rate, "% retained)\n\n")
}

# =============================================================================
# SAVE SUMMARY & CHECKPOINT
# =============================================================================

summary_file <- file.path(output_dir, "cellbender_filtering_summary.csv")
write.csv(filtering_summary, summary_file, row.names = FALSE)
cat("Filtering summary saved to:", summary_file, "\n")
print(filtering_summary)

saveRDS(filtered_seurat_objects, ckpt$cellbender)
cat("\nCheckpoint saved:", ckpt$cellbender, "\n")
cat("\n=== STEP 2 COMPLETE ===\n")
