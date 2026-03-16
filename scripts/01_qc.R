#!/usr/bin/env Rscript

# =============================================================================
# 01_qc.R  —  Quality Control
# Input : Cell Ranger filtered_feature_bc_matrix directories (per sample)
# Output: checkpoints/01_qc.rds  — named list of per-sample filtered Seurat objects
#         pipeline_outputs/plots/ — QC violin & scatter plots (before & after)
#         pipeline_outputs/qc_summary.csv
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

source("scripts/00_load_config.R")

cat("=== STEP 1: Quality Control ===\n\n")

# =============================================================================
# 1A  LOAD CELL RANGER DATA & CREATE SEURAT OBJECTS
# =============================================================================

cat("--- 1A: Loading Cell Ranger data ---\n")

seurat_objects <- list()

for (sample_name in names(cellranger_dirs)) {
  cat("Processing", sample_name, "...\n")
  cellranger_path <- cellranger_dirs[[sample_name]]

  if (!dir.exists(cellranger_path)) {
    cat("  Warning: Directory not found:", cellranger_path, "\n")
    next
  }

  tryCatch({
    count_matrix <- Read10X(cellranger_path)

    if (class(count_matrix) == "list") {
      cat("  Multiple assays detected — using Gene Expression\n")
      count_matrix <- count_matrix[["Gene Expression"]]
    }

    seurat_obj <- CreateSeuratObject(
      counts      = count_matrix,
      project     = sample_name,
      min.cells   = min_cells,
      min.features = min_features
    )

    # Attach sample + user-defined metadata
    seurat_obj@meta.data$sample <- sample_name
    for (field in names(sample_metadata[[sample_name]])) {
      seurat_obj@meta.data[[field]] <- sample_metadata[[sample_name]][[field]]
    }

    # Mitochondrial %
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)

    seurat_objects[[sample_name]] <- seurat_obj
    cat("  Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")

  }, error = function(e) {
    cat("  ERROR loading", sample_name, ":", e$message, "\n")
  })
}

cat("\nLoaded", length(seurat_objects), "samples\n\n")

# =============================================================================
# 1B  PRE-FILTER QC PLOTS
# =============================================================================

cat("--- 1B: Generating pre-filter QC plots ---\n")

plot_qc_violin <- function(obj_list, suffix) {
  plots <- lapply(names(obj_list), function(sn) {
    obj <- obj_list[[sn]]
    p1 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0) +
            ggtitle(paste(sn, "— nFeature")) + NoLegend()
    p2 <- VlnPlot(obj, features = "nCount_RNA",   pt.size = 0) +
            ggtitle(paste(sn, "— nCount"))   + NoLegend()
    p3 <- VlnPlot(obj, features = "percent.mt",   pt.size = 0) +
            ggtitle(paste(sn, "— MT%"))      + NoLegend()
    p1 | p2 | p3
  })
  wrap_plots(plots, ncol = 1)
}

plot_qc_scatter <- function(obj_list, suffix) {
  plots <- lapply(names(obj_list), function(sn) {
    obj <- obj_list[[sn]]
    p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
            ggtitle(paste(sn, "— Count vs Features")) + NoLegend()
    p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
            ggtitle(paste(sn, "— Count vs MT%")) + NoLegend()
    p1 | p2
  })
  wrap_plots(plots, ncol = 1)
}

pre_violin  <- plot_qc_violin(seurat_objects,  "pre")
pre_scatter <- plot_qc_scatter(seurat_objects, "pre")

ggsave(file.path(plots_dir, "qc_violin_pre_filter.pdf"),  pre_violin,
       width = 14, height = 4 * length(seurat_objects))
ggsave(file.path(plots_dir, "qc_scatter_pre_filter.pdf"), pre_scatter,
       width = 12, height = 4 * length(seurat_objects))

cat("Pre-filter QC plots saved\n")

# =============================================================================
# 1C  DOUBLET DETECTION (scDblFinder)
# =============================================================================

if (run_doublet_detection) {
  cat("\n--- 1C: Doublet detection (scDblFinder) ---\n")

  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("Package 'scDblFinder' not found. Install with:\n",
         "  BiocManager::install('scDblFinder')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' not found. Install with:\n",
         "  BiocManager::install('SingleCellExperiment')")
  }

  library(scDblFinder)
  library(SingleCellExperiment)

  for (sample_name in names(seurat_objects)) {
    cat("  Running scDblFinder on", sample_name, "...\n")
    tryCatch({
      sce <- as.SingleCellExperiment(seurat_objects[[sample_name]])
      sce <- scDblFinder(sce)
      seurat_objects[[sample_name]]@meta.data$doublet_score <- sce$scDblFinder.score
      seurat_objects[[sample_name]]@meta.data$doublet_class <- sce$scDblFinder.class
      n_doublets <- sum(sce$scDblFinder.class == "doublet")
      cat("    Detected", n_doublets, "doublets\n")
    }, error = function(e) {
      cat("    WARNING: scDblFinder failed for", sample_name, ":", e$message, "\n")
    })
  }
}

# =============================================================================
# 1D  APPLY FILTERS
# =============================================================================

cat("\n--- 1D: Applying QC filters ---\n")

qc_summary <- data.frame(
  Sample           = character(),
  Cells_Before     = integer(),
  Cells_After      = integer(),
  Cells_Removed    = integer(),
  Pct_Retained     = numeric(),
  Median_nFeature  = numeric(),
  Median_nCount    = numeric(),
  Median_MT_Pct    = numeric(),
  stringsAsFactors = FALSE
)

filtered_objects <- list()

for (sample_name in names(seurat_objects)) {
  obj    <- seurat_objects[[sample_name]]
  before <- ncol(obj)
  md     <- obj@meta.data

  # Build filter expression from config thresholds (NULL values disable that filter)
  keep <- rep(TRUE, ncol(obj))
  if (!is.null(qc_min_features))   keep <- keep & (md$nFeature_RNA >= qc_min_features)
  if (!is.null(qc_max_features))   keep <- keep & (md$nFeature_RNA <= qc_max_features)
  if (!is.null(qc_min_counts))     keep <- keep & (md$nCount_RNA   >= qc_min_counts)
  if (!is.null(qc_max_counts))     keep <- keep & (md$nCount_RNA   <= qc_max_counts)
  if (!is.null(qc_max_mt_percent)) keep <- keep & (md$percent.mt   <= qc_max_mt_percent)

  # Remove doublets if detection was run
  if (run_doublet_detection && "doublet_class" %in% colnames(md)) {
    keep <- keep & (md$doublet_class == "singlet")
  }

  cells_keep <- colnames(obj)[keep]
  obj_filt   <- subset(obj, cells = cells_keep)
  after      <- ncol(obj_filt)

  cat(sample_name, ": before =", before, "| after =", after,
      "| removed =", before - after, "\n")

  filtered_objects[[sample_name]] <- obj_filt

  qc_summary <- rbind(qc_summary, data.frame(
    Sample          = sample_name,
    Cells_Before    = before,
    Cells_After     = after,
    Cells_Removed   = before - after,
    Pct_Retained    = round(after / before * 100, 2),
    Median_nFeature = round(median(obj_filt$nFeature_RNA), 1),
    Median_nCount   = round(median(obj_filt$nCount_RNA),   1),
    Median_MT_Pct   = round(median(obj_filt$percent.mt),   2),
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# 1E  POST-FILTER QC PLOTS
# =============================================================================

cat("\n--- 1E: Generating post-filter QC plots ---\n")

post_violin  <- plot_qc_violin(filtered_objects,  "post")
post_scatter <- plot_qc_scatter(filtered_objects, "post")

ggsave(file.path(plots_dir, "qc_violin_post_filter.pdf"),  post_violin,
       width = 14, height = 4 * length(filtered_objects))
ggsave(file.path(plots_dir, "qc_scatter_post_filter.pdf"), post_scatter,
       width = 12, height = 4 * length(filtered_objects))

cat("Post-filter QC plots saved\n")

# =============================================================================
# 1F  SAVE QC SUMMARY & CHECKPOINT
# =============================================================================

summary_file <- file.path(output_dir, "qc_summary.csv")
write.csv(qc_summary, summary_file, row.names = FALSE)
cat("\nQC summary saved to:", summary_file, "\n")
print(qc_summary)

saveRDS(filtered_objects, ckpt$qc)
cat("\nCheckpoint saved:", ckpt$qc, "\n")
cat("\n=== STEP 1 COMPLETE ===\n")
