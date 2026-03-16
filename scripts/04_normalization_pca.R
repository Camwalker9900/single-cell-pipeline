#!/usr/bin/env Rscript

# =============================================================================
# 04_normalization_pca.R  —  Normalization, HVG Selection, Scaling & PCA
# Input : checkpoint from previous step (cell_cycle, cellbender, or qc)
#         — combined Seurat object OR named list of per-sample objects
# Output: checkpoints/04_norm_pca.rds  — normalized, scaled, PCA-reduced object
#         pipeline_outputs/plots/variable_features.pdf
#         pipeline_outputs/plots/pca_*.pdf
#
# NOTE: If cell cycle step was run, normalization and HVG have already been
#       done. This script detects that and skips those sub-steps to avoid
#       running them twice.
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

source("scripts/00_load_config.R")

cat("=== STEP 4: Normalization & PCA ===\n\n")

# =============================================================================
# LOAD BEST AVAILABLE UPSTREAM CHECKPOINT
# Priority: cell_cycle > cellbender > qc
# =============================================================================

prev_ckpt <- if      (file.exists(ckpt$cell_cycle))  ckpt$cell_cycle  else
             if      (file.exists(ckpt$cellbender))   ckpt$cellbender  else
             if      (file.exists(ckpt$qc))            ckpt$qc          else
             stop("No upstream checkpoint found. Run 01_qc.R first.")

cat("Loading checkpoint:", prev_ckpt, "\n")
input_data <- readRDS(prev_ckpt)

# =============================================================================
# MERGE LIST INTO COMBINED OBJECT (if coming straight from qc or cellbender)
# =============================================================================

if (is.list(input_data) && !inherits(input_data, "Seurat")) {
  cat("Merging", length(input_data), "samples...\n")

  sample_list <- input_data
  for (i in seq_along(sample_list)) {
    sn  <- names(sample_list)[i]
    obj <- sample_list[[i]]
    obj <- RenameCells(obj, new.names = paste0(sn, "_", colnames(obj)))
    sample_list[[i]] <- obj
    cat("  Renamed", ncol(obj), "cells for", sn, "\n")
  }

  combined_obj <- merge(sample_list[[1]],
                        y       = sample_list[-1],
                        add.cell.ids = NULL,
                        project = project_name)

  cat("Combined:", ncol(combined_obj), "cells,", nrow(combined_obj), "genes\n\n")
} else {
  combined_obj <- input_data
  cat("Combined object:", ncol(combined_obj), "cells,", nrow(combined_obj), "genes\n\n")
}

# =============================================================================
# NORMALIZATION & VARIABLE FEATURES
# Skip if cell cycle step already ran these
# =============================================================================

already_normalized <- prev_ckpt == ckpt$cell_cycle

if (!already_normalized) {
  cat("--- 4A: Normalizing data ---\n")
  combined_obj <- NormalizeData(combined_obj,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)

  cat("--- 4B: Finding variable features (top", n_variable_features, ") ---\n")
  combined_obj <- FindVariableFeatures(combined_obj,
                                       selection.method = "vst",
                                       nfeatures = n_variable_features)
} else {
  cat("--- 4A/4B: Normalization & HVG already done in Step 3 — skipping ---\n")
}

# Variable features plot
top10_vars <- head(VariableFeatures(combined_obj), 10)
cat("Top 10 variable genes:", paste(top10_vars, collapse = ", "), "\n")

var_plot <- VariableFeaturePlot(combined_obj)
ggsave(file.path(plots_dir, "variable_features.pdf"), var_plot, width = 12, height = 8)
cat("Variable features plot saved\n")

# =============================================================================
# SCALING (only if cell cycle step was NOT run — that step already scales)
# =============================================================================

if (!already_normalized) {
  cat("\n--- 4C: Scaling data ---\n")
  all_genes    <- rownames(combined_obj)
  combined_obj <- ScaleData(combined_obj, features = all_genes)
  cat("Scaling complete\n")
} else {
  cat("\n--- 4C: Data already scaled with CC regression in Step 3 — skipping ---\n")
}

# =============================================================================
# PCA
# =============================================================================

cat("\n--- 4D: PCA (", n_pca_components, "components) ---\n")
combined_obj <- RunPCA(combined_obj,
                       features = VariableFeatures(object = combined_obj),
                       npcs     = n_pca_components)

# PCA plots
cc_group_fields <- head(all_meta_fields, 2)

pca_plot1 <- DimPlot(combined_obj, reduction = "pca", group.by = cc_group_fields[1]) +
  ggtitle(paste("PCA by", tools::toTitleCase(cc_group_fields[1])))

if (length(cc_group_fields) >= 2) {
  pca_plot2 <- DimPlot(combined_obj, reduction = "pca", group.by = cc_group_fields[2]) +
    ggtitle(paste("PCA by", tools::toTitleCase(cc_group_fields[2])))

  if ("Phase" %in% colnames(combined_obj@meta.data)) {
    pca_plot3    <- DimPlot(combined_obj, reduction = "pca", group.by = "Phase") +
      ggtitle("PCA by Cell Cycle Phase")
    pca_combined <- (pca_plot1 | pca_plot2) / pca_plot3
  } else {
    pca_combined <- pca_plot1 | pca_plot2
  }
} else {
  pca_combined <- pca_plot1
}

ggsave(file.path(plots_dir, "pca_plots.pdf"), pca_combined, width = 12, height = 10)
cat("PCA plots saved\n")

elbow_plot <- ElbowPlot(combined_obj, ndims = n_pca_components)
ggsave(file.path(plots_dir, "pca_elbow_plot.pdf"), elbow_plot, width = 10, height = 6)
cat("Elbow plot saved\n")

# =============================================================================
# SAVE CHECKPOINT
# =============================================================================

saveRDS(combined_obj, ckpt$norm_pca)
cat("\nCheckpoint saved:", ckpt$norm_pca, "\n")
cat("\n=== STEP 4 COMPLETE ===\n")
