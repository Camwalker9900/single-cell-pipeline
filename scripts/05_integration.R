#!/usr/bin/env Rscript

# =============================================================================
# 05_integration.R  —  Harmony Batch Integration  (OPTIONAL)
# Input : checkpoints/04_norm_pca.rds
# Output: checkpoints/05_integration.rds  — object with harmony reduction added
#         pipeline_outputs/plots/harmony_integration.pdf
#
# If steps.run_integration is false in config this script is skipped.
# 06_clustering_umap.R will use clustering.reduction from config instead.
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

source("scripts/00_load_config.R")

cat("=== STEP 5: Harmony Integration ===\n\n")

# =============================================================================
# LOAD PREVIOUS CHECKPOINT
# =============================================================================

if (!file.exists(ckpt$norm_pca)) stop("norm_pca checkpoint not found: ", ckpt$norm_pca,
                                       "\nRun 04_normalization_pca.R first.")
cat("Loading checkpoint:", ckpt$norm_pca, "\n")
combined_obj <- readRDS(ckpt$norm_pca)
cat("Object:", ncol(combined_obj), "cells,", nrow(combined_obj), "genes\n\n")

# =============================================================================
# HARMONY INTEGRATION
# =============================================================================

cat("--- 5A: Running Harmony integration on '", harmony_group_vars, "' ---\n", sep = "")

combined_obj <- RunHarmony(combined_obj,
                           group.by.vars = harmony_group_vars,
                           reduction     = "pca",
                           dims          = 1:n_pca_components,
                           assay         = "RNA")

cat("Harmony integration complete\n")

# =============================================================================
# POST-HARMONY PLOTS
# =============================================================================

cat("\n--- 5B: Post-Harmony visualizations ---\n")

harmony_plots <- lapply(all_meta_fields, function(field) {
  DimPlot(combined_obj, reduction = "harmony", group.by = field) +
    ggtitle(paste("Harmony by", tools::toTitleCase(field)))
})

harmony_plots[["sample"]] <- DimPlot(combined_obj, reduction = "harmony", group.by = "sample") +
  ggtitle("Harmony by Sample")

harmony_combined <- wrap_plots(harmony_plots, ncol = 2)
ggsave(file.path(plots_dir, "harmony_integration.pdf"), harmony_combined, width = 12, height = 10)
cat("Harmony plots saved\n")

# =============================================================================
# SAVE CHECKPOINT
# =============================================================================

saveRDS(combined_obj, ckpt$integration)
cat("\nCheckpoint saved:", ckpt$integration, "\n")
cat("\n=== STEP 5 COMPLETE ===\n")
