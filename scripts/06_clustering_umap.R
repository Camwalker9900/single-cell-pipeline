#!/usr/bin/env Rscript

# =============================================================================
# 06_clustering_umap.R  —  UMAP & Clustering
# Input : checkpoints/05_integration.rds  (if integration was run)
#         checkpoints/04_norm_pca.rds     (if integration was skipped)
# Output: checkpoints/06_clustering.rds  — final fully-processed object
#         pipeline_outputs/plots/final_umap_*.pdf
#         pipeline_outputs/plots/quality_control_metrics.pdf
#         pipeline_outputs/analysis_summary.csv
#         pipeline_outputs/cluster_summary.csv
#
# The reduction used for neighbors + UMAP is set by:
#   clustering.reduction in config.yaml  ("harmony" or "pca")
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

source("scripts/00_load_config.R")

cat("=== STEP 6: UMAP & Clustering ===\n\n")

# =============================================================================
# LOAD BEST AVAILABLE UPSTREAM CHECKPOINT
# Priority: integration > norm_pca
# =============================================================================

prev_ckpt <- if (file.exists(ckpt$integration)) ckpt$integration else
             if (file.exists(ckpt$norm_pca))     ckpt$norm_pca    else
             stop("No upstream checkpoint found. Run 04_normalization_pca.R first.")

cat("Loading checkpoint:", prev_ckpt, "\n")
combined_obj <- readRDS(prev_ckpt)
cat("Object:", ncol(combined_obj), "cells,", nrow(combined_obj), "genes\n")

# Validate reduction exists
if (!clustering_reduction %in% names(combined_obj@reductions)) {
  available <- paste(names(combined_obj@reductions), collapse = ", ")
  stop("Reduction '", clustering_reduction, "' not found in object.\n",
       "Available reductions: ", available, "\n",
       "Update 'clustering.reduction' in config.yaml.")
}
cat("Clustering reduction:", clustering_reduction, "\n\n")

# =============================================================================
# UMAP
# =============================================================================

cat("--- 6A: Computing UMAP ---\n")
combined_obj <- RunUMAP(combined_obj,
                        reduction = clustering_reduction,
                        dims      = 1:n_pca_components)

# =============================================================================
# NEIGHBORS & CLUSTERING
# =============================================================================

cat("--- 6B: Finding neighbors ---\n")
combined_obj <- FindNeighbors(combined_obj,
                              reduction = clustering_reduction,
                              dims      = 1:n_pca_components)

cat("--- 6C: Clustering at multiple resolutions ---\n")
for (res in clustering_resolutions) {
  combined_obj <- FindClusters(combined_obj, resolution = res)
  n_clusters   <- length(levels(combined_obj))
  cat("  Resolution", res, ":", n_clusters, "clusters\n")
}

Idents(combined_obj) <- paste0("RNA_snn_res.", selected_resolution)
cat("Active resolution set to:", selected_resolution, "\n")

# =============================================================================
# UMAP VISUALIZATIONS
# =============================================================================

cat("\n--- 6D: UMAP visualizations ---\n")

umap_meta_plots <- lapply(all_meta_fields, function(field) {
  DimPlot(combined_obj, reduction = "umap", group.by = field) +
    ggtitle(paste("UMAP by", tools::toTitleCase(field)))
})

umap_sample   <- DimPlot(combined_obj, reduction = "umap", group.by = "sample") +
  ggtitle("UMAP by Sample")
umap_clusters <- DimPlot(combined_obj, reduction = "umap", label = TRUE) +
  ggtitle(paste0("UMAP Clusters (res ", selected_resolution, ")"))

umap_combined <- wrap_plots(c(umap_meta_plots, list(umap_sample, umap_clusters)), ncol = 2)
ggsave(file.path(plots_dir, "final_umap_visualizations.pdf"), umap_combined, width = 16, height = 12)
cat("Final UMAP plots saved\n")

# Cell cycle UMAP (only if cell cycle scoring was done)
if ("Phase" %in% colnames(combined_obj@meta.data)) {
  cat("Creating cell cycle UMAP plots...\n")
  cc_group_fields <- head(all_meta_fields, 2)

  cc_final1 <- DimPlot(combined_obj, reduction = "umap", group.by = "Phase") +
    ggtitle("Cell Cycle Phase")
  cc_final2 <- DimPlot(combined_obj, reduction = "umap", group.by = "Phase",
                        split.by = cc_group_fields[1]) +
    ggtitle(paste("Cell Cycle by", tools::toTitleCase(cc_group_fields[1])))

  if (length(cc_group_fields) >= 2) {
    cc_final3 <- DimPlot(combined_obj, reduction = "umap", group.by = "Phase",
                          split.by = cc_group_fields[2]) +
      ggtitle(paste("Cell Cycle by", tools::toTitleCase(cc_group_fields[2])))
    cc_final_combined <- cc_final1 / (cc_final2 | cc_final3)
  } else {
    cc_final_combined <- cc_final1 / cc_final2
  }

  ggsave(file.path(plots_dir, "cell_cycle_final_umap.pdf"), cc_final_combined,
         width = 16, height = 12)
  cat("Cell cycle UMAP plots saved\n")
}

# QC metrics on UMAP
qc_features <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          colnames(combined_obj@meta.data))

qc_plots    <- lapply(qc_features, function(f) {
  FeaturePlot(combined_obj, features = f, reduction = "umap") + ggtitle(f)
})
qc_combined <- wrap_plots(qc_plots, nrow = 1)
ggsave(file.path(plots_dir, "quality_control_metrics.pdf"), qc_combined, width = 18, height = 6)
cat("QC metrics UMAP plot saved\n")

# =============================================================================
# SAVE SUMMARIES
# =============================================================================

cat("\n--- 6E: Saving summaries ---\n")

# Analysis summary
summary_data <- data.frame(
  Metric = c("Total Cells", "Total Genes", "Samples",
             paste0("Unique ", tools::toTitleCase(all_meta_fields)),
             "Variable Features", "PCA Components", "Clustering Reduction",
             "Clustering Resolution",
             if ("Phase" %in% colnames(combined_obj@meta.data)) "Cell Cycle Regression" else NULL,
             "Integration Method"),
  Value  = c(ncol(combined_obj), nrow(combined_obj), length(unique(combined_obj$sample)),
             sapply(all_meta_fields, function(f) length(unique(combined_obj@meta.data[[f]]))),
             n_variable_features, n_pca_components, clustering_reduction,
             selected_resolution,
             if ("Phase" %in% colnames(combined_obj@meta.data)) "S.Score + G2M.Score" else NULL,
             if (file.exists(ckpt$integration)) "Harmony" else "None"),
  stringsAsFactors = FALSE
)

write.csv(summary_data, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)
cat("Analysis summary saved\n")

# Cluster summary — dynamically built from actual metadata values
cluster_levels  <- levels(combined_obj)
cluster_summary <- data.frame(
  Cluster    = cluster_levels,
  Cell_Count = as.vector(table(Idents(combined_obj))),
  stringsAsFactors = FALSE
)

for (field in all_meta_fields) {
  field_values <- unique(combined_obj@meta.data[[field]])
  field_table  <- table(combined_obj@meta.data[[field]], combined_obj@active.ident)
  for (val in field_values) {
    col_name <- paste0(tools::toTitleCase(field), "_", val)
    cluster_summary[[col_name]] <- if (val %in% rownames(field_table)) {
      as.vector(field_table[val, ])
    } else {
      0L
    }
  }
}
cluster_summary[is.na(cluster_summary)] <- 0

write.csv(cluster_summary, file.path(output_dir, "cluster_summary.csv"), row.names = FALSE)
cat("Cluster summary saved\n")

# =============================================================================
# SAVE FINAL CHECKPOINT & OBJECT
# =============================================================================

saveRDS(combined_obj, ckpt$clustering)
cat("\nFinal checkpoint saved:", ckpt$clustering, "\n")

final_obj_file <- file.path(output_dir, "final_object.rds")
saveRDS(combined_obj, final_obj_file)
cat("Final object saved:", final_obj_file, "\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n", rep("=", 60), "\n")
cat("PIPELINE COMPLETE\n")
cat(rep("=", 60), "\n")
cat("Total cells:", ncol(combined_obj), "\n")
cat("Total genes:", nrow(combined_obj), "\n")
cat("Active clusters:", length(levels(combined_obj)),
    "(resolution", selected_resolution, ")\n")
cat("Reduction used:", clustering_reduction, "\n")
cat("\nAll plots saved to:", plots_dir, "\n")
cat("All summaries saved to:", output_dir, "\n")

cat("\nCluster Distribution (res", selected_resolution, "):\n")
print(table(Idents(combined_obj)))

for (field in all_meta_fields) {
  cat(paste0("\n", tools::toTitleCase(field), " by Cluster:\n"))
  print(table(combined_obj@meta.data[[field]], Idents(combined_obj)))
}

cat("\n🎉 Ready for cell type annotation!\n")
cat("\n=== STEP 6 COMPLETE ===\n")
