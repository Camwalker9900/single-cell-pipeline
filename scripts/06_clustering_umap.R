#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

source("scripts/00_load_config.R")
source("scripts/00_pipeline_utils.R")

cat("=== STEP 6: Clustering & UMAP ===\n\n")

prev_ckpt <- get_best_checkpoint_path(
  list(ckpt$integration, ckpt$norm_pca),
  "No upstream checkpoint found. Run 04_normalization_pca.R first."
)
cat("Loading checkpoint:", prev_ckpt, "\n")
combined_obj <- readRDS(prev_ckpt)
cat("Object:", ncol(combined_obj), "cells,", nrow(combined_obj), "genes\n")

if (!clustering_reduction %in% names(combined_obj@reductions)) {
  available <- paste(names(combined_obj@reductions), collapse = ", ")
  stop("Reduction '", clustering_reduction, "' not found in object. Available reductions: ", available)
}

combined_obj <- maybe_join_layers(combined_obj, assay = "RNA")
cat("Clustering method:", clustering_method, "| reduction:", clustering_reduction, "\n\n")

combined_obj <- RunUMAP(combined_obj, reduction = clustering_reduction, dims = 1:clustering_dims)

if (identical(tolower(clustering_method), "choir")) {
  ensure_package("CHOIR", "Install the CHOIR package before using clustering.method: choir")
  library(CHOIR)

  cat("--- 6A: Running CHOIR clustering ---\n")
  combined_obj <- CHOIR(
    object = combined_obj,
    alpha = choir_alpha,
    n_cores = choir_n_cores,
    verbose = TRUE
  )

  if (!choir_cluster_column %in% colnames(combined_obj@meta.data)) {
    stop("Expected CHOIR cluster column not found: ", choir_cluster_column)
  }

  combined_obj$CHOIR_cluster <- as.character(combined_obj[[choir_cluster_column, drop = TRUE]])
  Idents(combined_obj) <- factor(combined_obj$CHOIR_cluster)
  active_ident_col <- "CHOIR_cluster"
  cluster_levels <- levels(Idents(combined_obj))
  cat("CHOIR identified", length(cluster_levels), "clusters\n")
} else {
  cat("--- 6A: Finding neighbors ---\n")
  combined_obj <- FindNeighbors(combined_obj, reduction = clustering_reduction, dims = 1:clustering_dims)

  cat("--- 6B: Clustering at multiple resolutions ---\n")
  for (res in clustering_resolutions) {
    combined_obj <- FindClusters(combined_obj, resolution = res)
    cat("  Resolution", res, ":", length(levels(combined_obj)), "clusters\n")
  }

  active_ident_col <- paste0("RNA_snn_res.", selected_resolution)
  Idents(combined_obj) <- active_ident_col
  cluster_levels <- levels(Idents(combined_obj))
  cat("Active resolution set to:", selected_resolution, "\n")
}

cat("\n--- 6C: UMAP visualizations ---\n")
umap_meta_plots <- lapply(all_meta_fields, function(field) {
  DimPlot(combined_obj, reduction = "umap", group.by = field) +
    ggtitle(paste("UMAP by", tools::toTitleCase(field)))
})

umap_sample <- DimPlot(combined_obj, reduction = "umap", group.by = "sample") +
  ggtitle("UMAP by Sample")
umap_clusters <- DimPlot(combined_obj, reduction = "umap", label = TRUE) +
  ggtitle(paste("UMAP Clusters (", active_ident_col, ")", sep = ""))

umap_combined <- wrap_plots(c(umap_meta_plots, list(umap_sample, umap_clusters)), ncol = 2)
ggsave(file.path(plots_dir, "final_umap_visualizations.pdf"), umap_combined, width = 16, height = 12)

if ("Phase" %in% colnames(combined_obj@meta.data)) {
  cc_group_fields <- head(all_meta_fields, 2)
  cc_final1 <- DimPlot(combined_obj, reduction = "umap", group.by = "Phase") + ggtitle("Cell Cycle Phase")
  cc_final2 <- DimPlot(combined_obj, reduction = "umap", group.by = "Phase", split.by = cc_group_fields[1]) +
    ggtitle(paste("Cell Cycle by", tools::toTitleCase(cc_group_fields[1])))

  cc_final_combined <- if (length(cc_group_fields) >= 2) {
    cc_final3 <- DimPlot(combined_obj, reduction = "umap", group.by = "Phase", split.by = cc_group_fields[2]) +
      ggtitle(paste("Cell Cycle by", tools::toTitleCase(cc_group_fields[2])))
    cc_final1 / (cc_final2 | cc_final3)
  } else {
    cc_final1 / cc_final2
  }

  ggsave(file.path(plots_dir, "cell_cycle_final_umap.pdf"), cc_final_combined, width = 16, height = 12)
}

qc_features <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"), colnames(combined_obj@meta.data))
qc_plots <- lapply(qc_features, function(f) FeaturePlot(combined_obj, features = f, reduction = "umap") + ggtitle(f))
if (length(qc_plots) > 0) {
  qc_combined <- wrap_plots(qc_plots, nrow = 1)
  ggsave(file.path(plots_dir, "quality_control_metrics.pdf"), qc_combined, width = 18, height = 6)
}

cat("\n--- 6D: Saving summaries ---\n")
summary_data <- data.frame(
  Metric = c(
    "Total Cells", "Total Genes", "Samples",
    paste0("Unique ", tools::toTitleCase(all_meta_fields)),
    "Variable Features", "PCA Components", "Clustering Reduction",
    "Clustering Method", "Active Identity Column",
    if (identical(tolower(clustering_method), "seurat")) "Clustering Resolution" else "CHOIR Alpha",
    if ("Phase" %in% colnames(combined_obj@meta.data)) "Cell Cycle Regression" else NULL,
    "Integration Method"
  ),
  Value = c(
    ncol(combined_obj), nrow(combined_obj), length(unique(combined_obj$sample)),
    sapply(all_meta_fields, function(f) length(unique(combined_obj@meta.data[[f]]))),
    n_variable_features, n_pca_components, clustering_reduction,
    clustering_method, active_ident_col,
    if (identical(tolower(clustering_method), "seurat")) selected_resolution else choir_alpha,
    if ("Phase" %in% colnames(combined_obj@meta.data)) "S.Score + G2M.Score" else NULL,
    if (file.exists(ckpt$integration)) "Harmony" else "None"
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_data, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

cluster_summary <- data.frame(
  Cluster = cluster_levels,
  Cell_Count = as.vector(table(Idents(combined_obj))),
  stringsAsFactors = FALSE
)
for (field in all_meta_fields) {
  field_values <- unique(combined_obj@meta.data[[field]])
  field_table <- table(combined_obj@meta.data[[field]], Idents(combined_obj))
  for (val in field_values) {
    col_name <- paste0(tools::toTitleCase(field), "_", safe_name(val))
    cluster_summary[[col_name]] <- if (val %in% rownames(field_table)) as.vector(field_table[val, ]) else 0L
  }
}
cluster_summary[is.na(cluster_summary)] <- 0
write.csv(cluster_summary, file.path(output_dir, "cluster_summary.csv"), row.names = FALSE)

saveRDS(combined_obj, ckpt$clustering)
saveRDS(combined_obj, file.path(output_dir, "final_object.rds"))

cat("\nPIPELINE CLUSTERING COMPLETE\n")
cat("Clusters:", length(cluster_levels), "\n")
cat("Checkpoint saved:", ckpt$clustering, "\n")
cat("\n=== STEP 6 COMPLETE ===\n")
