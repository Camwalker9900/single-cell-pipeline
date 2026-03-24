#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

source("scripts/00_load_config.R")
source("scripts/00_pipeline_utils.R")

cat("=== STEP 9: Downstream Tools ===\n\n")

prev_ckpt <- get_best_checkpoint_path(
  list(ckpt$annotation, ckpt$clustering),
  "No upstream checkpoint found. Run clustering or annotation first."
)
obj <- readRDS(prev_ckpt)
obj <- maybe_join_layers(obj, assay = "RNA")
DefaultAssay(obj) <- "RNA"

predicted_label_col <- if (annotation_prediction_col %in% colnames(obj@meta.data)) annotation_prediction_col else NULL

ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
get_embedding <- function(obj, reduction) {
  if (!reduction %in% names(obj@reductions)) stop("Reduction not present in object: ", reduction)
  Embeddings(obj, reduction)
}

# -----------------------------------------------------------------------------
# CellChat
# -----------------------------------------------------------------------------
if (isTRUE(cellchat_cfg$enabled)) {
  cat("--- 9A: CellChat ---\n")
  ensure_package("CellChat", "Install CellChat before enabling downstream.cellchat.enabled")
  library(CellChat)

  out_dir <- file.path(downstream_dir, "cellchat")
  ensure_dir(out_dir)

  label_col <- choose_first_present_col(obj, c(cellchat_cfg$label_col, predicted_label_col, choose_default_ident_col(obj, selected_resolution)))
  split_by <- cellchat_cfg$split_by
  db_choice <- tolower(cellchat_cfg$database %||% "human")
  min_cells <- cellchat_cfg$min_cells %||% 10

  if (is.null(label_col)) stop("Could not resolve a CellChat label column from config or object metadata")
  if (!is.null(split_by) && !split_by %in% colnames(obj@meta.data)) stop("CellChat split_by column not found: ", split_by)
  groups <- if (!is.null(split_by)) unique(as.character(obj@meta.data[[split_by]])) else "all"
  cellchat_list <- list()

  for (grp in groups) {
    sub_obj <- if (identical(grp, "all")) {
      obj
    } else {
      keep_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[split_by]]) == grp]
      subset(obj, cells = keep_cells)
    }
    labels <- as.character(sub_obj@meta.data[[label_col]])
    labels <- paste0("cluster_", safe_name(labels))
    names(labels) <- colnames(sub_obj)
    meta <- data.frame(labels = labels, row.names = colnames(sub_obj))
    data_input <- GetAssayData(sub_obj, assay = "RNA", layer = "data")

    cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
    cc@DB <- if (db_choice == "mouse") CellChatDB.mouse else CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- filterCommunication(cc, min.cells = min_cells)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cellchat_list[[as.character(grp)]] <- cc
    saveRDS(cc, file.path(out_dir, paste0("cellchat_", safe_name(grp), ".rds")))
    utils::write.csv(subsetCommunication(cc), file.path(out_dir, paste0("cellchat_", safe_name(grp), "_interactions.csv")), row.names = FALSE)
  }

  if (length(cellchat_list) >= 2) {
    merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    saveRDS(merged, file.path(out_dir, "cellchat_merged.rds"))
  }
}

# -----------------------------------------------------------------------------
# Slingshot
# -----------------------------------------------------------------------------
if (isTRUE(slingshot_cfg$enabled)) {
  cat("\n--- 9B: Slingshot ---\n")
  ensure_package("slingshot", "Install Bioconductor package 'slingshot'")
  ensure_package("SingleCellExperiment", "Install Bioconductor package 'SingleCellExperiment'")
  library(slingshot)
  library(SingleCellExperiment)

  out_dir <- file.path(downstream_dir, "slingshot")
  ensure_dir(out_dir)

  reduction <- slingshot_cfg$reduction %||% clustering_reduction
  cluster_col <- choose_first_present_col(obj, c(slingshot_cfg$cluster_col, predicted_label_col, choose_default_ident_col(obj, selected_resolution)))
  start_cluster <- slingshot_cfg$start_cluster
  end_clusters <- slingshot_cfg$end_clusters

  if (is.null(cluster_col)) stop("Could not resolve a Slingshot cluster column from config or object metadata")

  sce <- as.SingleCellExperiment(obj)
  reducedDims(sce)$REDUCED <- get_embedding(obj, reduction)
  colData(sce)$cluster_label <- as.character(obj@meta.data[[cluster_col]])

  sce <- slingshot(
    sce,
    clusterLabels = "cluster_label",
    reducedDim = "REDUCED",
    start.clus = start_cluster,
    end.clus = end_clusters
  )

  pt <- slingPseudotime(sce)
  if (!is.null(pt) && ncol(pt) > 0) {
    for (i in seq_len(ncol(pt))) {
      obj[[paste0("slingshot_pseudotime_", i)]] <- pt[, i]
    }
    pt_df <- data.frame(cell = rownames(pt), pt, stringsAsFactors = FALSE)
    write.csv(pt_df, file.path(out_dir, "pseudotime_by_cell.csv"), row.names = FALSE)
  }

  emb <- get_embedding(obj, reduction)
  png(file.path(out_dir, paste0(reduction, "_slingshot_overlay.png")), width = 1600, height = 1300, res = 200)
  plot(emb[, 1], emb[, 2], pch = 16, cex = 0.3, col = "grey60", xlab = paste0(reduction, "_1"), ylab = paste0(reduction, "_2"))
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  dev.off()
  saveRDS(sce, file.path(out_dir, "slingshot_sce.rds"))
}

# -----------------------------------------------------------------------------
# Monocle3
# -----------------------------------------------------------------------------
if (isTRUE(monocle_cfg$enabled)) {
  cat("\n--- 9C: Monocle3 ---\n")
  ensure_package("monocle3", "Install package 'monocle3' before enabling downstream.monocle.enabled")
  ensure_package("SeuratWrappers", "Install package 'SeuratWrappers' before enabling downstream.monocle.enabled")
  library(monocle3)
  library(SeuratWrappers)

  out_dir <- file.path(downstream_dir, "monocle")
  ensure_dir(out_dir)

  reduction <- monocle_cfg$reduction %||% "umap"
  root_col <- monocle_cfg$root_col
  root_values <- monocle_cfg$root_values
  label_col <- choose_first_present_col(obj, c(monocle_cfg$label_col, predicted_label_col, choose_default_ident_col(obj, selected_resolution)))

  cds <- as.cell_data_set(obj)
  for (col in colnames(obj@meta.data)) colData(cds)[[col]] <- obj@meta.data[[col]]
  if (reduction %in% names(obj@reductions)) cds@int_colData@listData$reducedDims$UMAP <- Embeddings(obj, reduction)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  if (!is.null(label_col) && label_col %in% colnames(colData(cds))) cds@clusters$UMAP$clusters <- factor(colData(cds)[[label_col]])
  cds <- learn_graph(cds)

  if (!is.null(root_col) && !is.null(root_values) && root_col %in% colnames(colData(cds))) {
    root_cells <- colnames(cds)[colData(cds)[[root_col]] %in% root_values]
    cds <- order_cells(cds, root_cells = root_cells)
  } else {
    cds <- order_cells(cds)
  }

  obj$monocle3_pseudotime <- pseudotime(cds)
  write.csv(data.frame(cell = colnames(cds), pseudotime = pseudotime(cds)), file.path(out_dir, "pseudotime_data.csv"), row.names = FALSE)
  ggsave(file.path(out_dir, "trajectory_pseudotime.pdf"), plot_cells(cds, color_cells_by = "pseudotime"), width = 10, height = 8)
  if (!is.null(label_col) && label_col %in% colnames(colData(cds))) {
    ggsave(file.path(out_dir, "trajectory_labels.pdf"), plot_cells(cds, color_cells_by = label_col, label_cell_groups = TRUE), width = 10, height = 8)
  }
  saveRDS(cds, file.path(out_dir, "monocle_cds.rds"))
}

# -----------------------------------------------------------------------------
# CytoTRACE2
# -----------------------------------------------------------------------------
if (isTRUE(cytotrace_cfg$enabled)) {
  cat("\n--- 9D: CytoTRACE2 ---\n")
  ensure_package("CytoTRACE2", "Install package 'CytoTRACE2' before enabling downstream.cytotrace.enabled")
  library(CytoTRACE2)

  out_dir <- file.path(downstream_dir, "cytotrace")
  ensure_dir(out_dir)

  expr <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "counts"))
  species <- cytotrace_cfg$species %||% "human"
  ncores <- cytotrace_cfg$n_cores %||% 1
  reduction <- cytotrace_cfg$plot_reduction %||% get_plot_reduction(obj)

  res <- CytoTRACE2::cytotrace2(expr, species = species, ncores = ncores)
  if ("CytoTRACE2_Score" %in% colnames(res)) obj$CytoTRACE2_Score <- res$CytoTRACE2_Score[colnames(obj)]
  if ("CytoTRACE2_Relative" %in% colnames(res)) obj$CytoTRACE2_Relative <- res$CytoTRACE2_Relative[colnames(obj)]
  if ("CytoTRACE2_Potency" %in% colnames(res)) obj$CytoTRACE2_Potency <- res$CytoTRACE2_Potency[colnames(obj)]
  write.csv(res, file.path(out_dir, "cytotrace_results.csv"), row.names = FALSE)

  if (!is.null(reduction) && "CytoTRACE2_Score" %in% colnames(obj@meta.data)) {
    p <- FeaturePlot(obj, features = "CytoTRACE2_Score", reduction = reduction) + ggtitle("CytoTRACE2 Score")
    ggsave(file.path(out_dir, "cytotrace_score.pdf"), p, width = 9, height = 7)
  }
}

# -----------------------------------------------------------------------------
# miloR
# -----------------------------------------------------------------------------
if (isTRUE(milo_cfg$enabled)) {
  cat("\n--- 9E: miloR ---\n")
  ensure_package("miloR", "Install Bioconductor package 'miloR' before enabling downstream.milo.enabled")
  ensure_package("SingleCellExperiment", "Install Bioconductor package 'SingleCellExperiment'")
  ensure_package("edgeR", "Install Bioconductor package 'edgeR'")
  library(miloR)
  library(SingleCellExperiment)

  out_dir <- file.path(downstream_dir, "milo")
  ensure_dir(out_dir)

  sample_col <- milo_cfg$sample_col %||% "sample"
  condition_col <- milo_cfg$condition_col %||% "condition"
  reduced_dim <- milo_cfg$reduction %||% "PCA"
  reduced_dim_slot <- if (identical(reduced_dim, "PCA") && "pca" %in% names(obj@reductions)) "PCA" else if (identical(reduced_dim, "pca") && "pca" %in% names(obj@reductions)) "pca" else reduced_dim
  d <- milo_cfg$d %||% min(30, n_pca_components)
  k <- milo_cfg$k %||% 30

  if (!sample_col %in% colnames(obj@meta.data)) stop("milo sample_col not found: ", sample_col)
  if (!condition_col %in% colnames(obj@meta.data)) stop("milo condition_col not found: ", condition_col)

  sce <- as.SingleCellExperiment(obj)
  if (identical(reduced_dim_slot, "PCA") && "pca" %in% names(obj@reductions)) {
    reducedDims(sce)$PCA <- Embeddings(obj, "pca")
  } else if (reduced_dim_slot %in% names(obj@reductions)) {
    reducedDims(sce)[[reduced_dim_slot]] <- Embeddings(obj, reduced_dim_slot)
  } else {
    stop("milo reduction not available in object: ", reduced_dim)
  }

  milo_obj <- Milo(sce)
  milo_obj <- buildGraph(milo_obj, k = k, d = d, reduced.dim = reduced_dim_slot)
  milo_obj <- makeNhoods(milo_obj, prop = milo_cfg$prop %||% 0.1, k = k, d = d, refined = TRUE, reduced_dims = reduced_dim_slot)

  meta_df <- data.frame(
    Sample = obj@meta.data[[sample_col]],
    Condition = obj@meta.data[[condition_col]],
    row.names = colnames(obj)
  )
  milo_obj <- countCells(milo_obj, meta.data = meta_df, sample = "Sample")

  design_df <- distinct(data.frame(Sample = obj@meta.data[[sample_col]], Condition = obj@meta.data[[condition_col]], stringsAsFactors = FALSE))
  rownames(design_df) <- design_df$Sample
  design <- stats::model.matrix(~ Condition, data = design_df)

  milo_res <- testNhoods(milo_obj, design = design, design.df = design_df)
  write.csv(milo_res, file.path(out_dir, "milo_results.csv"), row.names = FALSE)
  saveRDS(milo_obj, file.path(out_dir, "milo_object.rds"))
}

saveRDS(obj, ckpt$downstream)
saveRDS(obj, file.path(output_dir, "final_object_downstream.rds"))
cat("\nCheckpoint saved:", ckpt$downstream, "\n")
cat("=== STEP 9 COMPLETE ===\n")
