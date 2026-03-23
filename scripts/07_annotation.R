#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/00_load_config.R")
source("scripts/00_pipeline_utils.R")

cat("=== STEP 7: Annotation ===\n\n")

prev_ckpt <- get_best_checkpoint_path(
  list(ckpt$clustering),
  "Clustering checkpoint not found. Run 06_clustering_umap.R first."
)
obj <- readRDS(prev_ckpt)
obj <- maybe_join_layers(obj, assay = "RNA")

ident_col <- annotation_ident_col %||% choose_default_ident_col(obj, selected_resolution = selected_resolution)
if (!is.null(ident_col)) obj <- set_idents_from_column(obj, ident_col)
plot_reduction <- get_plot_reduction(obj)

cat("Using identity column:", ident_col %||% "active.ident", "\n")
cat("Annotation methods enabled:", paste(names(annotation_methods)[unlist(annotation_methods)], collapse = ", "), "\n\n")

if (isTRUE(annotation_methods$findallmarkers)) {
  cat("--- 7A: FindAllMarkers ---\n")
  fam_only_pos <- annotation_findall_cfg$only_pos
  if (is.null(fam_only_pos)) fam_only_pos <- TRUE
  fam_min_pct <- annotation_findall_cfg$min_pct %||% 0.25
  fam_logfc <- annotation_findall_cfg$logfc_threshold %||% 0.25
  fam_test <- annotation_findall_cfg$test_use %||% "wilcox"

  markers <- FindAllMarkers(
    obj,
    only.pos = fam_only_pos,
    min.pct = fam_min_pct,
    logfc.threshold = fam_logfc,
    test.use = fam_test
  )

  write.csv(markers, file.path(annotation_dir, "findallmarkers_all.csv"), row.names = FALSE)
  save_top_markers(markers, file.path(annotation_dir, "findallmarkers_top10.csv"), n = 10)
  cat("Saved marker tables to:", annotation_dir, "\n")
}

if (isTRUE(annotation_methods$module_scores)) {
  cat("\n--- 7B: Module scores ---\n")
  gene_sets <- read_gene_sets(annotation_gene_set_file, annotation_gene_sets)
  if (length(gene_sets) == 0) {
    stop("Module score method enabled but no gene sets were provided via annotation.module_scores.gene_sets or gene_set_file")
  }

  used_sets <- list()
  for (set_name in names(gene_sets)) {
    genes_present <- intersect(unique(gene_sets[[set_name]]), rownames(obj))
    if (length(genes_present) < 2) next

    score_prefix <- paste0("MS_", safe_name(set_name))
    obj <- AddModuleScore(
      obj,
      features = list(genes_present),
      name = paste0(score_prefix, "_"),
      assay = "RNA",
      seed = 1,
      search = FALSE
    )

    raw_col <- grep(paste0("^", score_prefix, "_"), colnames(obj@meta.data), value = TRUE)[1]
    obj[[score_prefix]] <- obj[[raw_col]]
    used_sets[[set_name]] <- genes_present

    if (!is.null(plot_reduction)) {
      p <- FeaturePlot(obj, features = score_prefix, reduction = plot_reduction) +
        ggtitle(paste("Module score:", set_name))
      ggsave(file.path(annotation_dir, paste0(score_prefix, "_", plot_reduction, ".pdf")), p, width = 8, height = 6)
    }
  }

  if (length(used_sets) == 0) {
    stop("None of the configured module score gene sets had at least two genes present in the object")
  }

  used_df <- do.call(rbind, lapply(names(used_sets), function(nm) data.frame(set = nm, gene = used_sets[[nm]], stringsAsFactors = FALSE)))
  write.csv(used_df, file.path(annotation_dir, "module_score_gene_sets_used.csv"), row.names = FALSE)

  ms_cols <- grep("^MS_", colnames(obj@meta.data), value = TRUE)
  ms_summary <- data.frame(
    score = ms_cols,
    mean = sapply(ms_cols, function(col) mean(obj@meta.data[[col]], na.rm = TRUE)),
    median = sapply(ms_cols, function(col) median(obj@meta.data[[col]], na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  write.csv(ms_summary, file.path(annotation_dir, "module_score_summary.csv"), row.names = FALSE)
  cat("Saved module score outputs to:", annotation_dir, "\n")
}

if (isTRUE(annotation_methods$label_transfer)) {
  cat("\n--- 7C: Label transfer ---\n")
  if (is.null(annotation_reference_rds) || !file.exists(annotation_reference_rds)) {
    stop("Label transfer enabled but annotation.label_transfer.reference_rds is missing or does not exist")
  }
  if (is.null(annotation_reference_label)) {
    stop("Label transfer enabled but annotation.label_transfer.reference_label_col is not set")
  }

  ref_obj <- readRDS(annotation_reference_rds)
  ref_obj <- maybe_join_layers(ref_obj, assay = annotation_reference_assay)
  if (!annotation_reference_label %in% colnames(ref_obj@meta.data)) {
    stop("Reference label column not found in reference object: ", annotation_reference_label)
  }

  common_genes <- intersect(rownames(obj), rownames(ref_obj))
  if (length(common_genes) < 200) {
    stop("Too few shared genes for label transfer: ", length(common_genes))
  }

  DefaultAssay(ref_obj) <- annotation_reference_assay
  DefaultAssay(obj) <- "RNA"

  if (!"pca" %in% names(ref_obj@reductions)) {
    ref_obj <- NormalizeData(ref_obj)
    ref_obj <- FindVariableFeatures(ref_obj, nfeatures = min(length(common_genes), n_variable_features))
    ref_obj <- ScaleData(ref_obj, verbose = FALSE)
    ref_obj <- RunPCA(ref_obj, npcs = max(annotation_transfer_dims), verbose = FALSE)
  }

  if (!"pca" %in% names(obj@reductions)) {
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, nfeatures = min(length(common_genes), n_variable_features))
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(annotation_transfer_dims), verbose = FALSE)
  }

  anchors <- FindTransferAnchors(
    reference = ref_obj,
    query = obj,
    reference.assay = annotation_reference_assay,
    query.assay = "RNA",
    features = common_genes,
    dims = annotation_transfer_dims,
    reference.reduction = "pca",
    normalization.method = "LogNormalize"
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref_obj[[annotation_reference_label, drop = TRUE]],
    dims = annotation_transfer_dims,
    weight.reduction = obj[["pca"]]
  )

  pred_col <- if ("predicted.id" %in% colnames(predictions)) "predicted.id" else colnames(predictions)[1]
  score_candidates <- intersect(c("prediction.score.max", annotation_score_col), colnames(predictions))
  score_col <- if (length(score_candidates) > 0) score_candidates[1] else NULL

  obj[[annotation_prediction_col]] <- predictions[[pred_col]]
  if (!is.null(score_col)) obj[[annotation_score_col]] <- predictions[[score_col]]

  write.csv(predictions, file.path(annotation_dir, "label_transfer_predictions.csv"), row.names = TRUE)
  saveRDS(anchors, file.path(annotation_dir, "label_transfer_anchors.rds"))

  if (!is.null(plot_reduction)) {
    p <- DimPlot(obj, reduction = plot_reduction, group.by = annotation_prediction_col, label = TRUE) +
      ggtitle("Transferred labels")
    ggsave(file.path(annotation_dir, paste0("label_transfer_", plot_reduction, ".pdf")), p, width = 10, height = 8)
  }
  cat("Saved label transfer outputs to:", annotation_dir, "\n")
}

if (annotation_tangram_enabled) {
  cat("\n--- 7D: Tangram ---\n")
  if (is.null(annotation_tangram_query_path) || identical(annotation_tangram_query_path, "")) {
    stop("Tangram enabled but annotation.tangram.query_path is not set")
  }
  if (!file.exists(annotation_tangram_script)) {
    stop("Tangram helper script not found: ", annotation_tangram_script)
  }

  tangram_ref_h5ad <- annotation_tangram_ref_h5ad
  if (is.null(tangram_ref_h5ad) || identical(tangram_ref_h5ad, "")) {
    ensure_package("SingleCellExperiment", "Install Bioconductor package 'SingleCellExperiment'")
    ensure_package("zellkonverter", "Install Bioconductor package 'zellkonverter' to export H5AD for Tangram")
    library(SingleCellExperiment)
    library(zellkonverter)

    export_obj <- obj
    if (!annotation_tangram_label_col %in% colnames(export_obj@meta.data)) {
      if (!is.null(ident_col)) {
        export_obj[[annotation_tangram_label_col]] <- export_obj[[ident_col, drop = TRUE]]
      } else {
        stop("Tangram label column not found in object metadata: ", annotation_tangram_label_col)
      }
    }

    sce <- as.SingleCellExperiment(export_obj)
    tangram_ref_h5ad <- file.path(exports_dir, "pipeline_reference_for_tangram.h5ad")
    zellkonverter::writeH5AD(sce, tangram_ref_h5ad)
  }

  tangram_out_prefix <- file.path(annotation_dir, "tangram")
  args <- c(
    annotation_tangram_script,
    "--ad_sc", tangram_ref_h5ad,
    "--out_prefix", tangram_out_prefix,
    "--celltype_col", annotation_tangram_label_col,
    "--method", "tangram"
  )
  if (identical(tolower(annotation_tangram_query_type), "xenium")) {
    args <- c(args, "--xenium_dir", annotation_tangram_query_path)
  } else {
    args <- c(args, "--query_h5ad", annotation_tangram_query_path)
  }

  status <- system2(annotation_tangram_python, args = args)
  if (!identical(status, 0L)) stop("Tangram helper script failed with exit code ", status)
  cat("Tangram outputs saved with prefix:", tangram_out_prefix, "\n")
}

saveRDS(obj, ckpt$annotation)
cat("\nCheckpoint saved:", ckpt$annotation, "\n")
cat("\n=== STEP 7 COMPLETE ===\n")
