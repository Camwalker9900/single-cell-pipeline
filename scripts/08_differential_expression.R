#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("scripts/00_load_config.R")
source("scripts/00_pipeline_utils.R")

cat("=== STEP 8: Differential Expression ===\n\n")

prev_ckpt <- get_best_checkpoint_path(
  list(ckpt$annotation, ckpt$clustering),
  "No upstream checkpoint found. Run clustering or annotation first."
)
obj <- readRDS(prev_ckpt)
obj <- maybe_join_layers(obj, assay = "RNA")
DefaultAssay(obj) <- "RNA"

if (length(de_comparisons) == 0) {
  stop("Differential expression step enabled but no differential_expression.comparisons were provided in config.yaml")
}

results_summary <- list()

for (cmp in de_comparisons) {
  cmp_name <- cmp$name %||% paste(cmp$ident_1, "vs", cmp$ident_2, sep = "_")
  cat("--- Running comparison:", cmp_name, "---\n")

  cmp_obj <- obj
  if (!is.null(cmp$subset_col) && !is.null(cmp$subset_values)) {
    if (!cmp$subset_col %in% colnames(cmp_obj@meta.data)) {
      stop("Subset column not found for comparison ", cmp_name, ": ", cmp$subset_col)
    }
    keep <- cmp_obj@meta.data[[cmp$subset_col]] %in% cmp$subset_values
    cmp_obj <- subset(cmp_obj, cells = rownames(cmp_obj@meta.data)[keep])
  }

  group_col <- cmp$group_by %||% choose_default_ident_col(cmp_obj, selected_resolution = selected_resolution)
  if (is.null(group_col)) stop("Could not determine group_by column for comparison ", cmp_name)
  cmp_obj <- set_idents_from_column(cmp_obj, group_col)

  test_use <- cmp$test_use %||% "wilcox"
  min_pct <- cmp$min_pct %||% 0.1
  logfc_threshold <- cmp$logfc_threshold %||% 0.25
  latent_vars <- cmp$latent_vars

  markers <- FindMarkers(
    cmp_obj,
    ident.1 = cmp$ident_1,
    ident.2 = cmp$ident_2,
    test.use = test_use,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    latent.vars = latent_vars
  )

  if (nrow(markers) > 0) {
    markers$gene <- rownames(markers)
    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[1]
    p_col <- intersect(c("p_val_adj", "p_val"), colnames(markers))[1]
    sig <- if (!is.na(fc_col) && !is.na(p_col)) {
      markers[abs(markers[[fc_col]]) >= logfc_threshold & markers[[p_col]] < 0.05, , drop = FALSE]
    } else {
      markers
    }
  } else {
    sig <- markers
  }

  cmp_dir <- file.path(de_dir, safe_name(cmp_name))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(markers, file.path(cmp_dir, "all_results.csv"), row.names = FALSE)
  write.csv(sig, file.path(cmp_dir, "significant_results.csv"), row.names = FALSE)

  results_summary[[cmp_name]] <- data.frame(
    comparison = cmp_name,
    group_by = group_col,
    ident_1 = cmp$ident_1,
    ident_2 = cmp$ident_2,
    total_genes = nrow(markers),
    significant_genes = nrow(sig),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, results_summary)
write.csv(summary_df, file.path(de_dir, "comparison_summary.csv"), row.names = FALSE)

saveRDS(obj, ckpt$de)
cat("\nCheckpoint saved:", ckpt$de, "\n")
cat("=== STEP 8 COMPLETE ===\n")
