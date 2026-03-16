#!/usr/bin/env Rscript

# =============================================================================
# 03_cell_cycle.R  —  Cell Cycle Scoring & Regression  (OPTIONAL)
# Input : checkpoint from previous step (cellbender or qc)
#         — combined Seurat object
# Output: checkpoints/03_cell_cycle.rds  — combined object with S.Score,
#           G2M.Score, Phase added; data scaled with CC regression
#         pipeline_outputs/plots/cell_cycle_*.pdf
#
# NOTE: This step merges all samples into one combined object and scales the
#       data. If you skip this step, 04_normalization_pca.R will scale without
#       cell cycle regression.
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)

source("scripts/00_load_config.R")

cat("=== STEP 3: Cell Cycle Scoring & Regression ===\n\n")

# =============================================================================
# LOAD PREVIOUS CHECKPOINT
# Prefer cellbender checkpoint if it exists, otherwise fall back to qc
# =============================================================================

prev_ckpt <- if (file.exists(ckpt$cellbender)) ckpt$cellbender else ckpt$qc

if (!file.exists(prev_ckpt)) stop("No upstream checkpoint found. Run 01_qc.R first.")
cat("Loading checkpoint:", prev_ckpt, "\n")

input_data <- readRDS(prev_ckpt)

# =============================================================================
# MERGE PER-SAMPLE LIST INTO COMBINED OBJECT (if not already merged)
# =============================================================================

if (is.list(input_data) && !inherits(input_data, "Seurat")) {
  cat("Merging", length(input_data), "samples into combined object...\n")

  sample_list <- input_data

  # Add unique cell prefixes to avoid barcode conflicts
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
  cat("Object already merged:", ncol(combined_obj), "cells\n\n")
}

# =============================================================================
# NORMALIZE (required before CellCycleScoring)
# =============================================================================

cat("--- 3A: Normalizing data ---\n")
combined_obj <- NormalizeData(combined_obj,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)

cat("--- 3B: Finding variable features ---\n")
combined_obj <- FindVariableFeatures(combined_obj,
                                     selection.method = "vst",
                                     nfeatures = n_variable_features)

# =============================================================================
# CELL CYCLE SCORING
# =============================================================================

cat("\n--- 3C: Cell Cycle Scoring ---\n")

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cat("S phase genes:", length(s.genes), "| G2/M phase genes:", length(g2m.genes), "\n")

combined_obj <- CellCycleScoring(combined_obj,
                                 s.features   = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident    = TRUE)

cat("Cell cycle phase distribution:\n")
print(table(combined_obj@meta.data$Phase))

# =============================================================================
# CELL CYCLE PLOTS (bar plots — before regression)
# =============================================================================

cat("\n--- 3D: Cell cycle distribution plots ---\n")

cc_group_fields <- head(all_meta_fields, 2)

cc_summary <- combined_obj@meta.data %>%
  group_by(across(all_of(c(cc_group_fields, "Phase")))) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(across(all_of(cc_group_fields))) %>%
  mutate(percentage = count / sum(count) * 100)

cc_plot1 <- ggplot(cc_summary, aes(x = .data[[cc_group_fields[1]]], y = percentage, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = paste("Cell Cycle Phase Distribution by", tools::toTitleCase(cc_group_fields[1])),
       y = "Percentage of Cells", x = tools::toTitleCase(cc_group_fields[1])) +
  theme_minimal() + scale_fill_brewer(palette = "Set2")

if (length(cc_group_fields) >= 2) {
  cc_plot2 <- ggplot(cc_summary, aes(x = .data[[cc_group_fields[2]]], y = percentage, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Cell Cycle Phase Distribution by", tools::toTitleCase(cc_group_fields[2])),
         y = "Percentage of Cells", x = tools::toTitleCase(cc_group_fields[2])) +
    theme_minimal() + scale_fill_brewer(palette = "Set2")

  cc_plot3 <- ggplot(cc_summary,
    aes(x = interaction(.data[[cc_group_fields[1]]], .data[[cc_group_fields[2]]]),
        y = percentage, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Cell Cycle by", tools::toTitleCase(cc_group_fields[1]),
                       "and", tools::toTitleCase(cc_group_fields[2])),
         y = "Percentage of Cells",
         x = paste0(tools::toTitleCase(cc_group_fields[1]), ".",
                    tools::toTitleCase(cc_group_fields[2]))) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")

  cc_combined_plot <- (cc_plot1 | cc_plot2) / cc_plot3
} else {
  cc_combined_plot <- cc_plot1
}

ggsave(file.path(plots_dir, "cell_cycle_phases.pdf"), cc_combined_plot, width = 12, height = 10)
cat("Cell cycle distribution plots saved\n")

# =============================================================================
# SCALE WITH CELL CYCLE REGRESSION
# =============================================================================

cat("\n--- 3E: Scaling data with cell cycle regression ---\n")
cat("This may take several minutes for large datasets...\n")

all_genes    <- rownames(combined_obj)
combined_obj <- ScaleData(combined_obj,
                          features       = all_genes,
                          vars.to.regress = c("S.Score", "G2M.Score"))

cat("Data scaling complete\n")

# =============================================================================
# SAVE CHECKPOINT
# =============================================================================

saveRDS(combined_obj, ckpt$cell_cycle)
cat("\nCheckpoint saved:", ckpt$cell_cycle, "\n")
cat("\n=== STEP 3 COMPLETE ===\n")
