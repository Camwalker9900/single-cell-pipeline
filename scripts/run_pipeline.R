#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

cfg_idx <- which(args == "--config")
if (length(cfg_idx) > 0) Sys.setenv(SCRNA_CONFIG = args[cfg_idx + 1])

step_idx <- which(args == "--step")
from_idx <- which(args == "--from")

run_single <- if (length(step_idx) > 0) args[step_idx + 1] else NULL
run_from   <- if (length(from_idx) > 0) args[from_idx + 1] else NULL

source("scripts/00_load_config.R")

all_steps <- c(
  "qc", "cellbender", "cell_cycle", "norm_pca", "integration",
  "clustering", "annotation", "differential_expression", "downstream",
  "export"
)

step_scripts <- list(
  qc                    = "scripts/01_qc.R",
  cellbender            = "scripts/02_cellbender.R",
  cell_cycle            = "scripts/03_cell_cycle.R",
  norm_pca              = "scripts/04_normalization_pca.R",
  integration           = "scripts/05_integration.R",
  clustering            = "scripts/06_clustering_umap.R",
  annotation            = "scripts/07_annotation.R",
  differential_expression = "scripts/08_differential_expression.R",
  downstream            = "scripts/09_downstream_tools.R",
  export                = "scripts/10_export.R"
)

step_enabled <- list(
  qc                    = TRUE,
  cellbender            = run_cellbender,
  cell_cycle            = run_cell_cycle,
  norm_pca              = TRUE,
  integration           = run_integration,
  clustering            = TRUE,
  annotation            = run_annotation,
  differential_expression = run_de,
  downstream            = run_downstream,
  export                = run_export
)

if (!is.null(run_single)) {
  if (!run_single %in% all_steps) {
    stop("Unknown step '", run_single, "'. Choose from: ", paste(all_steps, collapse = ", "))
  }
  steps_to_run <- run_single
} else if (!is.null(run_from)) {
  if (!run_from %in% all_steps) {
    stop("Unknown step '", run_from, "'. Choose from: ", paste(all_steps, collapse = ", "))
  }
  start_idx <- which(all_steps == run_from)
  steps_to_run <- all_steps[start_idx:length(all_steps)]
} else {
  steps_to_run <- all_steps
}

if (is.null(run_single)) {
  steps_to_run <- steps_to_run[sapply(steps_to_run, function(s) isTRUE(step_enabled[[s]]))]
}

cat(rep("=", 60), "\n")
cat("scRNA-seq Pipeline\n")
cat("Project:", project_name, "\n")
cat("Steps to run:", paste(steps_to_run, collapse = " -> "), "\n")
cat(rep("=", 60), "\n\n")

for (step in steps_to_run) {
  script <- step_scripts[[step]]
  cat("\n", rep("-", 60), "\n")
  cat("Running:", script, "\n")
  cat(rep("-", 60), "\n\n")

  tryCatch(
    source(script),
    error = function(e) {
      cat("\nERROR in step '", step, "':\n", e$message, "\n", sep = "")
      cat("Fix the error above and re-run from this step with:\n")
      cat("  Rscript scripts/run_pipeline.R --from", step, "\n\n")
      quit(status = 1)
    }
  )
}

cat("\n", rep("=", 60), "\n")
cat("All requested steps completed successfully.\n")
cat("Latest object checkpoint:", if (file.exists(ckpt$downstream)) ckpt$downstream else if (file.exists(ckpt$annotation)) ckpt$annotation else ckpt$clustering, "\n")
cat(rep("=", 60), "\n")
