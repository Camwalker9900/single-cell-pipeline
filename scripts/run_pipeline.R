#!/usr/bin/env Rscript

# =============================================================================
# run_pipeline.R  —  Pipeline Orchestrator
#
# Usage:
#   Rscript run_pipeline.R                          # run all steps per config
#   Rscript run_pipeline.R --step qc                # run one step only
#   Rscript run_pipeline.R --from norm_pca          # run from a step onwards
#   Rscript run_pipeline.R --config path/to/cfg.yaml
#
# Steps: qc | cellbender | cell_cycle | norm_pca | integration | clustering
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# --- Parse --config ---
cfg_idx <- which(args == "--config")
if (length(cfg_idx) > 0) Sys.setenv(SCRNA_CONFIG = args[cfg_idx + 1])

# --- Parse --step and --from ---
step_idx <- which(args == "--step")
from_idx <- which(args == "--from")

run_single <- if (length(step_idx) > 0) args[step_idx + 1] else NULL
run_from   <- if (length(from_idx) > 0) args[from_idx + 1] else NULL

# Load config to read step toggles
source("scripts/00_load_config.R")

# -----------------------------------------------------------------------------
# Step registry — ordered
# -----------------------------------------------------------------------------
all_steps <- c("qc", "cellbender", "cell_cycle", "norm_pca", "integration", "clustering")

step_scripts <- list(
  qc          = "scripts/01_qc.R",
  cellbender  = "scripts/02_cellbender.R",
  cell_cycle  = "scripts/03_cell_cycle.R",
  norm_pca    = "scripts/04_normalization_pca.R",
  integration = "scripts/05_integration.R",
  clustering  = "scripts/06_clustering_umap.R"
)

# Which steps are enabled by config toggles
step_enabled <- list(
  qc          = TRUE,            # always required
  cellbender  = run_cellbender,
  cell_cycle  = run_cell_cycle,
  norm_pca    = TRUE,            # always required
  integration = run_integration,
  clustering  = TRUE             # always required
)

# -----------------------------------------------------------------------------
# Determine which steps to actually run
# -----------------------------------------------------------------------------
if (!is.null(run_single)) {
  if (!run_single %in% all_steps) {
    stop("Unknown step '", run_single, "'. Choose from: ", paste(all_steps, collapse = ", "))
  }
  steps_to_run <- run_single

} else if (!is.null(run_from)) {
  if (!run_from %in% all_steps) {
    stop("Unknown step '", run_from, "'. Choose from: ", paste(all_steps, collapse = ", "))
  }
  start_idx    <- which(all_steps == run_from)
  steps_to_run <- all_steps[start_idx:length(all_steps)]

} else {
  steps_to_run <- all_steps
}

# Apply config toggles (optional steps can still be force-run via --step)
if (is.null(run_single)) {
  steps_to_run <- steps_to_run[sapply(steps_to_run, function(s) isTRUE(step_enabled[[s]]))]
}

# -----------------------------------------------------------------------------
# Execute
# -----------------------------------------------------------------------------
cat(rep("=", 60), "\n")
cat("scRNA-seq Pipeline\n")
cat("Project:", project_name, "\n")
cat("Steps to run:", paste(steps_to_run, collapse = " → "), "\n")
cat(rep("=", 60), "\n\n")

for (step in steps_to_run) {
  script <- step_scripts[[step]]
  cat("\n", rep("-", 60), "\n")
  cat("▶ Running:", script, "\n")
  cat(rep("-", 60), "\n\n")

  tryCatch(
    source(script),
    error = function(e) {
      cat("\n❌ ERROR in step '", step, "':\n", e$message, "\n", sep = "")
      cat("Fix the error above and re-run from this step with:\n")
      cat("  Rscript run_pipeline.R --from", step, "\n\n")
      quit(status = 1)
    }
  )
}

cat("\n", rep("=", 60), "\n")
cat("✅ All steps completed successfully!\n")
cat("Final object: pipeline_outputs/final_object.rds\n")
cat(rep("=", 60), "\n")
