# AGENTS.md

## Purpose

This repository is a checkpoint-based scRNA-seq pipeline in R centered on Seurat objects. It now covers:

- QC and optional CellBender filtering
- optional cell-cycle regression
- normalization, PCA, and optional Harmony integration
- clustering with either Seurat or CHOIR
- optional annotation methods
- optional differential expression
- optional downstream analyses including CellChat, Slingshot, Monocle3, CytoTRACE2, and miloR

## Entry Points

Use the current repo layout, not the stale examples in older docs.

- Full pipeline: `Rscript scripts/run_pipeline.R`
- Single step: `Rscript scripts/run_pipeline.R --step clustering`
- Resume: `Rscript scripts/run_pipeline.R --from annotation`
- Alternate config: `Rscript scripts/run_pipeline.R --config path/to/config.yaml`
- Install baseline dependencies: `Rscript scripts/requirements.R`

## Step Order

1. `scripts/01_qc.R`
2. `scripts/02_cellbender.R` optional
3. `scripts/03_cell_cycle.R` optional
4. `scripts/04_normalization_pca.R`
5. `scripts/05_integration.R` optional
6. `scripts/06_clustering_umap.R`
7. `scripts/07_annotation.R` optional
8. `scripts/08_differential_expression.R` optional
9. `scripts/09_downstream_tools.R` optional

The orchestrator sources each step script directly. Shared config and paths come from [`scripts/00_load_config.R`](/home/cam/github/single-cell-pipeline/scripts/00_load_config.R), and shared helper logic lives in [`scripts/00_pipeline_utils.R`](/home/cam/github/single-cell-pipeline/scripts/00_pipeline_utils.R).

## Data Contracts

- Steps 1 and 2 work on a named list of per-sample Seurat objects.
- Steps 3 onward work on a single merged Seurat object.
- Checkpoints are the real control plane. Downstream scripts choose the best available upstream checkpoint by file existence.
- New checkpoints now continue through:
  - `checkpoints/07_annotation.rds`
  - `checkpoints/08_differential_expression.rds`
  - `checkpoints/09_downstream_tools.rds`

## Clustering

[`scripts/06_clustering_umap.R`](/home/cam/github/single-cell-pipeline/scripts/06_clustering_umap.R) supports two modes:

- `clustering.method: seurat`
  Uses `FindNeighbors` plus multi-resolution `FindClusters`, then sets identities to `RNA_snn_res.<selected_resolution>`.
- `clustering.method: choir`
  Runs CHOIR on the merged object, copies the chosen CHOIR metadata column into `CHOIR_cluster`, and uses that as the active identity.

UMAP still uses `clustering.reduction` regardless of clustering backend.

## Annotation

[`scripts/07_annotation.R`](/home/cam/github/single-cell-pipeline/scripts/07_annotation.R) can run any enabled combination of:

- `FindAllMarkers`
- `AddModuleScore` using inline gene sets or a CSV/TSV/YAML gene-set file
- Seurat label transfer from a configured reference RDS
- Tangram through the Python helper [`scripts/tangram_annotate.py`](/home/cam/github/single-cell-pipeline/scripts/tangram_annotate.py)

Important:

- Tangram is only meaningful when a spatial query is configured. This pipeline remains scRNA-first; Tangram acts as an optional external mapping hook, not a direct annotation method on the scRNA object itself.
- Label transfer expects a reference Seurat object and a valid metadata label column.

## Differential Expression

[`scripts/08_differential_expression.R`](/home/cam/github/single-cell-pipeline/scripts/08_differential_expression.R) runs only the comparisons declared in `config.yaml` under `differential_expression.comparisons`.

Each comparison can specify:

- a name
- a `group_by` column
- `ident_1` and `ident_2`
- an optional metadata subset before testing
- `test_use`, `min_pct`, and `logfc_threshold`

Outputs land in `pipeline_outputs/differential_expression/<comparison_name>/`.

## Downstream Tools

[`scripts/09_downstream_tools.R`](/home/cam/github/single-cell-pipeline/scripts/09_downstream_tools.R) is modular. Each tool has its own config block and runs only when enabled.

Implemented hooks:

- CellChat
- Slingshot
- Monocle3
- CytoTRACE2
- miloR

Notes:

- These are optional by design because package availability and project assumptions vary a lot.
- Several tools require additional metadata discipline, especially `condition`, `sample`, root-state definitions, or label columns.

## Config Guidance

Treat [`config.yaml`](/home/cam/github/single-cell-pipeline/config.yaml) as the authoritative interface.

High-value sections for new work:

- `steps.*`
- `clustering.*`
- `annotation.*`
- `differential_expression.*`
- `downstream.*`

When adding a new method or analysis branch, update all of:

- [`scripts/00_load_config.R`](/home/cam/github/single-cell-pipeline/scripts/00_load_config.R)
- [`scripts/run_pipeline.R`](/home/cam/github/single-cell-pipeline/scripts/run_pipeline.R) if it is a new stage
- [`config.yaml`](/home/cam/github/single-cell-pipeline/config.yaml)
- this file if the data contract or stage model changed

## Operational Caveats

- Optional steps are controlled by config, but stale checkpoints still matter. If you disable a step after previously generating its checkpoint, downstream code may still see that file unless the stage explicitly prefers another path.
- `--step ...` still force-runs a named stage even if its config toggle is false.
- CHOIR, CellChat, Monocle3, CytoTRACE2, and Tangram are not installed automatically by the baseline requirements script.
- Tangram requires Python-side dependencies and an external spatial query input.
- This repo is script-driven; verify changes with parse checks and targeted stage runs rather than assuming a test suite exists.

## Verification

Useful checks after edits:

- `Rscript -e "for (f in Sys.glob('scripts/*.R')) parse(file=f)"`
- `python3 -m py_compile scripts/tangram_annotate.py`
- `Rscript scripts/run_pipeline.R --step clustering --config <cfg>`
- `Rscript scripts/run_pipeline.R --step annotation --config <cfg>`
