# scRNA-seq Analysis Pipeline

A modular, checkpoint-based single-cell RNA-seq pipeline built around Seurat. It supports:

- QC and optional CellBender filtering
- optional cell-cycle regression
- normalization, PCA, and optional Harmony integration
- clustering with either Seurat or CHOIR
- optional annotation methods
- optional differential expression
- optional downstream analysis hooks for CellChat, Slingshot, Monocle3, CytoTRACE2, and miloR
- optional export to AnnData (`.h5ad`)

The pipeline is configured through a single root-level file: `config.yaml`.

## Pipeline Steps

| # | Script | Required | Description |
|---|--------|----------|-------------|
| 01 | `scripts/01_qc.R` | Yes | Load per-sample `seurat_rds` inputs or Cell Ranger outputs, QC filtering, optional doublet detection, QC plots |
| 02 | `scripts/02_cellbender.R` | Optional | Filter QC-passed cells against CellBender barcodes |
| 03 | `scripts/03_cell_cycle.R` | Optional | Merge samples, score cell cycle, regress `S.Score` and `G2M.Score` |
| 04 | `scripts/04_normalization_pca.R` | Yes | Normalize, find variable features, scale, PCA |
| 05 | `scripts/05_integration.R` | Optional | Run Harmony integration |
| 06 | `scripts/06_clustering_umap.R` | Yes | Run UMAP and cluster with Seurat or CHOIR |
| 07 | `scripts/07_annotation.R` | Optional | Run `FindAllMarkers`, module scores, label transfer, and/or Tangram hook |
| 08 | `scripts/08_differential_expression.R` | Optional | Run configured differential-expression comparisons |
| 09 | `scripts/09_downstream_tools.R` | Optional | Run configured CellChat, Slingshot, Monocle3, CytoTRACE2, and/or miloR analyses |
| 10 | `scripts/10_export.R` | Optional | Export the best available object checkpoint to AnnData (`.h5ad`) |

## Quick Start

```bash
# 1. Clone
git clone https://github.com/Camwalker9900/single-cell-pipeline.git
cd scrna-pipeline

# 2. Install baseline dependencies
Rscript scripts/requirements.R

# 3. Edit the root config
nano config.yaml

# 4. Run the pipeline
Rscript scripts/run_pipeline.R
```

Replace the repository URL above with the real GitHub remote before publishing.

## Running The Pipeline

```bash
# Full pipeline (respects step toggles in config.yaml)
Rscript scripts/run_pipeline.R

# Run a single stage
Rscript scripts/run_pipeline.R --step qc
Rscript scripts/run_pipeline.R --step cellbender
Rscript scripts/run_pipeline.R --step cell_cycle
Rscript scripts/run_pipeline.R --step norm_pca
Rscript scripts/run_pipeline.R --step integration
Rscript scripts/run_pipeline.R --step clustering
Rscript scripts/run_pipeline.R --step annotation
Rscript scripts/run_pipeline.R --step differential_expression
Rscript scripts/run_pipeline.R --step downstream
Rscript scripts/run_pipeline.R --step export

# Resume from a stage onward
Rscript scripts/run_pipeline.R --from clustering
Rscript scripts/run_pipeline.R --from annotation

# Use a different config file
Rscript scripts/run_pipeline.R --config path/to/config.yaml
```

## Configuration

All project-specific behavior lives in `config.yaml`. For cluster runs, it is reasonable to keep additional project configs under `configs/` and submit batch scripts from `jobs/`.

Important sections:

- `project`: project name
- `steps`: enable or disable optional stages
- `paths`: output, checkpoint, and CellBender directories
- `samples`: one entry per sample, using either `seurat_rds` or `cellranger_path`, plus optional CellBender filenames and metadata
- `seurat`: object creation defaults
- `qc`: per-cell QC thresholds
- `analysis`: feature-selection and PCA settings
- `integration`: Harmony settings
- `clustering`: choose `method: seurat` or `method: choir`
- `annotation`: enable and configure marker discovery, module scores, label transfer, and Tangram
- `differential_expression`: explicit DE comparisons to run
- `downstream`: per-tool configs for CellChat, Slingshot, Monocle3, CytoTRACE2, and miloR
- `export`: optional `.h5ad` export settings

### Minimal Required Config

At minimum, fill in:

- `paths.cellbender_output_dir` if you will use CellBender
- `samples[*].seurat_rds` or `samples[*].cellranger_path`
- `samples[*].metadata`
- `clustering.reduction`
  Use `harmony` if integration is enabled and `pca` otherwise.

### Annotation Inputs

The annotation stage supports multiple independent methods:

- `findallmarkers`
  No extra input beyond a valid clustering identity column.
- `module_scores`
  Provide either inline `annotation.module_scores.gene_sets` or a `gene_set_file` in CSV, TSV, or YAML format.
- `label_transfer`
  Provide `annotation.label_transfer.reference_rds` and `reference_label_col`.
- `tangram`
  Provide an external spatial query through `annotation.tangram.query_path` plus Python dependencies. This is a spatial mapping hook, not a native scRNA annotation pass on the Seurat object itself.

### Differential Expression Inputs

DE does nothing unless you define comparisons under `differential_expression.comparisons`.

Example:

```yaml
differential_expression:
  comparisons:
    - name: "cluster0_control_vs_treatment"
      subset_col: "RNA_snn_res.0.4"
      subset_values: ["0"]
      group_by: "condition"
      ident_1: "control"
      ident_2: "treatment"
      test_use: "wilcox"
      min_pct: 0.1
      logfc_threshold: 0.25
```

## Checkpoints

Each stage writes an `.rds` checkpoint into `checkpoints/`.

```text
checkpoints/
├── 01_qc.rds
├── 02_cellbender.rds
├── 03_cell_cycle.rds
├── 04_norm_pca.rds
├── 05_integration.rds
├── 06_clustering.rds
├── 07_annotation.rds
├── 08_differential_expression.rds
└── 09_downstream_tools.rds
```

`10_export` does not create a new checkpoint; it writes files under `exports/`.

Important behavior:

- Early steps use a named list of per-sample Seurat objects.
- Step 3 or 4 merges samples into a single Seurat object.
- Downstream and export stages choose the best available upstream checkpoint by file existence.
- If you keep stale checkpoints from disabled stages, downstream behavior may still pick them up.

## Output Structure

Baseline outputs live under `pipeline_outputs/`.

Typical top-level structure:

```text
pipeline_outputs/
├── plots/
├── annotation/
├── differential_expression/
├── downstream/
├── exports/
├── qc_summary.csv
├── cellbender_filtering_summary.csv
├── analysis_summary.csv
├── cluster_summary.csv
├── final_object.rds
└── final_object_downstream.rds   # only if downstream stage runs
```

## Dependencies

`scripts/requirements.R` installs the baseline CRAN and Bioconductor packages used by the core pipeline.

These optional packages are not installed automatically there:

- `CHOIR`
- `CellChat`
- `monocle3`
- `SeuratWrappers`
- `CytoTRACE2`

You only need those if you enable the corresponding optional methods.

Tangram also requires Python-side dependencies such as `scanpy`, `pandas`, and `tangram`. The export step requires a working Python environment with `anndata` plus the helper scripts in `scripts/`.

## Verification

Useful low-cost checks after changes:

```bash
Rscript -e "for (f in Sys.glob('scripts/*.R')) parse(file=f)"
python3 -m py_compile scripts/tangram_annotate.py
Rscript scripts/run_pipeline.R --step clustering --config config.yaml
```

## Notes

- The pipeline assumes you run commands from the repository root.
- `--step ...` force-runs a named stage even if the corresponding config toggle is false.
- The safest way to add new methods is to update `config.yaml`, `scripts/00_load_config.R`, and `scripts/run_pipeline.R` together.
