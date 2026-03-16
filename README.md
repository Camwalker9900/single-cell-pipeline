# scRNA-seq Analysis Pipeline

A modular, reusable single-cell RNA-seq analysis pipeline built on **Seurat**, **Harmony**, and **scDblFinder**. Clone this repo, fill in `config/config.yaml`, and run.

---

## Pipeline Steps

| # | Script | Required | Description |
|---|--------|----------|-------------|
| 01 | `01_qc.R` | ✅ | Load CellRanger outputs, QC filtering, doublet detection, pre/post plots |
| 02 | `02_cellbender.R` | optional | Ambient RNA removal using CellBender barcodes |
| 03 | `03_cell_cycle.R` | optional | Cell cycle scoring (S/G2M) + regression during scaling |
| 04 | `04_normalization_pca.R` | ✅ | Normalize, HVG selection, scale, PCA |
| 05 | `05_integration.R` | optional | Harmony batch correction |
| 06 | `06_clustering_umap.R` | ✅ | UMAP, multi-resolution clustering, plots, summaries |

Optional steps are toggled in `config/config.yaml` under `steps:`.

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/your-org/scrna-pipeline.git
cd scrna-pipeline

# 2. Install dependencies (once)
Rscript envs/requirements.R

# 3. Edit config for your project
#    — set your paths, samples, and metadata
nano config/config.yaml

# 4. Run
Rscript run_pipeline.R
```

---

## Configuration

Everything project-specific lives in **`config/config.yaml`**. There are four things to fill in:

### 1. Step toggles
```yaml
steps:
  run_cellbender:  true   # set false to skip
  run_cell_cycle:  true
  run_integration: true
```

### 2. Paths
```yaml
paths:
  cellbender_output_dir: "/data/my_project/cellbender"
  output_dir: "pipeline_outputs"
  checkpoints_dir: "checkpoints"
```

### 3. Samples — one entry per sample
```yaml
samples:
  - name: "ctrl_rep1"
    cellranger_path: "/data/cellranger/ctrl_rep1/outs/filtered_feature_bc_matrix"
    cellbender_h5:  "ctrl_rep1_output_filtered.h5"
    cellbender_csv: "ctrl_rep1_output_cell_barcodes.csv"
    metadata:
      batch: "batch01"
      condition: "control"
```

Metadata fields (e.g. `batch`, `condition`) are **completely flexible** — name them whatever fits your experiment. The pipeline uses them dynamically for plots, summaries, and integration grouping.

### 4. Clustering reduction
```yaml
clustering:
  reduction: "harmony"   # use "pca" if integration is skipped
```

---

## Running Individual Steps

```bash
# Run full pipeline (respects config toggles)
Rscript run_pipeline.R

# Run a single step only
Rscript run_pipeline.R --step qc
Rscript run_pipeline.R --step cellbender
Rscript run_pipeline.R --step cell_cycle
Rscript run_pipeline.R --step norm_pca
Rscript run_pipeline.R --step integration
Rscript run_pipeline.R --step clustering

# Resume from a step (if a previous step failed)
Rscript run_pipeline.R --from norm_pca

# Use a different config file
Rscript run_pipeline.R --config path/to/other_config.yaml
```

---

## Checkpoint System

Each step saves its output to `checkpoints/` as an `.rds` file. Every step auto-detects which upstream checkpoint is most recent, so **you can skip any optional step and the pipeline routes around it automatically**.

```
checkpoints/
├── 01_qc.rds
├── 02_cellbender.rds   (only if run)
├── 03_cell_cycle.rds   (only if run)
├── 04_norm_pca.rds
├── 05_integration.rds  (only if run)
└── 06_clustering.rds
```

---

## Output Structure

```
pipeline_outputs/
├── final_object.rds           # Final Seurat object
├── qc_summary.csv
├── cellbender_filtering_summary.csv
├── analysis_summary.csv
├── cluster_summary.csv
└── plots/
    ├── qc_violin_pre_filter.pdf
    ├── qc_scatter_pre_filter.pdf
    ├── qc_violin_post_filter.pdf
    ├── qc_scatter_post_filter.pdf
    ├── variable_features.pdf
    ├── cell_cycle_phases.pdf        (if Step 03 ran)
    ├── pca_plots.pdf
    ├── pca_elbow_plot.pdf
    ├── harmony_integration.pdf      (if Step 05 ran)
    ├── final_umap_visualizations.pdf
    ├── cell_cycle_final_umap.pdf    (if Step 03 ran)
    └── quality_control_metrics.pdf
```

---

## QC Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `qc.min_features` | 200 | Min genes per cell |
| `qc.max_features` | 6000 | Max genes per cell |
| `qc.min_counts` | 500 | Min UMIs per cell |
| `qc.max_counts` | 30000 | Max UMIs per cell |
| `qc.max_mt_percent` | 20 | Max mitochondrial % |
| `qc.run_doublet_detection` | true | Run scDblFinder |

Set any threshold to `null` to disable it.

---

## Adding More Steps

1. Create `modules/R/07_your_step.R`
2. Add `source("modules/R/00_load_config.R")` at the top — all config variables are immediately available
3. Load from `ckpt$clustering` (or whichever step precedes yours)
4. Save your output and register a new `ckpt$your_step` in `00_load_config.R`
5. Add an entry to `step_scripts` in `run_pipeline.R`
