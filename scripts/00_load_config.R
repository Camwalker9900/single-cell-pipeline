#!/usr/bin/env Rscript

if (!requireNamespace("yaml", quietly = TRUE)) stop("Install 'yaml': install.packages('yaml')")
library(yaml)

source("scripts/00_pipeline_utils.R")

CONFIG_PATH <- Sys.getenv("SCRNA_CONFIG", unset = "config.yaml")
if (!file.exists(CONFIG_PATH)) {
  stop("Config file not found: ", CONFIG_PATH,
       "\nRun from the project root, or set SCRNA_CONFIG env variable.")
}

cfg <- yaml::read_yaml(CONFIG_PATH)

project_name <- cfg$project$name %||% "scRNA_Project"

run_cellbender  <- isTRUE(cfg$steps$run_cellbender)
run_cell_cycle  <- isTRUE(cfg$steps$run_cell_cycle)
run_integration <- isTRUE(cfg$steps$run_integration)
run_annotation  <- isTRUE(cfg$steps$run_annotation)
run_de          <- isTRUE(cfg$steps$run_differential_expression)
run_downstream  <- isTRUE(cfg$steps$run_downstream_tools)

cellbender_dir  <- cfg$paths$cellbender_output_dir %||% "cellbender_outputs"
output_dir      <- cfg$paths$output_dir %||% "pipeline_outputs"
checkpoints_dir <- cfg$paths$checkpoints_dir %||% "checkpoints"
plots_dir       <- file.path(output_dir, "plots")
annotation_dir  <- file.path(output_dir, "annotation")
de_dir          <- file.path(output_dir, "differential_expression")
downstream_dir  <- file.path(output_dir, "downstream")
exports_dir     <- file.path(output_dir, "exports")

for (d in c(output_dir, checkpoints_dir, plots_dir, annotation_dir, de_dir, downstream_dir, exports_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

ckpt <- list(
  qc          = file.path(checkpoints_dir, "01_qc.rds"),
  cellbender  = file.path(checkpoints_dir, "02_cellbender.rds"),
  cell_cycle  = file.path(checkpoints_dir, "03_cell_cycle.rds"),
  norm_pca    = file.path(checkpoints_dir, "04_norm_pca.rds"),
  integration = file.path(checkpoints_dir, "05_integration.rds"),
  clustering  = file.path(checkpoints_dir, "06_clustering.rds"),
  annotation  = file.path(checkpoints_dir, "07_annotation.rds"),
  de          = file.path(checkpoints_dir, "08_differential_expression.rds"),
  downstream  = file.path(checkpoints_dir, "09_downstream_tools.rds")
)

sample_names <- sapply(cfg$samples, `[[`, "name")

cellranger_dirs <- setNames(
  lapply(cfg$samples, `[[`, "cellranger_path"),
  sample_names
)

cb_files <- setNames(
  lapply(cfg$samples, function(s) file.path(cellbender_dir, s$cellbender_h5)),
  sample_names
)

barcode_files <- setNames(
  lapply(cfg$samples, function(s) file.path(cellbender_dir, s$cellbender_csv)),
  sample_names
)

sample_metadata <- setNames(
  lapply(cfg$samples, `[[`, "metadata"),
  sample_names
)

all_meta_fields <- unique(unlist(lapply(sample_metadata, names)))

min_cells    <- cfg$seurat$min_cells %||% 3
min_features <- cfg$seurat$min_features %||% 200
mt_pattern   <- cfg$seurat$mt_pattern %||% "^MT-"

qc_min_features       <- cfg$qc$min_features
qc_max_features       <- cfg$qc$max_features
qc_min_counts         <- cfg$qc$min_counts
qc_max_counts         <- cfg$qc$max_counts
qc_max_mt_percent     <- cfg$qc$max_mt_percent
run_doublet_detection <- isTRUE(cfg$qc$run_doublet_detection)

n_variable_features <- cfg$analysis$n_variable_features %||% 3000
n_pca_components    <- cfg$analysis$n_pca_components %||% 30

harmony_group_vars <- cfg$integration$harmony_group_vars %||% "batch"

clustering_method        <- cfg$clustering$method %||% "seurat"
clustering_reduction     <- cfg$clustering$reduction %||% if (run_integration) "harmony" else "pca"
clustering_resolutions   <- cfg$clustering$resolutions %||% c(0.2, 0.4, 0.6, 0.8, 1.0)
selected_resolution      <- cfg$clustering$selected_resolution %||% clustering_resolutions[[1]]
clustering_dims          <- cfg$clustering$dims %||% n_pca_components
choir_alpha              <- cfg$clustering$choir$alpha %||% 0.05
choir_n_cores            <- cfg$clustering$choir$n_cores %||% 1
choir_cluster_column     <- paste0("CHOIR_clusters_", choir_alpha)

annotation_cfg             <- cfg$annotation %||% list()
annotation_methods         <- annotation_cfg$methods %||% list()
annotation_ident_col       <- annotation_cfg$ident_col
annotation_findall_cfg     <- annotation_cfg$findallmarkers %||% list()
annotation_module_cfg      <- annotation_cfg$module_scores %||% list()
annotation_label_cfg       <- annotation_cfg$label_transfer %||% list()
annotation_tangram_cfg     <- annotation_cfg$tangram %||% list()

annotation_gene_sets       <- annotation_module_cfg$gene_sets
annotation_gene_set_file   <- annotation_module_cfg$gene_set_file

annotation_reference_rds   <- annotation_label_cfg$reference_rds
annotation_reference_assay <- annotation_label_cfg$reference_assay %||% "RNA"
annotation_reference_label <- annotation_label_cfg$reference_label_col
annotation_prediction_col  <- annotation_label_cfg$output_label_col %||% "predicted_label"
annotation_score_col       <- annotation_label_cfg$output_score_col %||% "prediction_score"
annotation_transfer_dims   <- annotation_label_cfg$dims %||% seq_len(min(30, n_pca_components))

annotation_tangram_enabled    <- isTRUE(annotation_methods$tangram)
annotation_tangram_query_type <- annotation_tangram_cfg$query_type %||% "xenium"
annotation_tangram_query_path <- annotation_tangram_cfg$query_path
annotation_tangram_python     <- annotation_tangram_cfg$python_bin %||% "python3"
annotation_tangram_script     <- annotation_tangram_cfg$script_path %||% "scripts/tangram_annotate.py"
annotation_tangram_label_col  <- annotation_tangram_cfg$label_col %||% annotation_reference_label %||% annotation_prediction_col
annotation_tangram_ref_h5ad   <- annotation_tangram_cfg$reference_h5ad

de_cfg          <- cfg$differential_expression %||% list()
de_comparisons  <- de_cfg$comparisons %||% list()

downstream_cfg       <- cfg$downstream %||% list()
cellchat_cfg         <- downstream_cfg$cellchat %||% list()
slingshot_cfg        <- downstream_cfg$slingshot %||% list()
monocle_cfg          <- downstream_cfg$monocle %||% list()
cytotrace_cfg        <- downstream_cfg$cytotrace %||% list()
milo_cfg             <- downstream_cfg$milo %||% list()

cat("Config loaded:", CONFIG_PATH, "\n")
cat("Project:", project_name, "\n")
cat("Samples:", paste(sample_names, collapse = ", "), "\n")
cat("Steps — CellBender:", run_cellbender,
    "| Cell Cycle:", run_cell_cycle,
    "| Integration:", run_integration,
    "| Annotation:", run_annotation,
    "| DE:", run_de,
    "| Downstream:", run_downstream, "\n")
cat("Clustering method:", clustering_method, "| reduction:", clustering_reduction, "\n\n")
