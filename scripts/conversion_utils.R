#!/usr/bin/env Rscript

# =============================================================================
# Seurat v5 <-> AnnData Conversion Utilities
#
# Core functions for converting between Seurat v5 objects and AnnData via
# an intermediate directory format. Matrices are stored as Matrix Market (.mtx),
# metadata as CSV, and a manifest.json tracks the structure.
#
# Intermediate directory layout:
#   manifest.json
#   obs.csv
#   assay_<name>/
#     features.tsv
#     barcodes.tsv
#     var.csv           (feature metadata)
#     var_features.txt   (variable features)
#     counts.mtx         (raw counts)
#     data.mtx           (normalized, if available)
#   reductions/
#     <name>/
#       embeddings.csv
#       loadings.csv     (if available)
#       stdev.txt        (if available)
#   graphs/
#     <name>.mtx
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# Seurat v5 -> Intermediate Directory
# -----------------------------------------------------------------------------
seurat_to_intermediate <- function(obj, out_dir,
                                   assays = NULL,
                                   join_layers = TRUE,
                                   include_graphs = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  library(Seurat)
  library(Matrix)

  if (is.null(assays)) assays <- DefaultAssay(obj)

  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)

  manifest <- list(
    format_version = "1.0",
    source = "seurat_v5",
    seurat_version = as.character(packageVersion("Seurat")),
    default_assay = DefaultAssay(obj),
    exported_assays = assays,
    n_cells = ncol(obj),
    project_name = obj@project.name
  )

  # ---- Cell metadata (obs) ----
  obs <- obj@meta.data
  obs$cell_id <- rownames(obs)
  obs <- obs[, c("cell_id", setdiff(colnames(obs), "cell_id")), drop = FALSE]

  # Active identity
  obs$active_ident <- as.character(Idents(obj))

  # Track factor columns for round-trip fidelity
  factor_cols <- list()
  for (col in colnames(obj@meta.data)) {
    if (is.factor(obj@meta.data[[col]])) {
      factor_cols[[col]] <- levels(obj@meta.data[[col]])
    }
  }
  manifest$factor_columns <- factor_cols

  write.csv(obs, file.path(out_dir, "obs.csv"), row.names = FALSE)
  cat("  Exported cell metadata:", ncol(obs) - 1, "columns,", nrow(obs), "cells\n")

  # ---- Per-assay export ----
  manifest$assay_info <- list()

  for (assay_name in assays) {
    cat("  Exporting assay:", assay_name, "\n")
    assay_dir <- file.path(out_dir, paste0("assay_", assay_name))
    dir.create(assay_dir, recursive = TRUE)

    assay_obj <- obj[[assay_name]]
    available_layers <- SeuratObject::Layers(assay_obj)

    # Join split layers if needed
    has_split <- any(grepl("\\.", available_layers))
    if (join_layers && has_split) {
      cat("    Joining split layers\n")
      obj <- JoinLayers(obj, assay = assay_name)
      assay_obj <- obj[[assay_name]]
      available_layers <- SeuratObject::Layers(assay_obj)
    }

    # Skip scale.data (dense, large, recomputable)
    export_layers <- setdiff(available_layers, "scale.data")

    # Feature and barcode names
    features <- rownames(assay_obj)
    barcodes <- colnames(assay_obj)
    writeLines(features, file.path(assay_dir, "features.tsv"))
    writeLines(barcodes, file.path(assay_dir, "barcodes.tsv"))

    # Feature metadata
    feat_meta <- tryCatch(assay_obj[[]], error = function(e) data.frame())
    if (!is.null(feat_meta) && ncol(feat_meta) > 0) {
      feat_meta$feature_id <- rownames(feat_meta)
      feat_meta <- feat_meta[, c("feature_id", setdiff(colnames(feat_meta), "feature_id")),
                             drop = FALSE]
      write.csv(feat_meta, file.path(assay_dir, "var.csv"), row.names = FALSE)
    }

    # Variable features
    vf <- VariableFeatures(obj, assay = assay_name)
    if (length(vf) > 0) {
      writeLines(vf, file.path(assay_dir, "var_features.txt"))
    }

    # Export each layer as MTX
    layer_info <- list()
    for (layer_name in export_layers) {
      mat <- LayerData(obj, assay = assay_name, layer = layer_name)
      if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")

      mtx_file <- file.path(assay_dir, paste0(layer_name, ".mtx"))
      Matrix::writeMM(mat, mtx_file)
      cat("    Layer", layer_name, ":", nrow(mat), "x", ncol(mat),
          "(", length(mat@x), "nonzero )\n")

      layer_info[[layer_name]] <- list(
        n_features = nrow(mat),
        n_cells = ncol(mat),
        n_nonzero = length(mat@x)
      )
    }

    manifest$assay_info[[assay_name]] <- list(
      class = class(assay_obj)[1],
      n_features = length(features),
      layers = layer_info,
      n_variable_features = length(vf)
    )
  }

  # ---- Reductions ----
  red_names <- names(obj@reductions)
  manifest$reductions <- list()

  if (length(red_names) > 0) {
    red_dir <- file.path(out_dir, "reductions")
    dir.create(red_dir, recursive = TRUE)

    for (red_name in red_names) {
      cat("  Exporting reduction:", red_name, "\n")
      red <- obj@reductions[[red_name]]
      red_subdir <- file.path(red_dir, red_name)
      dir.create(red_subdir, recursive = TRUE)

      # Embeddings (cells x dims)
      emb <- Embeddings(red)
      write.csv(emb, file.path(red_subdir, "embeddings.csv"))

      # Loadings (features x dims) — may not exist
      has_loadings <- FALSE
      tryCatch({
        load_mat <- Loadings(red)
        if (length(load_mat) > 0 && nrow(load_mat) > 0) {
          write.csv(load_mat, file.path(red_subdir, "loadings.csv"))
          has_loadings <- TRUE
        }
      }, error = function(e) NULL)

      # Standard deviations — may not exist
      has_stdev <- FALSE
      tryCatch({
        sdev <- Stdev(red)
        if (length(sdev) > 0) {
          writeLines(as.character(sdev), file.path(red_subdir, "stdev.txt"))
          has_stdev <- TRUE
        }
      }, error = function(e) NULL)

      manifest$reductions[[red_name]] <- list(
        assay = slot(red, "assay.used"),
        key = Key(red),
        n_dims = ncol(emb),
        has_loadings = has_loadings,
        has_stdev = has_stdev
      )
    }
  }

  # ---- Graphs (best-effort) ----
  graph_names <- names(obj@graphs)
  manifest$graphs <- list()

  if (include_graphs && length(graph_names) > 0) {
    graph_dir <- file.path(out_dir, "graphs")
    dir.create(graph_dir, recursive = TRUE)

    for (gn in graph_names) {
      tryCatch({
        g <- obj@graphs[[gn]]
        if (!inherits(g, "dgCMatrix")) g <- as(g, "dgCMatrix")
        Matrix::writeMM(g, file.path(graph_dir, paste0(gn, ".mtx")))
        cat("  Exported graph:", gn, "\n")
        manifest$graphs[[gn]] <- list(n_cells = nrow(g))
      }, error = function(e) {
        cat("  Warning: could not export graph", gn, ":", e$message, "\n")
      })
    }
  }

  # ---- Write manifest ----
  writeLines(
    jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    file.path(out_dir, "manifest.json")
  )

  cat("Intermediate export complete:", out_dir, "\n")
  invisible(manifest)
}


# -----------------------------------------------------------------------------
# Intermediate Directory -> Seurat v5
# -----------------------------------------------------------------------------
intermediate_to_seurat <- function(intermediate_dir) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  library(Seurat)
  library(Matrix)

  manifest <- jsonlite::fromJSON(
    file.path(intermediate_dir, "manifest.json"),
    simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE
  )

  default_assay <- manifest$default_assay
  assay_names <- manifest$exported_assays

  # ---- Cell metadata ----
  obs <- read.csv(file.path(intermediate_dir, "obs.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  rownames(obs) <- obs$cell_id
  cell_ids <- obs$cell_id
  obs$cell_id <- NULL

  # Extract active identity before removing from metadata
  active_ident <- NULL
  if ("active_ident" %in% colnames(obs)) {
    active_ident <- obs$active_ident
    names(active_ident) <- rownames(obs)
    obs$active_ident <- NULL
  }

  # Restore factor columns — defensive against mangled levels from h5ad round-trips
  factor_cols <- manifest$factor_columns
  if (!is.null(factor_cols)) {
    for (col in names(factor_cols)) {
      if (col %in% colnames(obs)) {
        lvls <- unlist(factor_cols[[col]])
        test <- factor(obs[[col]], levels = lvls)
        if (any(is.na(test) & !is.na(obs[[col]]))) {
          # Stored levels don't match data — fall back to unique values
          obs[[col]] <- factor(obs[[col]])
        } else {
          obs[[col]] <- test
        }
      }
    }
  }

  # ---- Build default assay ----
  cat("  Building assay:", default_assay, "\n")
  assay_dir <- file.path(intermediate_dir, paste0("assay_", default_assay))
  features <- readLines(file.path(assay_dir, "features.tsv"))

  # Counts (required)
  counts_path <- file.path(assay_dir, "counts.mtx")
  if (!file.exists(counts_path)) {
    mtx_files <- list.files(assay_dir, pattern = "\\.mtx$", full.names = TRUE)
    if (length(mtx_files) == 0) stop("No matrix files found for assay: ", default_assay)
    counts_path <- mtx_files[1]
    cat("    No counts.mtx found, using:", basename(counts_path), "\n")
  }

  counts <- readMM(counts_path)
  counts <- as(counts, "dgCMatrix")
  rownames(counts) <- features
  colnames(counts) <- cell_ids

  # Create Seurat object with counts
  project_name <- manifest$project_name %||% "ConvertedObject"
  obj <- CreateSeuratObject(counts = counts, assay = default_assay,
                            project = project_name)

  # Overwrite meta.data with the full obs (CreateSeuratObject adds its own columns)
  # Preserve the auto-computed columns, then merge with imported obs
  auto_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
  for (col in colnames(obs)) {
    if (!col %in% auto_cols) {
      obj@meta.data[[col]] <- obs[[col]]
    } else if (col %in% colnames(obs)) {
      obj@meta.data[[col]] <- obs[[col]]
    }
  }

  # Data layer (normalized) — if available
  data_path <- file.path(assay_dir, "data.mtx")
  if (file.exists(data_path)) {
    cat("    Adding data layer\n")
    data_mat <- readMM(data_path)
    data_mat <- as(data_mat, "dgCMatrix")
    rownames(data_mat) <- features
    colnames(data_mat) <- cell_ids
    LayerData(obj, assay = default_assay, layer = "data") <- data_mat
  }

  # Variable features
  vf_path <- file.path(assay_dir, "var_features.txt")
  if (file.exists(vf_path)) {
    vf <- readLines(vf_path)
    vf <- vf[vf != ""]
    VariableFeatures(obj, assay = default_assay) <- vf
    cat("    Restored", length(vf), "variable features\n")
  }

  # ---- Additional assays ----
  for (assay_name in setdiff(assay_names, default_assay)) {
    cat("  Building assay:", assay_name, "\n")
    extra_dir <- file.path(intermediate_dir, paste0("assay_", assay_name))
    if (!dir.exists(extra_dir)) {
      cat("    Warning: directory not found, skipping\n")
      next
    }

    extra_features <- readLines(file.path(extra_dir, "features.tsv"))
    extra_counts_path <- file.path(extra_dir, "counts.mtx")
    if (!file.exists(extra_counts_path)) {
      mtx_files <- list.files(extra_dir, pattern = "\\.mtx$", full.names = TRUE)
      if (length(mtx_files) == 0) next
      extra_counts_path <- mtx_files[1]
    }

    extra_counts <- readMM(extra_counts_path)
    extra_counts <- as(extra_counts, "dgCMatrix")
    rownames(extra_counts) <- extra_features
    colnames(extra_counts) <- cell_ids

    extra_assay <- CreateAssay5Object(counts = extra_counts)

    # Data layer
    extra_data_path <- file.path(extra_dir, "data.mtx")
    if (file.exists(extra_data_path)) {
      extra_data <- readMM(extra_data_path)
      extra_data <- as(extra_data, "dgCMatrix")
      rownames(extra_data) <- extra_features
      colnames(extra_data) <- cell_ids
      LayerData(extra_assay, layer = "data") <- extra_data
    }

    # Variable features
    extra_vf_path <- file.path(extra_dir, "var_features.txt")
    if (file.exists(extra_vf_path)) {
      VariableFeatures(extra_assay) <- readLines(extra_vf_path)
    }

    obj[[assay_name]] <- extra_assay
  }

  # ---- Reductions ----
  red_dir <- file.path(intermediate_dir, "reductions")
  if (dir.exists(red_dir) && !is.null(manifest$reductions)) {
    for (red_name in names(manifest$reductions)) {
      cat("  Restoring reduction:", red_name, "\n")
      red_subdir <- file.path(red_dir, red_name)
      emb_path <- file.path(red_subdir, "embeddings.csv")
      if (!file.exists(emb_path)) next

      emb <- as.matrix(read.csv(emb_path, row.names = 1, check.names = FALSE))

      red_info <- manifest$reductions[[red_name]]

      # Loadings
      load_mat <- matrix(nrow = 0, ncol = 0)
      load_path <- file.path(red_subdir, "loadings.csv")
      if (file.exists(load_path)) {
        load_mat <- as.matrix(read.csv(load_path, row.names = 1, check.names = FALSE))
      }

      # Stdev
      sdev <- numeric(0)
      stdev_path <- file.path(red_subdir, "stdev.txt")
      if (file.exists(stdev_path)) {
        sdev_raw <- readLines(stdev_path)
        sdev <- as.numeric(sdev_raw[sdev_raw != ""])
      }

      red_key <- red_info$key %||% paste0(toupper(red_name), "_")
      red_assay <- red_info$assay %||% default_assay

      # Build DimReduc — only pass loadings if non-empty
      if (nrow(load_mat) > 0 && ncol(load_mat) > 0) {
        red_obj <- CreateDimReducObject(
          embeddings = emb,
          loadings = load_mat,
          stdev = sdev,
          key = red_key,
          assay = red_assay
        )
      } else {
        red_obj <- CreateDimReducObject(
          embeddings = emb,
          stdev = sdev,
          key = red_key,
          assay = red_assay
        )
      }

      obj[[red_name]] <- red_obj
    }
  }

  # ---- Graphs (best-effort) ----
  graphs_dir <- file.path(intermediate_dir, "graphs")
  if (dir.exists(graphs_dir) && !is.null(manifest$graphs)) {
    for (gn in names(manifest$graphs)) {
      mtx_path <- file.path(graphs_dir, paste0(gn, ".mtx"))
      if (!file.exists(mtx_path)) next
      cat("  Restoring graph:", gn, "\n")
      tryCatch({
        g <- readMM(mtx_path)
        g <- as(g, "dgCMatrix")
        rownames(g) <- cell_ids
        colnames(g) <- cell_ids
        obj@graphs[[gn]] <- as.Graph(g)
      }, error = function(e) {
        cat("    Warning: could not restore graph", gn, ":", e$message, "\n")
      })
    }
  }

  # ---- Active identity ----
  if (!is.null(active_ident)) {
    Idents(obj) <- active_ident
  }

  cat("Seurat object rebuilt:", ncol(obj), "cells,", nrow(obj), "features\n")
  obj
}
