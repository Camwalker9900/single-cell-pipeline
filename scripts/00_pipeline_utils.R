#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

ensure_package <- function(pkg, hint = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Required package not installed: ", pkg)
    if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
    stop(msg, call. = FALSE)
  }
}

maybe_join_layers <- function(obj, assay = NULL) {
  assay <- assay %||% DefaultAssay(obj)
  if (!exists("JoinLayers", mode = "function")) return(obj)
  tryCatch(JoinLayers(obj, assay = assay), error = function(e) obj)
}

merge_sample_objects <- function(input_data, project_name) {
  if (!is.list(input_data) || inherits(input_data, "Seurat")) return(input_data)
  if (length(input_data) == 0) stop("Input sample list is empty", call. = FALSE)

  sample_list <- input_data
  for (i in seq_along(sample_list)) {
    sn <- names(sample_list)[i]
    sample_list[[i]] <- RenameCells(sample_list[[i]], new.names = paste0(sn, "_", colnames(sample_list[[i]])))
  }

  merge(sample_list[[1]], y = sample_list[-1], add.cell.ids = NULL, project = project_name)
}

read_best_checkpoint <- function(paths, missing_message = "No upstream checkpoint found.") {
  existing <- Filter(file.exists, paths)
  if (length(existing) == 0) stop(missing_message, call. = FALSE)
  readRDS(existing[[1]])
}

get_best_checkpoint_path <- function(paths, missing_message = "No upstream checkpoint found.") {
  existing <- Filter(file.exists, paths)
  if (length(existing) == 0) stop(missing_message, call. = FALSE)
  existing[[1]]
}

set_idents_from_column <- function(obj, ident_col) {
  if (is.null(ident_col)) return(obj)
  if (!ident_col %in% colnames(obj@meta.data)) {
    stop("Identity column not found in metadata: ", ident_col, call. = FALSE)
  }
  Idents(obj) <- obj[[ident_col, drop = TRUE]]
  obj
}

choose_default_ident_col <- function(obj, selected_resolution = NULL) {
  candidate_cols <- c(
    if (!is.null(selected_resolution)) paste0("RNA_snn_res.", selected_resolution),
    "CHOIR_cluster",
    "seurat_clusters"
  )
  candidate_cols <- unique(candidate_cols[!is.na(candidate_cols)])
  for (col in candidate_cols) {
    if (col %in% colnames(obj@meta.data)) return(col)
  }
  NULL
}

choose_first_present_col <- function(obj, candidates) {
  for (col in candidates) {
    if (!is.null(col) && !is.na(col) && col %in% colnames(obj@meta.data)) return(col)
  }
  NULL
}

read_gene_sets <- function(path = NULL, inline_sets = NULL) {
  if (!is.null(inline_sets)) {
    out <- inline_sets
    if (is.null(names(out)) || any(names(out) == "")) {
      stop("Inline gene sets must be a named list", call. = FALSE)
    }
    return(lapply(out, unique))
  }

  if (is.null(path) || identical(path, "")) return(list())
  if (!file.exists(path)) stop("Gene set file not found: ", path, call. = FALSE)

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yaml", "yml")) {
    ensure_package("yaml", "Install with install.packages('yaml')")
    cfg <- yaml::read_yaml(path)
    if (!is.list(cfg) || is.null(names(cfg))) {
      stop("YAML gene set file must be a named mapping of set -> genes", call. = FALSE)
    }
    return(lapply(cfg, function(x) unique(as.character(unlist(x)))))
  }

  sep <- if (ext == "tsv") "\t" else ","
  df <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
  set_col <- intersect(c("set", "cluster", "label", "name"), colnames(df))[1]
  gene_col <- intersect(c("gene", "feature", "marker"), colnames(df))[1]
  if (is.na(set_col) || is.na(gene_col)) {
    stop("Gene set table must contain columns like 'set' and 'gene'", call. = FALSE)
  }

  split(df[[gene_col]], df[[set_col]]) |> lapply(function(x) unique(as.character(stats::na.omit(x))))
}

save_top_markers <- function(markers, out_file, n = 10) {
  if (nrow(markers) == 0) return(invisible(NULL))
  gene_col <- intersect(c("gene", "Gene", "feature"), colnames(markers))[1]
  score_col <- intersect(c("avg_log2FC", "avg_logFC", "avg_diff"), colnames(markers))[1]
  cluster_col <- intersect(c("cluster", "ident"), colnames(markers))[1]
  if (is.na(gene_col) || is.na(score_col) || is.na(cluster_col)) return(invisible(NULL))

  top_markers <- do.call(
    rbind,
    lapply(split(markers, markers[[cluster_col]]), function(df) df[order(df[[score_col]], decreasing = TRUE), , drop = FALSE][seq_len(min(n, nrow(df))), , drop = FALSE])
  )
  utils::write.csv(top_markers, out_file, row.names = FALSE)
}

get_plot_reduction <- function(obj, preferred = c("umap", "harmony", "pca")) {
  for (red in preferred) {
    if (red %in% names(obj@reductions)) return(red)
  }
  NULL
}
