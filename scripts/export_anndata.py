#!/usr/bin/env python3
"""Export an AnnData (.h5ad) object to the intermediate directory format.

Called by the R conversion wrapper via system2(). Writes the intermediate
directory structure that intermediate_to_seurat() can read to reconstruct
a Seurat v5 object.

Matrices are transposed from AnnData (cells x features) to Seurat
orientation (features x cells) before writing as Matrix Market files.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad


def save_mtx(mat, path):
    """Save a sparse matrix in Matrix Market format (features x cells)."""
    if not sp.issparse(mat):
        mat = sp.csc_matrix(mat)
    sio.mmwrite(path, mat.tocsc())


def export_anndata(h5ad_path, out_dir):
    """Export AnnData to intermediate directory format."""

    adata = ad.read_h5ad(h5ad_path)
    print(f"Loaded AnnData: {adata.n_obs} cells x {adata.n_vars} features")

    if os.path.isdir(out_dir):
        import shutil
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Determine assay name from provenance or default to RNA
    assay_name = "RNA"
    if "seurat_conversion" in adata.uns:
        assay_name = adata.uns["seurat_conversion"].get("default_assay", "RNA")

    # ---- Cell metadata (obs) ----
    obs = adata.obs.copy()
    obs.insert(0, "cell_id", obs.index)
    obs.to_csv(os.path.join(out_dir, "obs.csv"), index=False)
    print(f"  Exported obs: {obs.shape[0]} cells, {obs.shape[1] - 1} columns")

    barcodes = adata.obs_names.tolist()
    features = adata.var_names.tolist()

    # ---- Default assay ----
    assay_dir = os.path.join(out_dir, f"assay_{assay_name}")
    os.makedirs(assay_dir, exist_ok=True)

    # Feature and barcode names
    with open(os.path.join(assay_dir, "features.tsv"), "w") as f:
        f.write("\n".join(features))
    with open(os.path.join(assay_dir, "barcodes.tsv"), "w") as f:
        f.write("\n".join(barcodes))

    # Feature metadata (var)
    var = adata.var.copy()
    var.insert(0, "feature_id", var.index)
    var.to_csv(os.path.join(assay_dir, "var.csv"), index=False)

    # Variable features
    if "highly_variable" in adata.var.columns:
        vf = adata.var_names[adata.var["highly_variable"]].tolist()
        if vf:
            with open(os.path.join(assay_dir, "var_features.txt"), "w") as f:
                f.write("\n".join(vf))
            print(f"  Variable features: {len(vf)}")

    # X -> counts (transpose to features x cells)
    X_t = adata.X.T if sp.issparse(adata.X) else sp.csc_matrix(adata.X.T)
    save_mtx(X_t, os.path.join(assay_dir, "counts.mtx"))
    layer_info = {
        "counts": {
            "n_features": X_t.shape[0],
            "n_cells": X_t.shape[1],
            "n_nonzero": int(X_t.nnz) if sp.issparse(X_t) else int(np.count_nonzero(X_t))
        }
    }
    print(f"  Exported X as counts: {X_t.shape[0]} x {X_t.shape[1]}")

    # Named layers that belong to the default assay
    for layer_name, layer_data in adata.layers.items():
        # Skip counts (already exported as X) and prefixed layers (other assays)
        if layer_name == "counts":
            continue
        if "_" in layer_name and layer_name.split("_")[0] != assay_name:
            # Might be a prefixed layer from another assay
            # Only export non-prefixed or same-assay layers
            parts = layer_name.split("_", 1)
            if parts[0] != assay_name:
                continue
            clean_name = parts[1]
        else:
            clean_name = layer_name

        mat_t = layer_data.T if sp.issparse(layer_data) else sp.csc_matrix(layer_data.T)
        if mat_t.shape[0] == len(features) and mat_t.shape[1] == len(barcodes):
            save_mtx(mat_t, os.path.join(assay_dir, f"{clean_name}.mtx"))
            layer_info[clean_name] = {
                "n_features": mat_t.shape[0],
                "n_cells": mat_t.shape[1],
                "n_nonzero": int(mat_t.nnz) if sp.issparse(mat_t) else 0
            }
            print(f"  Exported layer '{clean_name}': {mat_t.shape[0]} x {mat_t.shape[1]}")

    # ---- Reductions from obsm ----
    manifest_reds = {}
    red_dir = os.path.join(out_dir, "reductions")

    for key, emb in adata.obsm.items():
        red_name = key[2:] if key.startswith("X_") else key
        red_subdir = os.path.join(red_dir, red_name)
        os.makedirs(red_subdir, exist_ok=True)

        # Column names follow Seurat convention: KEY_1, KEY_2, ...
        col_prefix = red_name.upper() + "_"
        col_names = [f"{col_prefix}{i + 1}" for i in range(emb.shape[1])]
        emb_df = pd.DataFrame(emb, index=barcodes, columns=col_names)
        emb_df.to_csv(os.path.join(red_subdir, "embeddings.csv"))

        red_meta = {
            "assay": assay_name,
            "key": col_prefix,
            "n_dims": int(emb.shape[1]),
            "has_loadings": False,
            "has_stdev": False,
        }

        # Check for loadings in varm
        varm_key = f"{red_name}_loadings"
        if varm_key in adata.varm:
            loadings = adata.varm[varm_key]
            load_cols = [f"{col_prefix}{i + 1}" for i in range(loadings.shape[1])]
            load_df = pd.DataFrame(loadings, index=features, columns=load_cols)
            load_df.to_csv(os.path.join(red_subdir, "loadings.csv"))
            red_meta["has_loadings"] = True

        # Check for partial loadings in uns
        if red_name in adata.uns and "loadings" in adata.uns[red_name]:
            uns_load = adata.uns[red_name]["loadings"]
            uns_features = adata.uns[red_name].get("loadings_features", None)
            if uns_features is not None and len(uns_features) > 0:
                load_cols = [f"{col_prefix}{i + 1}" for i in range(uns_load.shape[1])]
                load_df = pd.DataFrame(uns_load, index=uns_features, columns=load_cols)
                load_df.to_csv(os.path.join(red_subdir, "loadings.csv"))
                red_meta["has_loadings"] = True

        # Stdev from uns
        if red_name in adata.uns and "stdev" in adata.uns[red_name]:
            stdev = adata.uns[red_name]["stdev"]
            with open(os.path.join(red_subdir, "stdev.txt"), "w") as f:
                f.write("\n".join(str(float(s)) for s in stdev))
            red_meta["has_stdev"] = True

        manifest_reds[red_name] = red_meta
        print(f"  Exported reduction '{red_name}': {emb.shape}")

    # ---- Graphs from obsp ----
    manifest_graphs = {}
    graphs_dir = os.path.join(out_dir, "graphs")

    for key, graph in adata.obsp.items():
        os.makedirs(graphs_dir, exist_ok=True)
        save_mtx(graph, os.path.join(graphs_dir, f"{key}.mtx"))
        manifest_graphs[key] = {"n_cells": int(graph.shape[0])}
        print(f"  Exported graph '{key}': {graph.shape}")

    # ---- Factor columns from uns ----
    # h5ad may mangle string arrays (especially single-element ones), so be defensive
    factor_columns = {}
    if "seurat_factor_columns" in adata.uns:
        raw = adata.uns["seurat_factor_columns"]
        for k in raw:
            v = raw[k]
            if isinstance(v, str):
                factor_columns[k] = [v]
            elif isinstance(v, np.ndarray):
                if v.ndim == 0:
                    factor_columns[k] = [str(v.item())]
                else:
                    factor_columns[k] = [str(x) for x in v]
            elif hasattr(v, "__iter__"):
                factor_columns[k] = [str(x) for x in v]
            else:
                factor_columns[k] = [str(v)]

    # ---- Manifest ----
    n_vf = int(adata.var.get("highly_variable", pd.Series(dtype=bool)).sum())
    manifest = {
        "format_version": "1.0",
        "source": "anndata",
        "default_assay": assay_name,
        "exported_assays": [assay_name],
        "n_cells": int(adata.n_obs),
        "project_name": "",
        "assay_info": {
            assay_name: {
                "class": "Assay5",
                "n_features": len(features),
                "layers": layer_info,
                "n_variable_features": n_vf,
            }
        },
        "reductions": manifest_reds,
        "graphs": manifest_graphs,
        "factor_columns": factor_columns,
    }

    # Carry over provenance if available
    if "seurat_conversion" in adata.uns:
        sc = adata.uns["seurat_conversion"]
        manifest["project_name"] = sc.get("project_name", "")

    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    print(f"\nIntermediate export complete: {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export AnnData (.h5ad) to Seurat intermediate directory"
    )
    parser.add_argument("--input", required=True,
                        help="Input .h5ad file path")
    parser.add_argument("--output", required=True,
                        help="Output intermediate directory path")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    export_anndata(args.input, args.output)
