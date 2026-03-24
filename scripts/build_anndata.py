#!/usr/bin/env python3
"""Build an AnnData (.h5ad) object from a Seurat intermediate directory.

Called by the R conversion wrapper via system2(). Reads the intermediate
directory format produced by seurat_to_intermediate() and writes a
standards-compliant .h5ad file.

AnnData mapping conventions:
  - adata.X          = raw counts (default assay)
  - adata.layers['counts'] = raw counts (explicit copy)
  - adata.layers['data']   = normalized data (if available)
  - adata.layers['<assay>_<layer>'] = extra assay layers
  - adata.obs        = cell metadata
  - adata.var        = feature metadata
  - adata.var['highly_variable'] = variable features flag
  - adata.obsm['X_<red>']  = reduction embeddings
  - adata.varm['<red>_loadings'] = reduction loadings (if full-rank)
  - adata.uns['<red>']     = stdev and partial loadings
  - adata.obsp['<name>']   = graph objects
  - adata.uns['seurat_conversion'] = provenance metadata
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


def load_mtx(path):
    """Load a Matrix Market file as CSC sparse matrix."""
    return sp.csc_matrix(sio.mmread(path))


def build_anndata(intermediate_dir, output_path):
    """Build AnnData from intermediate directory and write to .h5ad."""

    manifest_path = os.path.join(intermediate_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        print(f"ERROR: manifest.json not found in {intermediate_dir}", file=sys.stderr)
        sys.exit(1)

    with open(manifest_path) as f:
        manifest = json.load(f)

    default_assay = manifest["default_assay"]
    exported_assays = manifest["exported_assays"]

    # ---- Cell metadata (obs) ----
    obs = pd.read_csv(os.path.join(intermediate_dir, "obs.csv"))
    obs = obs.set_index("cell_id")
    obs.index.name = None

    # ---- Default assay ----
    assay_dir = os.path.join(intermediate_dir, f"assay_{default_assay}")
    if not os.path.isdir(assay_dir):
        print(f"ERROR: assay directory not found: {assay_dir}", file=sys.stderr)
        sys.exit(1)

    # Feature names
    features_path = os.path.join(assay_dir, "features.tsv")
    with open(features_path) as f:
        features = [line.strip() for line in f if line.strip()]

    # Feature metadata (var)
    var_path = os.path.join(assay_dir, "var.csv")
    if os.path.exists(var_path):
        var = pd.read_csv(var_path)
        var = var.set_index("feature_id")
        var.index.name = None
    else:
        var = pd.DataFrame(index=features)

    # Variable features -> highly_variable flag
    vf_path = os.path.join(assay_dir, "var_features.txt")
    if os.path.exists(vf_path):
        with open(vf_path) as f:
            var_features = {line.strip() for line in f if line.strip()}
        var["highly_variable"] = var.index.isin(var_features)
    else:
        var["highly_variable"] = False

    # Load layers from manifest
    assay_info = manifest.get("assay_info", {}).get(default_assay, {})
    layer_names = list(assay_info.get("layers", {}).keys())

    X = None
    layers = {}

    for layer_name in layer_names:
        mtx_path = os.path.join(assay_dir, f"{layer_name}.mtx")
        if not os.path.exists(mtx_path):
            continue

        # Transpose: Seurat is features x cells, AnnData is cells x features
        mat = load_mtx(mtx_path).T.tocsr()

        if layer_name == "counts":
            X = mat
            layers["counts"] = mat.copy()
        else:
            layers[layer_name] = mat

    if X is None:
        if layers:
            # Use the first available layer as X
            first_name = next(iter(layers))
            X = layers[first_name]
            print(f"  No counts layer found, using '{first_name}' as X")
        else:
            print("ERROR: no matrix data found for default assay", file=sys.stderr)
            sys.exit(1)

    # Build the AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers)
    print(f"  Built AnnData: {adata.n_obs} cells x {adata.n_vars} features")
    print(f"  Layers: {list(adata.layers.keys())}")

    # ---- Additional assays -> prefixed layers ----
    for assay_name in exported_assays:
        if assay_name == default_assay:
            continue

        extra_dir = os.path.join(intermediate_dir, f"assay_{assay_name}")
        if not os.path.isdir(extra_dir):
            continue

        extra_info = manifest.get("assay_info", {}).get(assay_name, {})
        extra_layers = list(extra_info.get("layers", {}).keys())

        for layer_name in extra_layers:
            mtx_path = os.path.join(extra_dir, f"{layer_name}.mtx")
            if not os.path.exists(mtx_path):
                continue

            mat = load_mtx(mtx_path).T.tocsr()

            # Dimensions may differ from default assay (different feature set)
            # Only add if cell count matches; skip if feature count differs
            if mat.shape[0] == adata.n_obs:
                prefixed = f"{assay_name}_{layer_name}"
                # If feature count matches, store directly in layers
                if mat.shape[1] == adata.n_vars:
                    adata.layers[prefixed] = mat
                    print(f"  Added layer: {prefixed}")
                else:
                    # Store in uns as we can't add mismatched dimensions to layers
                    print(f"  Skipping layer {prefixed}: "
                          f"feature count {mat.shape[1]} != {adata.n_vars}")

    # ---- Reductions -> obsm / varm / uns ----
    reductions = manifest.get("reductions", {})
    for red_name, red_info in reductions.items():
        red_dir = os.path.join(intermediate_dir, "reductions", red_name)
        emb_path = os.path.join(red_dir, "embeddings.csv")
        if not os.path.exists(emb_path):
            continue

        emb = pd.read_csv(emb_path, index_col=0)
        obsm_key = f"X_{red_name}"
        adata.obsm[obsm_key] = emb.values.astype(np.float32)
        print(f"  Added obsm['{obsm_key}']: {emb.shape}")

        # Loadings -> varm (if feature dimension matches) or uns
        load_path = os.path.join(red_dir, "loadings.csv")
        if os.path.exists(load_path):
            loadings = pd.read_csv(load_path, index_col=0)
            if loadings.shape[0] == adata.n_vars:
                varm_key = f"{red_name}_loadings"
                adata.varm[varm_key] = loadings.values.astype(np.float32)
                print(f"  Added varm['{varm_key}']: {loadings.shape}")
            else:
                # Partial loadings (e.g., only variable features) -> store in uns
                if red_name not in adata.uns:
                    adata.uns[red_name] = {}
                adata.uns[red_name]["loadings"] = loadings.values.astype(np.float32)
                adata.uns[red_name]["loadings_features"] = loadings.index.tolist()
                print(f"  Added uns['{red_name}']['loadings']: {loadings.shape} "
                      f"(partial, {loadings.shape[0]} features)")

        # Stdev -> uns
        stdev_path = os.path.join(red_dir, "stdev.txt")
        if os.path.exists(stdev_path):
            with open(stdev_path) as f:
                stdev = [float(x) for x in f.read().strip().split("\n") if x.strip()]
            if red_name not in adata.uns:
                adata.uns[red_name] = {}
            adata.uns[red_name]["stdev"] = np.array(stdev, dtype=np.float32)
            adata.uns[red_name]["variance"] = np.array(stdev, dtype=np.float32) ** 2

    # ---- Graphs -> obsp ----
    graphs_dir = os.path.join(intermediate_dir, "graphs")
    graph_manifest = manifest.get("graphs", {})
    if os.path.isdir(graphs_dir) and graph_manifest:
        for graph_name in graph_manifest:
            mtx_path = os.path.join(graphs_dir, f"{graph_name}.mtx")
            if os.path.exists(mtx_path):
                g = load_mtx(mtx_path).tocsr()
                adata.obsp[graph_name] = g
                print(f"  Added obsp['{graph_name}']: {g.shape}")

    # ---- Provenance metadata ----
    adata.uns["seurat_conversion"] = {
        "source": manifest.get("source", "seurat_v5"),
        "seurat_version": manifest.get("seurat_version", "unknown"),
        "default_assay": default_assay,
        "exported_assays": exported_assays,
        "project_name": manifest.get("project_name", ""),
        "format_version": manifest.get("format_version", "1.0"),
    }

    # Store factor column info for round-tripping as a JSON string.
    # h5ad mangles nested dicts with string arrays (especially single-element
    # arrays get split into characters), so we serialize as JSON to preserve
    # the structure exactly.
    factor_columns = manifest.get("factor_columns", {})
    if factor_columns:
        adata.uns["seurat_factor_columns_json"] = json.dumps(factor_columns)

    # ---- Write .h5ad ----
    adata.write_h5ad(output_path)

    print(f"\nAnnData written: {output_path}")
    print(f"  Shape: {adata.n_obs} x {adata.n_vars}")
    print(f"  Layers: {list(adata.layers.keys())}")
    print(f"  Obsm: {list(adata.obsm.keys())}")
    if adata.obsp:
        print(f"  Obsp: {list(adata.obsp.keys())}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build AnnData (.h5ad) from Seurat intermediate directory"
    )
    parser.add_argument("--input", required=True,
                        help="Path to intermediate directory")
    parser.add_argument("--output", required=True,
                        help="Output .h5ad file path")
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: intermediate directory not found: {args.input}",
              file=sys.stderr)
        sys.exit(1)

    build_anndata(args.input, args.output)
