#!/usr/bin/env python3
"""
Run Tangram-based annotation against a spatial/Xenium query using a reference h5ad.

This helper is intentionally separate from the R pipeline because Tangram lives in
Python and is only meaningful when an external spatial query is available.
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def load_xenium(xenium_dir: Path):
    h5_path = xenium_dir / "cell_feature_matrix.h5"
    cells_path = xenium_dir / "cells.parquet"
    if not h5_path.exists():
        raise FileNotFoundError(f"Missing Xenium matrix: {h5_path}")
    if not cells_path.exists():
        raise FileNotFoundError(f"Missing Xenium cells table: {cells_path}")

    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()
    cells = pd.read_parquet(cells_path).set_index("cell_id")
    common = adata.obs_names.intersection(cells.index)
    adata = adata[common, :].copy()
    adata.obsm["spatial"] = cells.loc[common, ["x_centroid", "y_centroid"]].to_numpy()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    return adata


def run_tangram(ad_sc, ad_sp, celltype_col):
    import tangram as tg

    shared_genes = sorted(set(ad_sc.var_names).intersection(set(ad_sp.var_names)))
    if len(shared_genes) < 200:
        raise RuntimeError(f"Too few shared genes for Tangram: {len(shared_genes)}")

    ad_sc = ad_sc[:, shared_genes].copy()
    ad_sp = ad_sp[:, shared_genes].copy()

    if "counts" in ad_sc.layers:
        ad_sc.X = ad_sc.layers["counts"].copy()
    if "counts" in ad_sp.layers:
        ad_sp.X = ad_sp.layers["counts"].copy()

    sc.pp.normalize_total(ad_sc, target_sum=1e4)
    sc.pp.log1p(ad_sc)
    sc.pp.normalize_total(ad_sp, target_sum=1e4)
    sc.pp.log1p(ad_sp)

    tg.pp_adatas(ad_sc, ad_sp, genes=shared_genes)
    training_genes = ad_sc.uns.get("training_genes", [])
    if len(training_genes) < 200:
        raise RuntimeError(f"Too few Tangram training genes after preprocessing: {len(training_genes)}")

    ad_map = tg.map_cells_to_space(
        ad_sc,
        ad_sp,
        mode="clusters",
        cluster_label=celltype_col,
        density_prior="rna_count_based",
    )
    tg.project_cell_annotations(ad_map, ad_sp, annotation=celltype_col)
    probs = ad_sp.obsm["tangram_ct_pred"]
    ad_sp.obs["tangram_celltype"] = probs.idxmax(axis=1).values
    ad_sp.obs["tangram_confidence"] = probs.max(axis=1).values
    return ad_sp, ad_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ad_sc", required=True)
    parser.add_argument("--xenium_dir")
    parser.add_argument("--query_h5ad")
    parser.add_argument("--out_prefix", required=True)
    parser.add_argument("--celltype_col", required=True)
    parser.add_argument("--method", default="tangram")
    args = parser.parse_args()

    if bool(args.xenium_dir) == bool(args.query_h5ad):
        raise ValueError("Provide exactly one of --xenium_dir or --query_h5ad")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ad_sc = sc.read_h5ad(args.ad_sc)
    if args.celltype_col not in ad_sc.obs.columns:
        raise ValueError(f"Reference label column not found: {args.celltype_col}")

    if args.xenium_dir:
        ad_sp = load_xenium(Path(args.xenium_dir))
    else:
        ad_sp = sc.read_h5ad(args.query_h5ad)
        if "counts" not in ad_sp.layers:
            ad_sp.layers["counts"] = ad_sp.X.copy()

    ad_sp, ad_map = run_tangram(ad_sc, ad_sp, args.celltype_col)

    obs_out = out_prefix.with_suffix(".tangram_predictions.csv")
    adata_out = out_prefix.with_suffix(".tangram_annotated.h5ad")
    map_out = out_prefix.with_suffix(".tangram_map.h5ad")

    ad_sp.obs.to_csv(obs_out)
    ad_sp.write_h5ad(adata_out)
    ad_map.write_h5ad(map_out)
    print(f"Saved: {obs_out}")
    print(f"Saved: {adata_out}")
    print(f"Saved: {map_out}")


if __name__ == "__main__":
    main()
