# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import os
import sys
import argparse
import warnings


import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv


def fix_loom_barcodes(obs_names: np.ndarray) -> np.ndarray:
    fixed = np.array([str(bc).split(":")[-1] for bc in obs_names], dtype=object)

    fixed2 = []
    for bc in fixed:
        if bc.endswith("x"):
            fixed2.append(bc[:-1] + "-1")
        else:
            fixed2.append(bc)
    return np.array(fixed2, dtype=object)


def main():
    parser = argparse.ArgumentParser(
        description="Run scVelo velocity on CITE-seq (loom + Seurat-exported UMAP/cluster CSV) and save plots."
    )
    parser.add_argument("--loom", required=True, help="Path to velocyto .loom file")
    parser.add_argument("--umap_csv", required=True, help="Path to UMAP CSV exported from Seurat (row names = barcodes)")
    parser.add_argument("--clusters_csv", required=True, help="Path to clusters CSV exported from Seurat (row names = barcodes)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--plot", default="velocity_umap_clusters.png", help="Output plot filename")
    parser.add_argument("--mode", default="stochastic", choices=["stochastic", "deterministic", "dynamical"],
                        help="Velocity mode")
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--cluster_col", default="seurat_clusters",
                        help="Column name in clusters_csv to use (default: seurat_clusters)")
    parser.add_argument("--save_h5ad", default=None,
                        help="If set, save merged velocity AnnData to this filename (in outdir)")
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=DeprecationWarning)

    os.makedirs(args.outdir, exist_ok=True)
    scv.settings.figdir = args.outdir
    scv.settings.verbosity = 3
    scv.settings.set_figure_params(figsize=(6, 6), dpi=120)

    # 1) Load loom
    print(f"[INFO] Reading loom: {args.loom}")
    adata_velo = sc.read_loom(args.loom, sparse=True)

    # 2) Fix barcodes
    adata_velo.obs_names = fix_loom_barcodes(adata_velo.obs_names)
    if "CellID" in adata_velo.obs.columns:
        adata_velo.obs["CellID"] = adata_velo.obs_names

    # 3) Load UMAP / clusters CSV
    print(f"[INFO] Reading UMAP CSV: {args.umap_csv}")
    umap = pd.read_csv(args.umap_csv, index_col=0)

    print(f"[INFO] Reading clusters CSV: {args.clusters_csv}")
    clu = pd.read_csv(args.clusters_csv, index_col=0)

    if args.cluster_col not in clu.columns:
        print(f"[ERROR] cluster_col '{args.cluster_col}' not found in clusters_csv columns: {list(clu.columns)}", file=sys.stderr)
        sys.exit(10)

    # UMAP check
    if umap.shape[1] < 2:
        print(f"[ERROR] UMAP CSV seems to have <2 columns: shape={umap.shape}", file=sys.stderr)
        sys.exit(11)

    # 4) Intersect cells (loom vs CSV index)
    common = np.intersect1d(adata_velo.obs_names, umap.index.values)
    common = np.intersect1d(common, clu.index.values)
    print(f"[INFO] Common cells: {len(common)}")

    if len(common) == 0:
        print("[ERROR] No common cells after barcode fixing.", file=sys.stderr)
        print(f"[DEBUG] loom example: {adata_velo.obs_names[:5]}", file=sys.stderr)
        print(f"[DEBUG] umap index example: {umap.index.values[:5]}", file=sys.stderr)
        print(f"[DEBUG] clusters index example: {clu.index.values[:5]}", file=sys.stderr)
        sys.exit(2)

    # 5) Subset and attach UMAP/cluster
    adata_velo = adata_velo[common].copy()

    # UMAP save
    adata_velo.obsm["X_umap"] = umap.loc[common, umap.columns[:2]].values

    # cluster
    adata_velo.obs["seurat_clusters"] = pd.Categorical(clu.loc[common, args.cluster_col].astype(str).values)

    # 6) Velocity
    print("[INFO] Preprocessing & velocity estimation...")
    scv.pp.filter_and_normalize(
        adata_velo,
    )
    scv.pp.moments(adata_velo, n_pcs=args.n_pcs)

    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata_velo)
        scv.tl.velocity(adata_velo, mode="dynamical")
    else:
        scv.tl.velocity(adata_velo, mode=args.mode)

    scv.tl.velocity_graph(adata_velo)

    # 7) Plot (save only)
    print("[INFO] Saving plot...")
    scv.pl.velocity_embedding_stream(
        adata_velo,
        basis="umap",
        color="seurat_clusters",
        legend_loc="right",
        save=args.plot,
        show=False
    )
    plot_path = os.path.join(args.outdir, "scvelo_" + args.plot) if not args.plot.startswith("scvelo_") else os.path.join(args.outdir, args.plot)
    print(f"[INFO] Done. Plot saved under figdir: {args.outdir}")
    print(f"[INFO] (Note) scVelo typically prefixes 'scvelo_' to filenames. Look for: {plot_path}")

    # 8) Save merged AnnData (optional)
    if args.save_h5ad:
        out_h5ad = os.path.join(args.outdir, args.save_h5ad)
        print(f"[INFO] Saving AnnData: {out_h5ad}")
        adata_velo.write(out_h5ad)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
