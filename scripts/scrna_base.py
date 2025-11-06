#!/usr/bin/env python3
import os
from os.path import join
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")  # ensure headless plotting
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import median_abs_deviation
from BCBio import GFF
import scanpy as sc


# ---------------------------
# Helpers
# ---------------------------
def get_protein_id(feature):
    try:
        protein_id = feature.sub_features[0].qualifiers.get("protein_id")
        return protein_id[0] if protein_id else ""
    except Exception:
        return ""

def get_product_name(feature):
    try:
        product = feature.sub_features[0].qualifiers.get("product")
        return product[0] if product else ""
    except Exception:
        return ""

def parse_gff_build_maps(gff_file):
    protein_ids, feature_ids, feature_types, product_names = [], [], [], []
    with open(gff_file) as in_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                if feature.type == "gene":
                    protein_ids.append(get_protein_id(feature))
                    product_names.append(get_product_name(feature))
                    feature_ids.append(feature.qualifiers["ID"][0])
                    feature_types.append(feature.qualifiers.get("gene_biotype", [""])[0])
    feature_types_map = dict(zip(feature_ids, feature_types))
    protein_ids_map = dict(zip(feature_ids, protein_ids))
    product_names_map = dict(zip(feature_ids, product_names))
    return feature_types_map, protein_ids_map, product_names_map

def calculate_qc_metrics(adata):
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["rRNA", "tRNA"], inplace=True, percent_top=[20], log1p=True
    )

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def top_genes_for_pc(loadings_col, genes, n=30):
    abs_vals = np.abs(loadings_col)
    idx = np.argsort(abs_vals)[::-1][:n]
    return pd.DataFrame({
        "gene": genes[idx],
        "loading": loadings_col[idx],
        "abs_loading": abs_vals[idx],
    })


# ---------------------------
# Main pipeline
# ---------------------------
def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.figdir = args.outdir  # where scanpy will save figures
    sc.settings.autoshow = False

    # --- GFF annotations ---
    feature_types_map, protein_ids_map, product_names_map = parse_gff_build_maps(args.gff)

    # --- Read data ---
    adata = sc.read_10x_mtx(args.input, var_names="gene_ids", cache=True)
    print(adata)

    # --- Annotate genes ---
    adata.var["gene_type"] = adata.var_names.map(feature_types_map)
    adata = adata[:, ~adata.var["gene_type"].isna()].copy()
    adata.var["protein_id"]    = adata.var_names.map(protein_ids_map)
    adata.var["product_names"] = adata.var_names.map(product_names_map)

    # Mark rRNA and tRNA genes
    adata.var["rRNA"] = adata.var["gene_type"].str.contains("rRNA", case=False, na=False)
    adata.var["tRNA"] = adata.var["gene_type"].str.contains("tRNA", case=False, na=False)

    # --- QC (initial) ---
    calculate_qc_metrics(adata)
    print(f"Total number of cells (raw): {adata.n_obs}")

    # --- Filter low counts ---
    adata_filter = adata.copy()
    sc.pp.filter_cells(adata_filter, min_counts=args.min_counts)
    print(f"Total number of cells (min_counts≥{args.min_counts}): {adata_filter.n_obs}")

    # --- Identify outliers ---
    # These columns exist after calculate_qc_metrics(log1p=True, percent_top=[20])
    adata_filter.obs["outlier"] = (
        is_outlier(adata_filter, "log1p_total_counts", 5)
        | is_outlier(adata_filter, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata_filter, "pct_counts_in_top_20_genes", 5)
        | is_outlier(adata_filter, "pct_counts_rRNA", 5)
    )
    print(adata_filter.obs["outlier"].value_counts())

    # --- Remove outliers and tRNA-high cells (if available) ---
    adata_filter = adata_filter[~adata_filter.obs["outlier"]].copy()
    if "pct_counts_tRNA" in adata_filter.obs.columns:
        adata_filter = adata_filter[adata_filter.obs["pct_counts_tRNA"] < 10].copy()

    # --- Recalculate QC metrics after filtering ---
    calculate_qc_metrics(adata_filter)

    # --- Remove rRNA and tRNA genes from matrix (keeps obs QC cols intact) ---
    adata_filter = adata_filter[:, ~(adata_filter.var["rRNA"] | adata_filter.var["tRNA"])].copy()
    print(f"Final shape after filtering: {adata_filter.shape}")

    # ======================================================
    # Layers: counts, norm, log1p
    # ======================================================
    adata_filter.layers["counts"] = adata_filter.X.copy()  # raw counts layer

    norm_res = sc.pp.normalize_total(adata_filter, target_sum=args.target_sum, inplace=False)
    adata_filter.layers["norm"] = norm_res["X"].copy()

    adata_filter.layers["log1p"] = np.log1p(adata_filter.layers["norm"])
    adata_filter.X = adata_filter.layers["log1p"].copy()   # use for downstream

    # ======================================================
    # PCA & plots
    # ======================================================
    # Scale (sparse-safe)
    if sp.issparse(adata_filter.X):
        sc.pp.scale(adata_filter, max_value=10, zero_center=False)
    else:
        sc.pp.scale(adata_filter, max_value=10)

    # PCA (match zero_center path)
    if sp.issparse(adata_filter.X):
        sc.tl.pca(adata_filter, svd_solver="arpack", zero_center=False)
    else:
        sc.tl.pca(adata_filter, svd_solver="arpack")

    # Save PCA variance ratio and PCA scatter
    sc.pl.pca_variance_ratio(adata_filter, log=True, show=False, save="_pca_variance.png")
    sc.pl.pca(
        adata_filter,
        color=["total_counts", "n_genes_by_counts", "pct_counts_rRNA", "pct_counts_tRNA"],
        show=False, save="_pca_scatter.png"
    )

    # Identify genes that impact PC1/PC2 the most & save
    loadings = adata_filter.varm["PCs"]  # (n_genes, n_pcs)
    genes = adata_filter.var_names
    top_pc1 = top_genes_for_pc(loadings[:, 0], genes, n=args.top_loading_genes)
    top_pc2 = top_genes_for_pc(loadings[:, 1], genes, n=args.top_loading_genes)
    top_pc1.to_csv(join(args.outdir, "top_genes_PC1.csv"), index=False)
    top_pc2.to_csv(join(args.outdir, "top_genes_PC2.csv"), index=False)

    # Bar plots for loadings
    sc.pl.pca_loadings(adata_filter, components=[1, 2], n_points=args.top_loading_genes,
                       show=False, save="_pca_loadings_1_2.png")

    # ======================================================
    # Neighbors / UMAP / Leiden
    # ======================================================
    sc.pp.neighbors(adata_filter, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.umap(adata_filter)
    sc.pl.umap(
        adata_filter,
        color=["total_counts", "n_genes_by_counts", "pct_counts_rRNA", "pct_counts_tRNA"],
        show=False, save="_umap_qc.png"
    )

    sc.tl.leiden(adata_filter, resolution=args.resolution)
    sc.pl.pca(adata_filter, color=["leiden"], show=False, save="_pca_leiden.png")
    sc.pl.umap(adata_filter, color=["leiden"], legend_loc="on data", show=False, save="_umap_leiden.png")

    # ======================================================
    # Differential expression & heatmaps
    # ======================================================
    sc.tl.rank_genes_groups(adata_filter, groupby="leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata_filter, n_genes=10, sharey=False, show=False, save="_rank_genes.png")

    # Dendrogram must match current leiden categories (recompute defensively)
    adata_filter.obs["leiden"] = adata_filter.obs["leiden"].astype("category")
    adata_filter.obs["leiden"] = adata_filter.obs["leiden"].cat.remove_unused_categories()
    adata_filter.uns.pop("dendrogram_leiden", None)
    sc.tl.dendrogram(adata_filter, groupby="leiden")

    sc.pl.rank_genes_groups_heatmap(
        adata_filter, n_genes=5, groupby="leiden", show=False, show_gene_labels=True,
        save="_rank_genes_heatmap.png"
    )

    # Also save the AnnData object
    adata_out = join(args.outdir, "adata_filtered.h5ad")
    adata_filter.write(adata_out)
    print(f"\n✓ Done. Results saved in: {args.outdir}")
    print(f"   - AnnData: {adata_out}")
    print("   - Figures: saved by Scanpy into that directory (filenames start with the suffixes given to `save=`).")
    print("   - Tables: top_genes_PC1.csv, top_genes_PC2.csv")


# ---------------------------
# CLI
# ---------------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Scanpy pipeline with GFF annotation, layers (counts/norm/log1p), PCA/UMAP/Leiden, and saved plots."
    )
    p.add_argument("--input", required=True,
                   help="Path to 10X mtx directory (folder containing matrix.mtx, barcodes.tsv, features.tsv).")
    p.add_argument("--gff", required=True, help="GFF annotation file.")
    p.add_argument("--outdir", required=True, help="Output directory for plots, tables, and h5ad.")
    p.add_argument("--min-counts", type=int, default=20, help="Minimum counts per cell for filtering.")
    p.add_argument("--target-sum", type=float, default=1e4, help="Library size target for normalize_total.")
    p.add_argument("--n-neighbors", type=int, default=10, help="Neighbors for graph construction.")
    p.add_argument("--n-pcs", type=int, default=30, help="Number of PCs for neighbors/UMAP.")
    p.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution.")
    p.add_argument("--top-loading-genes", type=int, default=30, help="Top genes to report for PC1/PC2 and loading bars.")
    args = p.parse_args()
    main(args)
