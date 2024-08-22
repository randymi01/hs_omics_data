# test reading in 10x data with AnnData
import pandas as pd
import numpy as np
import scanpy as sc
from glob import glob
import os
import scvi
import pickle
from scipy import sparse
import anndata as ad
import itertools

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("Spectral"))
sns.set_style("whitegrid")

from matplotlib.ticker import StrMethodFormatter
# ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))

import requests
from tqdm import tqdm

data_dir = "data"

sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=250, frameon=False, figsize=(5, 3), facecolor="white")

def run_umap(path):
    adata = sc.read_h5ad(path)

    sc.pp.normalize_total(adata)
    
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        flavor = "seurat",
        n_top_genes = 3000,
        batch_key = "GSM" 
    )

    # move after highly variable?
    adata = adata[:, adata.var['highly_variable']]

    sc.pp.regress_out(adata, ["pct_counts_mt"])
        # results in genes with zero mean

    sc.pp.scale(adata, max_value = 10)

    sc.pp.pca(adata)

    # evaluate pca performance
    adata.uns['pca']['variance_ratio'].sum()
    sc.pl.pca_variance_ratio(adata)


    sc.external.pp.harmony_integrate(adata, key = "GSM", basis = 'X_pca')

    sc.pp.neighbors(adata, use_rep = "X_pca_harmony")

    # umap uses generated neighbors
    sc.tl.umap(adata)

    adata.write_h5ad(path[:-5] + '_harmony_umap.h5ad')

    return adata


def plot_umap(hs_combined, filename):
    fig, ax = plt.subplots(figsize = (10,10))
    uc = hs_combined.obsm['X_umap']
    sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = hs_combined.obs["batch"], legend = "full", s = 2, ax = ax)
    ax.set_title(filename)
    ax.legend(markerscale = 5)
    fig.savefig("figures" + '/' + filename + '.png', dpi = 300)


# alternate way to save plot
rc_params = {
'figure.dpi': 300,
'figure.figsize': (7, 7),
'figure.facecolor': 'white'
}
#with plt.rc_context(rc = rc_params):
#    sc.pl.func(show = False)
#    plt.savefig("path", bbox_inches = "tight")

# with plt.rc_context(rc = rc_params):
#     sc.pl.umap(adata, color = ["batch"], color_map = "viridis", show = False)
#     plt.legend(title = "GSE", loc = 'center left', bbox_to_anchor = (1,0.5))
#     plt.title("Harmony Batch Integration All HS Disease Samples")
#     plt.savefig("t.png", bbox_inches = "tight")


