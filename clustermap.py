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
sc.settings.set_figure_params(dpi=200, frameon=False, figsize=(5, 3), facecolor="white")

"""
with open("combined_hs_counts.pkl", "rb") as file:
    hs_combined = pickle.load(file)

# select highly varied genes on normalized
# top genes that overlap
# seurat_v3 expects raw counts
sc.pp.highly_variable_genes(
    hs_combined,
    flavor = "seurat_v3",
    n_top_genes = 3000,
    batch_key = "batch",
    # inplace subset
    subset = True
)

# normalize data
sc.pp.normalize_total(hs_combined)
sc.pp.log1p(hs_combined)

sc.tl.pca(hs_combined)

# integrate
sc.external.pp.harmony_integrate(hs_combined, key = "batch", basis = 'X_pca')

sc.pp.neighbors(hs_combined, use_rep = "X_pca_harmony")

# umap uses generated neighbors
sc.tl.umap(hs_combined)

sc.tl.leiden(hs_combined)

sc.tl.louvain(hs_combined)

sc.pl.umap(hs_combined, color=["GSM", "leiden", "louvain"], save = True, wspace = 0.25)

with open("combined_hs_processed.pkl", "wb") as file:
    pickle.dump(hs_combined, file)
"""
def pca(hs_combined):
    sc.pp.highly_variable_genes(
        hs_combined,
        flavor = "seurat_v3",
        n_top_genes = 3000,
        batch_key = "batch",
        # inplace subset
        subset = True
    )
    # normalize data
    sc.pp.regress_out(hs_combined, ["pct_counts_mt"])
    # results in genes with zero mean

    sc.pp.scale(hs_combined, max_value = 10)
    #sc.pp.normalize_total(hs_combined)
    #sc.pp.log1p(hs_combined)
    sc.pp.pca(hs_combined)


def make_umap(hs_combined, **kwargs):
    # harmony integrate
    sc.external.pp.harmony_integrate(hs_combined, key = "batch", basis = 'X_pca', **kwargs)

    # mnn integrate
    #
    #

    sc.pp.neighbors(hs_combined, use_rep = "X_pca_harmony")

    # umap uses generated neighbors
    sc.tl.umap(hs_combined)


def plot_umap(hs_combined, filename):
    fig, ax = plt.subplots(figsize = (10,10))
    uc = hs_combined.obsm['X_umap']
    sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = hs_combined.obs["batch"], legend = "full", s = 1, ax = ax)
    ax.set_title(filename)
    ax.legend(markerscale = 5)
    fig.savefig("figures" + '/' + filename + '.png', dpi = 300)

hyper_parameters = {'epsilon_harmony' : [1e-5,1e-6], 'lamb' : [0.75,0.5,0.25]}

def main():
    with open('combined_hs_no_cell_type.pkl', 'rb') as file:
        hs_combined = pickle.load(file)
    pca(hs_combined)
#    for values in itertools.product(*hyper_parameters.values()):
#        kwargs = dict(zip(hyper_parameters.keys(), values))
#        filename = f"umap_epsilon_{kwargs['epsilon_harmony']}_lambda_{kwargs['lamb']}"
#        make_umap(hs_combined, **kwargs)
#        plot_umap(hs_combined, filename)
