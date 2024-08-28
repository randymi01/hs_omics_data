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

    sc.pp.highly_variable_genes(
        adata,
        flavor = "seurat_v3",
        n_top_genes = 4000
    )

    sc.pp.normalize_total(adata)
    
    sc.pp.log1p(adata)

    # load the grch38 cell cycle genes pickle
    sc.tl.score_genes_cell_cycle(adata, s_genes = s_genes, g2m_genes = g2m_genes, use_raw = False)

    # move after highly variable?
    adata = adata[:, adata.var['highly_variable']]

    sc.pp.regress_out(adata, ["pct_counts_mt", "pct_counts_ribo", "S_score", "G2M_score"])

        # results in genes with zero mean

    sc.pp.scale(adata)

    sc.pp.pca(adata)

    # evaluate pca performance
    adata.uns['pca']['variance_ratio'].sum()
    sc.pl.pca_variance_ratio(adata)


    sc.external.pp.harmony_integrate(adata, key = "GSM", basis = 'X_pca')

    sc.pp.neighbors(adata, use_rep = "X_pca_harmony")

    # umap uses generated neighbors
    sc.tl.umap(adata)

    sc.tl.leiden(adata, 0.2, key_added = "leiden_0-2")

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


#### SNIPPETS ####

# read in marker list and find genes that exist within current list:
mlist = pd.read_excel('path-to-mlist.xlsx')
marker_genes = mlist[mlist.columns[mlist.columns.str.endswith('gene')]].to_dict(orient = "list")

marker_genes_exist = {}

for key, value in marker_genes.items():
    gene_list = []
    for gene in value:
            if gene in adata.var['symbol'].values:
                    gene_list.append(gene)
    marker_genes_exist[key] = gene_list
marker_genes_exist

# using rank genes
sc.tl.rank_genes_groups(adata, "leiden_0-2", key_added = "leiden_0-2_rank-genes", n_genes = 3)

with plt.rc_context():
    matplotlib.rcParams['patch.edgecolor'] = 'black'
    fig, ax = plt.subplots(figsize = (60,15))
    sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_0-2", key = "leiden_0-2_rank-genes", show = False, ax = ax)
    plt.tight_layout()
    fig.savefig("GSE154775/leiden_rank_genes_0-2-clust.png", bbox_inches = "tight", dpi = 300)
    plt.close()


# using marker list
with plt.rc_context():
    matplotlib.rcParams['patch.edgecolor'] = 'black'
    fig, ax = plt.subplots(figsize = (100,15))
    sc.pl.dotplot(adata, marker_genes_exist, groupby="leiden_0-2", standard_scale="var", show = False, ax = ax, use_raw = False)
    plt.tight_layout()
    fig.savefig("GSE154775/leiden_marker_genes_0-2-clust.png", bbox_inches = "tight", dpi = 300)
    plt.close()

# setting marker groups
adata.obs.loc[adata.obs['leiden_0-2'].isin(['1', '0', '4', '2', '15', '14', '13', '12', '3']), 'label'] = 'Keratinocyte'


