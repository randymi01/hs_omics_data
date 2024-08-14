"""
1. Load Haniffa healthy, remove extra var cols
2. apply same qc procedures (not doublet since seems to be none)
3. reduce dim of healthy
4. find common genes with combined_hs then run pca, neighbors and umap on common features
5. ingest into combined_hs
6. generate plots

"""
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


healthy = sc.read_h5ad("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/healthy_dl_raw.h5ad")
healthy.var['symbol'] = healthy.var.index
healthy.var['gene_id'] = healthy.var['gene_ids-SKN8090524']
healthy.var.index = healthy.var['gene_id']
healthy.var = healthy.var[['symbol']]

# apply QC
# ---------------------------------------------------------------
healthy.var["mt"] = healthy.var["symbol"].str.startswith("MT-")

# ribosomal genes
healthy.var["ribo"] = healthy.var["symbol"].str.startswith(("RPS","RPL"))

# hemoglobin genes
healthy.var["hb"] = healthy.var["symbol"].str.contains("^HB[^(P)]")
    
# calculate qc metrics
sc.pp.calculate_qc_metrics(
    healthy,
    qc_vars = ["mt","ribo","hb"],
    inplace = True
)

# filter cells and genes
sc.pp.filter_cells(healthy, min_genes = 300)
sc.pp.filter_cells(healthy, max_genes = 6000)

# remove cells with mitocondrial dna > 20%
healthy = healthy[healthy.obs["pct_counts_mt"] < 20]

# remove cells with ribosomal dna > 55%
healthy = healthy[healthy.obs["pct_counts_ribo"] < 55]

# remove cells with rna counts > 40000
healthy = healthy[healthy.obs["total_counts"] < 40000]

# remove genes with less than 3 cells
sc.pp.filter_genes(healthy, min_cells = 3)

# remove doublets using scanpy.pp.scrublet
# not enough memory for healthy
# sc.pp.scrublet(healthy, verbose=True)

# reduce dim of healthy
# --------------------------------------------------------------------
sc.pp.normalize_total(healthy)
sc.pp.log1p(healthy)

# pick 5000 here since haniffa is 500k cells
sc.pp.highly_variable_genes(
    healthy,
    flavor = "cell_ranger",
    n_top_genes = 5000,
)

healthy = healthy[:, healthy.var['highly_variable']]

sc.pp.regress_out(healthy, ["pct_counts_mt"])

sc.pp.scale(healthy, max_value = 10)

# Label transfer
# --------------------------------------------------------------------------
adata = sc.read_h5ad("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/hs_data_expanded_2/hvg_r_2_umap.h5ad")
common_genes = adata.var.index.intersection(healthy.var.index)

# find common genes with combined_data
adata = adata[:, common_genes]
healthy = healthy[:, common_genes]

# calculate pca of healthy
sc.pp.pca(healthy)
sc.pp.neighbors(healthy)
sc.tl.umap(healthy)

# ingest
sc.tl.ingest(adata, healthy, obs="final_clustering")

# plot clustering
folder_name = 'v1'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# add centroid clusters
def plot_umap(hs_combined, filename):
    fig, ax = plt.subplots(figsize = (10,10))
    centroids = {}

    cluster_names = hs_combined.obs["final_clustering"].unique()

    for cluster in cluster_names:
        temp = hs_combined[hs_combined.obs["final_clustering"] == cluster,]
        centroids[cluster] = temp.obsm['X_umap'].mean(axis = 0)
        
    plt.tight_layout()
    uc = hs_combined.obsm['X_umap']
    sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = hs_combined.obs["final_clustering"], legend = "full", s = 5, ax = ax)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 5)
    ax.set_title(filename)

    for cluster, (x_umap, y_umap) in centroids.items():
        plt.text(x_umap, y_umap, cluster, fontsize=12, ha='center')

    fig.savefig(filename + '.png', dpi = 300, bbox_inches='tight')
    plt.close()

# make plots
plot_umap(adata, folder_name + '/' + "Label_Transfer_New_Data_2_Aligned_Preprocessing")
plot_umap(healthy, folder_name + '/' + "Healthy_Label_Clustering_Aligned")

# add back in dropped vars from larger complete dataset


# reasign native gene names
df = pd.read_csv("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/gene_names_native.csv")
df.index = df["ENSEMBL_ID"]
adata.var['symbol'] = df.loc[adata.var.index, ]["SYMBOL"]

# write aligned healthy and adata
adata.write_h5ad(folder_name + '/' + 'label_transfer_combined_v1.h5ad')
healthy.write_h5ad(folder_name + '/' + 'healthy_v1.h5ad')