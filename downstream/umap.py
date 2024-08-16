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

def plot_umap(hs_combined, filename):
    fig, ax = plt.subplots(figsize = (10,10))
    centroids = {}

    cluster_names = hs_combined.obs["final_clustering"].unique()

    for cluster in cluster_names:
        temp = hs_combined[hs_combined.obs["final_clustering"] == cluster,]
        centroids[cluster] = temp.obsm['X_umap_harmony'].mean(axis = 0)
        
    plt.tight_layout()
    uc = hs_combined.obsm['X_umap_harmony']
    sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = hs_combined.obs["final_clustering"], legend = "full", s = 5, ax = ax)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 5)
    ax.set_title(filename)

    for cluster, (x_umap, y_umap) in centroids.items():
        plt.text(x_umap, y_umap, cluster, fontsize=12, ha='center')

    fig.savefig(filename + '.png', dpi = 300, bbox_inches='tight')
    plt.close()

def plot_umap_cell_type(hs_combined, celltype, filename):
    fig, ax = plt.subplots(figsize = (10,10))
        
    plt.tight_layout()
    uc = hs_combined.obsm['X_umap_harmony']

    temp = hs_combined[hs_combined.obs["final_clustering"] == celltype,:]
    uct = temp.obsm['X_umap_harmony']

    sns.scatterplot(x = uc[:,0], y = uc[:,1], color = 'grey', alpha = 0.5, s = 5, ax = ax)
    sns.scatterplot(x = uct[:,0], y = uct[:,1], color = 'yellow', s = 5, ax = ax)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 5)
    ax.set_title(filename)

    fig.savefig(filename + '.png', dpi = 300, bbox_inches='tight')
    plt.close()

    
def plot_gene_umap(hs_combined, gene: str, filename):
    fig, ax = plt.subplots(figsize = (10,10))
    gene = gene.upper()
    temp = hs_combined[:, hs_combined.var['symbol'] == gene]
    
    # counts data is already scaled to zero mean, move to zero min, val is sd + min(x)
    temp.X = temp.X - min(temp.X)
        
    uc = hs_combined.obsm['X_umap']
    sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = temp.X.T[0], legend = "full", s = 5, ax = ax)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 5)
    ax.set_title(filename)


    fig.savefig(filename + '.png', dpi = 300, bbox_inches='tight')
    plt.close()

def plot_gene_umap_sc(adata, gene: str):
    # assume you've stored it there
    gene = gene.upper()
    centroids = adata.uns['centroids']
    rc_params = {
        'figure.dpi': 300,
        'figure.figsize': (7, 7),
        'figure.facecolor': 'white'
    }
    with plt.rc_context(rc = rc_params):
        sc.settings.set_figure_params(dpi = 300, frameon=False, figsize = (7,7), facecolor = "white")
        sc.pl.umap(adata, color = gene, show = False, size = 5, title = f"{gene} Expression HS UMAP")
        for cluster, (x_umap, y_umap) in centroids.items():
            plt.text(x_umap, y_umap, cluster, fontsize=12, ha='center', color='grey')
        plt.savefig(f"{gene}_HS.png", dpi = 300, bbox_inches = 'tight')

    
