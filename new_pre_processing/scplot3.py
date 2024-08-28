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

import matplotlib

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("Spectral"))
sns.set_style("whitegrid")

from matplotlib.ticker import StrMethodFormatter
# ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))

import requests
from tqdm import tqdm

def plot_umap_clustering(hs_combined, filename):
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
    
def plot_umap(hs_combined, filename, colors, rep = 'X_umap'):
    uc = hs_combined.obsm[rep]
    
    if type(colors) == list:
      fig, axes = plt.subplots(1, len(colors), figsize = (5 * len(colors), 5))
      for ax, color in zip(axes,colors):
        sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = hs_combined.obs[color], legend = "full", s = 5, ax = ax)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 5)
        ax.set_title(filename + f': {color}')
    if type(colors) == str:
      fig, ax = plt.subplots(figsize = (10,10))
      sns.scatterplot(x = uc[:,0], y = uc[:,1], hue = hs_combined.obs[colors], legend = "full", s = 5, ax = ax)
      ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 5)
      ax.set_title(filename + f': {colors}')
    plt.tight_layout()
    fig.savefig(filename + '.png', dpi = 300, bbox_inches='tight')
    plt.close()

def plot_umap_sc(adata, color, filename = None, legend_loc = "on data"):
  with plt.rc_context():
      matplotlib.rcParams['patch.edgecolor'] = 'black'
      fig, ax = plt.subplots(figsize = (7,7))
      sc.pl.umap(adata, color = color, legend_loc = legend_loc, s = 3, show = False, ax = ax)
      plt.tight_layout()
      if filename:
        fig.savefig(f"{filename}.png", bbox_inches = "tight", dpi = 300)
      else:
        plt.show()
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

def plot_gene_umap_sc(adata, gene: str, show_centroids = False, folder_name: str = None, use_rep: str = None, show: bool = False):
    # assume centroids stored in adata.uns
    if use_rep:
        current_rep = adata.obsm["X_umap"].copy()
        adata.obsm["X_umap"] = adata.obsm[use_rep].copy()
    try:    
        gene = gene.upper()
        centroids = {}
        if show_centroids:
            centroids = adata.uns['centroids']
        with plt.rc_context():
            sc.settings.set_figure_params(dpi = 300, frameon=False, figsize = (7,7), facecolor = "white")
            sc.pl.umap(adata, color = gene, show = show, size = 5, title = f"{gene} Expression HS UMAP")
            for cluster, (x_umap, y_umap) in centroids.items():
                plt.text(x_umap, y_umap, cluster, fontsize=12, ha='center')
            if folder_name:
                plt.savefig(f"{folder_name}/{gene}_HS.png")
            else:
                plt.savefig(f"{gene}_HS.png")
    except:
        print("plotting error")
    finally:
        if use_rep:
            adata.obsm["X_umap"] = current_rep.copy()
       
def plot_genes_umap_sc(adata, genes: list[str], title = None, show_centroids = False, show = False):
    # assume you've stored it there
    genes = [gene.upper() for gene in genes]
    
    centroids = {}
    if show_centroids:
        centroids = adata.uns['centroids']
    
    with plt.rc_context():
        n_col = 3 if len(genes) >= 3 else len(genes)
        n_row = int(np.ceil(len(genes)/3))
        fig, axes = plt.subplots(n_row, n_col, figsize = (n_col * 2.5, n_row * 2.5))
        axes = axes.flatten()

        for i, gene in enumerate(genes):
            sc.pl.umap(adata, color = gene, show = False, size = 5, title = f"{gene} Expression HS UMAP", ax = axes[i])
            axes[i].set_title(axes[i].get_title(), fontsize = 10)
        
        for idx in range(len(genes), n_row*n_col):
            axes[idx].set_visible(False)

        for cluster, (x_umap, y_umap) in centroids.items():
            plt.text(x_umap, y_umap, cluster, fontsize=12, ha='center')
        plt.tight_layout()
        if title:
            fig.suptitle(title, fontsize=16, y=1.05, x = 0.0, ha='left')  # y adjusts the position of the suptitle
            plt.savefig(f"{title}_HS.png", bbox_inches='tight')
        if show:
          plt.show()

def sc_plot(scanpy_plotting_func, *args, filename = "", **kwargs):
    with plt.rc_context():
        matplotlib.rcParams['patch.edgecolor'] = 'black'
        fig, ax = plt.subplots(figsize = (7,7))
        scanpy_plotting_func(*args, **kwargs, ax = ax, show = False)
        if filename:
            ax.set_title(filename)
            plt.savefig(f"{filename}.png", bbox_inches = 'tight', dpi = 300)
        else:
            plt.show()
        plt.close()
    
