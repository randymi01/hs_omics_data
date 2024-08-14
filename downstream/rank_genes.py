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

adata = sc.read_h5ad("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/hs_data_expanded_2/label_transfer/v1/label_transfer_combined_v2.h5ad")

adata.var['gene-id'] = adata.var.index
adata.var.index = adata.var['symbol']

sc.tl.rank_genes(adata, 'final_clustering')

adata.write_h5ad("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/hs_data_expanded_2/label_transfer/v1/label_transfer_combined_v2_ranked_genes.h5ad")

rc_params = {
'figure.dpi': 300,
'figure.figsize': (7, 7),
'figure.facecolor': 'white'
}
with plt.rc_context(rc = rc_params):
    sc.pl.rank_genes_groups(adata, show = False)
    plt.savefig("rank_genes_hs.png", dpi = 300, bbox_inches = "tight")
