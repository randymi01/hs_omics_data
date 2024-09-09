# test reading in 10x data with AnnData
import pandas as pd
import numpy as np
import scanpy as sc
from glob import glob
import os
import pickle
from scipy import sparse
import anndata as ad
import scrublet as scr
from pathlib import Path

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("Spectral"))
sns.set_style("whitegrid")

from matplotlib.ticker import StrMethodFormatter
# ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))

import requests
from tqdm import tqdm

# data_dir = "C:/Users/randymi/Desktop/seuratwork/data"

data_dir = "/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/data/"

gsm_dict = {}

for root, folders, _ in os.walk(data_dir):
    for i in folders:
        if i.startswith("GSM"):
            gsm_dict[i.split("_")[0]] = str(Path(os.path.join(root,i)))

# linux has different delimiters:
def detect_delimiter(path):
    if '/' in path:
        return '/'
    elif '\\' in path:
        return '\\'
    else:
        return None  # No delimiter found

# returns full path from GSM number
def gsm_path(gsm):
    return gsm_dict["GSM" + str(gsm)]
df = pd.read_excel("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/sample_list_healthy.xlsx")

hc_gsm = df.loc[df['HC Data'],["GSE","GSM"]]

hc_gs_zip = list(zip(list(hc_gsm.to_dict()['GSE'].values()), list(hc_gsm.to_dict()['GSM'].values())))


# read in symbols for GSE173706

symbols = pd.read_csv("symbol.csv", index_col = 0)
symbols = symbols["symbol"]


# reads in a dir from GSM, returns adata
def read_dir(gse_code, gsm_code):

    data_path = gsm_path(gsm_code[3:])

    if gse_code == "GSE173706":
        adata = ad.read_csv(os.path.join(data_path, os.path.basename(data_path) + '.csv')).T
        adata.var["symbol"] = symbols[adata.var.index]
        adata.obs['cell_id'] = adata.obs.index.copy()

    else:
        matrix_path = os.path.join(data_path, 'matrix.mtx')
        barcodes_path = os.path.join(data_path, 'barcodes.tsv')
        features_path = os.path.join(data_path, 'features.tsv')
        
        adata = sc.read_mtx(matrix_path)
        adata_bc = pd.read_csv(barcodes_path, header=None, delimiter = "\t")
        adata_features = pd.read_csv(features_path, header=None, delimiter = '\t')
        
        adata = adata.T

        adata.obs['cell_id'] = adata_bc[0].values

        # gene ids
        adata.var['symbol'] = adata_features[1].values

        # make ensembl id the gene index since it is unique
        adata.var_names = adata_features[0].values

    # Each cell gets unique name
    adata.obs_names = [f"cell_{i:d}" for i in range(adata.n_obs)]
    
    adata.obs["GSM"] = gsm_code[3:]

    # GSE is batch key
    adata.obs["batch"] = gse_code[3:]
    
    return adata


def before_scrublet_step(gse, gsm):
    # preprocess healthy data

    adata = read_dir(gse, gsm)

    print(f"Preqc cells: {adata.n_obs}")

    # mitochondrial genes
    adata.var["mt"] = adata.var["symbol"].str.startswith("MT-")

    # ribosomal genes
    adata.var["ribo"] = adata.var["symbol"].str.startswith(("RPS","RPL"))

    # hemoglobin genes
    adata.var["hb"] = adata.var["symbol"].str.contains("^HB[^(P)]")

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(
    adata,
    qc_vars = ["mt","ribo","hb"],
    inplace = True,
    )

    # filter cells and genes
    sc.pp.filter_cells(adata, min_genes = 300)

    # remove cells with mitocondrial dna > 20%
    adata = adata[adata.obs["pct_counts_mt"] < 20]

    # remove cells with ribosomal dna > 55%
    adata = adata[adata.obs["pct_counts_ribo"] < 55]

    # remove genes expressed in less than 3 cells 
    sc.pp.filter_genes(adata, min_cells = 3)

    # remove genes with total counts under 3
    sc.pp.filter_genes(adata, min_counts = 3)

    return adata

def post_scrublet_step(adata):
    # moved these steps to after scrublet since they might be removing doublets which interferes with the expect doublet ratio
    sc.pp.filter_cells(adata, max_genes = 6000)

    # remove cells with rna counts > 40000
    adata = adata[adata.obs["total_counts"] < 40000]

