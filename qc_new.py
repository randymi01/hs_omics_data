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

# make dictionary of gsm code and path
gsm_dict = {}

for root, folders, _ in os.walk(data_dir):
    for i in folders:
        if i.startswith("GSM"):
            gsm_dict[i.split("_")[0]] = os.path.join(root,i)

# linux has different delimiters:
def detect_delimiter(path):
    if '/' in path:
        return '/'
    elif '\\' in path:
        return '\\'
    else:
        return None  # No delimiter found


gse_dict = {key : value.split(detect_delimiter(value))[1] for key, value in gsm_dict.items()}

# returns full path from GSM number
def gsm_path(gsm):
    return gsm_dict["GSM" + str(gsm)]

# saves as pickle
def save_pickle(adata_obj, gsm_code):
    output_dir = "preprocessing/pickles"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
        with open("preprocessing/pickles/{}.pkl".format(gsm_code), "wb") as File:
            pickle.dump(adata_obj, File)
    except Exception as e:
        print(e)

# reads in a dir from GSM, returns adata
def read_dir(gsm_code):

    data_path = gsm_path(gsm_code)

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
    
    adata.obs["GSM"] = gsm_code

    # GSE is batch key
    adata.obs["batch"] = gse_dict["GSM"+str(gsm_code)][3:]
    
    return adata


# applys preprocessing filters
def preprocess(gsm):

    adata = read_dir(gsm)

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

    ## filter cells and genes
    sc.pp.filter_cells(adata, min_genes = 300)
    sc.pp.filter_cells(adata, max_genes = 6000)

    # remove cells with mitocondrial dna > 20%
    adata = adata[adata.obs["pct_counts_mt"] < 20]

    # remove cells with ribosomal dna > 55%
    adata = adata[adata.obs["pct_counts_ribo"] < 55]

    # remove cells with rna counts > 40000
    adata = adata[adata.obs["total_counts"] < 40000]

    # remove genes with less than 3 cells
    sc.pp.filter_genes(adata, min_cells = 3)

    # remove doublets using scanpy.pp.scrublet
    sc.pp.scrublet(adata, verbose=True)

    adata = adata[~adata.obs["predicted_doublet"]]

    save_pickle(adata, sample_code)

def preprocess_adata(adata):

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

    ## filter cells and genes
    sc.pp.filter_cells(adata, min_genes = 300)
    sc.pp.filter_cells(adata, max_genes = 6000)

    # remove cells with mitocondrial dna > 20%
    adata = adata[adata.obs["pct_counts_mt"] < 20]

    # remove cells with ribosomal dna > 55%
    adata = adata[adata.obs["pct_counts_ribo"] < 55]

    # remove cells with rna counts > 40000
    adata = adata[adata.obs["total_counts"] < 40000]

    # remove genes with less than 3 cells
    sc.pp.filter_genes(adata, min_cells = 3)

    # remove doublets using scanpy.pp.scrublet
    sc.pp.scrublet(adata, verbose=True)

    adata = adata[~adata.obs["predicted_doublet"]]

    return adata

data_path = 'data/haniffa_hs7'

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

adata.obs["GSM"] = 0
# GSE is batch key
adata.obs["batch"] = 0

adata = preprocess_adata(adata)

with open('preprocessing/pickles/haniffa.pkl', 'wb') as file:
    pickle.dump(adata, file)


if __name__ == "__main__":
    with open("errors.txt","w") as ErrorFile:
        for sample in tqdm(list(gsm_dict.keys())):
            sample_code = int(sample[3:])
            preprocess(sample_code)
            print(sample_code)