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
df = pd.read_excel("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/sample_list.xlsx")

hc_gsm = df.loc[df['Using'],["GSE","GSM"]]

hc_gs_zip = list(zip(list(hc_gsm.to_dict()['GSE'].values()), list(hc_gsm.to_dict()['GSM'].values())))

print(hc_gs_zip)
print(len(hc_gs_zip))
# read in symbols for GSE173706

#symbols = pd.read_csv("symbol.csv", index_col = 0)
#symbols = symbols["symbol"]


# reads in a dir from GSM, returns adata
def read_dir(gse_code, gsm_code):

    data_path = gsm_path(gsm_code[3:])

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


def before_scrublet_step(adata):
    # preprocess healthy data

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


def scplot(func, *args, figsize = (7,7), path = "", **kwargs):
    rc_params = {
        'figure.dpi': 300,
        'figure.figsize': (7, 7),
        'figure.facecolor': 'white'
        }
    show = True
    with plt.rc_context(rc = rc_params):
        has_title = False
        if "title" in kwargs:
            has_title = True
            title = kwargs["title"]
            del kwargs["title"]
            show = False

        sc.settings.set_figure_params(dpi = 300, frameon=False, figsize = figsize, facecolor = "white")
        func(*args, **kwargs, show = show)
        
        if has_title:
            plt.title(title)
            if path:
                plt.savefig(os.path.join(path, f"{title}.png"), dpi = 300, bbox_inches = 'tight')
            else:
                plt.savefig(f"{title}.png", dpi = 300, bbox_inches = 'tight')
                
from scipy.interpolate import interp1d

multiplet_rate = np.array([0.40, 0.80, 1.60, 2.30, 3.10, 3.90, 4.60, 5.40, 6.10, 6.90, 7.60])/100
cells_loaded = np.array([800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000])

multiplet_rate_interp = interp1d(cells_loaded, multiplet_rate, kind = 'linear')


for GSE, GSM in tqdm(hc_gs_zip):
    try:
        path = os.path.join("preprocessing",GSE,GSM)
    
        if not os.path.exists(path):
            os.makedirs(path)
            
        # does first part of qc
        adata = read_dir(*(GSE, GSM))
    
        adata.write_h5ad(os.path.join(path, f"{GSM}_raw.h5ad"))
    
        preqc_cells = adata.n_obs
    
        adata = before_scrublet_step(adata)
    
        scplot(sc.pl.violin, adata, ["pct_counts_mt","pct_counts_ribo","pct_counts_hb"], title = "GSM5277170_violin_qc", path = path)
    
        # set expected doublet_rate accordingly
        doublet_rate = multiplet_rate_interp(adata.n_obs)
    
    
        scrub = scr.Scrublet(adata.X, expected_doublet_rate = doublet_rate)
        adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets(min_gene_variability_pctl=85, n_prin_comps=30)
        scrub.plot_histogram()
        plt.savefig(os.path.join(path,f"{GSM} scrublet hist.png"))
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
        scrub.plot_embedding('UMAP', order_points=True)
        plt.savefig(os.path.join(path,f"{GSM} scrublet umap.png"))
    
        adata = adata[~adata.obs["predicted_doublets"]]
    
        # post scrublet steps
        sc.pp.filter_cells(adata, max_genes = 6000)
    
        # remove cells with rna counts > 40000
        adata = adata[adata.obs["total_counts"] < 40000]
    
        postqc_cells = adata.n_obs
    
        pd.Series({"doublet rate" : doublet_rate, "preqc cells" : preqc_cells, "postqc cells" : postqc_cells}).to_csv(os.path.join(path, "args.csv"), header = False)
    
        adata.write_h5ad(os.path.join(path, f"{GSM}_pp.h5ad"))
    except Exception as e:
        print(e)
        print("\n")
        print("ERROR: ", GSE, GSM)

