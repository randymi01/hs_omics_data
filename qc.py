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

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("Spectral"))
sns.set_style("whitegrid")

from matplotlib.ticker import StrMethodFormatter
# ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))

import requests
from tqdm import tqdm

data_dir = "data"

# Function is beyond broken
# scanpy.read_10x_mtx("data_extract/GSE154775/GSM4679492_HS1")

# GSE 154775 - is GRCH37 everything else is 38

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

# returns anndata object of 10x data
def read_dir(gsm_code, **kwargs):
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
    
    # kwargs contain sample wide data
    for key, value in kwargs.items():
        adata.obs[key.upper()] = value
    
    return adata


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


def preprocess(sample_code, **kwargs):
    HS1 = read_dir(sample_code, **kwargs)

    output_dir = "preprocessing/GSM{}".format(str(sample_code))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ## add gene labels to features

    # mitochondrial genes
    HS1.var["mt"] = HS1.var["symbol"].str.startswith("MT-")

    # ribosomal genes
    HS1.var["ribo"] = HS1.var["symbol"].str.startswith(("RPS","RPL"))

    # hemoglobin genes
    HS1.var["hb"] = HS1.var["symbol"].str.contains("^HB[^(P)]")

    # output cell counts 
    qc_counts = HS1.var[["mt","ribo","hb"]].apply(sum, axis = 0)
    qc_counts["total pre-qc cells"] = HS1.n_obs
    qc_counts["total pre-qc genes"] = HS1.n_vars
    

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(
    HS1,
    qc_vars = ["mt","ribo","hb"],
    inplace = True,
    )

    # generate qc violin plot
    ax = sc.pl.violin(
    HS1,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    jitter=0.4,
    multi_panel=True,
    show = False,
    ax = None
    )
    
    # Set the title for the FacetGrid
    ax.fig.suptitle("GSM {} (n_cells = {})".format(str(sample_code), HS1.n_obs), fontsize=16)  # Adjust fontsize as needed

    ax.axes.flat[0].set_title('Unique Genes Count')
    ax.axes.flat[1].set_title('Total Gene Count')
    ax.axes.flat[2].set_title('Percent Count in Mt Genes')
    ax.axes.flat[3].set_title('Percent Count in Rb Genes')
    ax.axes.flat[0].set_ylabel('')
    ax.axes.flat[1].set_ylabel('')
    ax.axes.flat[2].set_ylabel('')
    ax.axes.flat[3].set_ylabel('')
    # Adjust layout to make room for the title
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    ax.savefig(os.path.join(output_dir,"{}_pre-qc_violin.png".format(str(sample_code))))
    plt.close()


    ## filter cells and genes
    sc.pp.filter_cells(HS1, min_genes = 300)
    sc.pp.filter_cells(HS1, max_genes = 6000)

    # remove cells with mitocondrial dna > 20%
    HS1 = HS1[HS1.obs["pct_counts_mt"] < 20]

    # remove cells with ribosomal dna > 55%
    HS1 = HS1[HS1.obs["pct_counts_ribo"] < 55]

    # remove cells with rna counts > 40000
    HS1 = HS1[HS1.obs["total_counts"] < 40000]

    # remove genes with less than 3 cells
    sc.pp.filter_genes(HS1, min_cells = 3)

    # remove doublets using scanpy.pp.scrublet
    sc.pp.scrublet(HS1, verbose=True)

    print("Percentage Doublets Detected: {}".format(
        sum(HS1.obs["predicted_doublet"])/HS1.n_obs))

    HS1 = HS1[~HS1.obs["predicted_doublet"]]


    # normalize and find varied features after joining do not do it in qc

    # normalize and store raw counts
    # HS1.layers["counts"] = HS1.X.copy()

    # normalize to median counts then lognormalize x' = log(x+1)
    #sc.pp.normalize_total(HS1)
    #sc.pp.log1p(HS1)
    
    # now HS1.X is the lognormalized rna counts

    # get most varied features
    #sc.pp.highly_variable_genes(HS1, n_top_genes=2000, batch_key='GSM', flavor = "seurat_v3")

    #sc.pl.highly_variable_genes(HS1, show = False)
    #plt.savefig(os.path.join(output_dir,"{}_variable_features.png".format(str(sample_code))))

    # generate qc violin plot
    ax = sc.pl.violin(
    HS1,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    jitter=0.4,
    multi_panel=True,
    show = False,
    ax = None
    )
    
    # Set the title for the FacetGrid
    ax.fig.suptitle("GSM {} (n_cells = {})".format(str(sample_code), HS1.n_obs), fontsize=16)  # Adjust fontsize as needed

    ax.axes.flat[0].set_title('Unique Genes Count')
    ax.axes.flat[1].set_title('Total Gene Count')
    ax.axes.flat[2].set_title('Percent Count in Mt Genes')
    ax.axes.flat[3].set_title('Percent Count in Rb Genes')
    ax.axes.flat[0].set_ylabel('')
    ax.axes.flat[1].set_ylabel('')
    ax.axes.flat[2].set_ylabel('')
    ax.axes.flat[3].set_ylabel('')
    # Adjust layout to make room for the title
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    ax.savefig(os.path.join(output_dir,"{}_post-qc_violin.png".format(str(sample_code))))

    sc.pl.scatter(HS1, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show = False)
    plt.savefig(os.path.join(output_dir,"{}_post-qc_mt_count.png".format(str(sample_code))))
    plt.close()
    
    qc_counts["total post-qc cells"] = HS1.n_obs
    qc_counts["total post-qc genes"] = HS1.n_vars
    qc_counts.to_csv(os.path.join(output_dir,"qc_cell_counts.csv"),header=False)

    save_pickle(HS1, sample_code)

    # save gsm info as html
    #response = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM{}".format(str(sample_code)))

    # 
    #if response.status_code == 200:
    #    # Open a file in write mode
    #    with open(os.path.join(output_dir,'sample_info.html'), 'w', encoding='utf-8') as file:
    #        # Write the content of the response to the file
    #        file.write(response.text)
    #    print('Sample info saved to output.html')
    #else:
    #    print(f'Failed to retrieve the URL. Status code: {response.status_code}')
        
# get list of pickled GSMs
pickled_gsms = [int(i[:-4]) for i in os.listdir('preprocessing/pickles')]

print(pickled_gsms)

if __name__ == "__main__":
    with open("errors.txt","w") as ErrorFile:
        for sample in tqdm(list(gsm_dict.keys())):
            sample_code = int(sample[3:])
            if sample_code not in pickled_gsms:
                preprocess(sample_code)