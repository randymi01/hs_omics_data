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


def chunk_list(lst, n):
    return [lst[i:i + n] for i in range(0, len(lst), n)]

# takes list of gsm codes, loop merge pickles three at a time
def merge_anndata(codes):
    first_chunk = True
    chunk_size = 3
    for chunk in tqdm(chunk_list(codes, chunk_size)):
        adatas = []
        for code in chunk:
            with open('preprocessing/pickles/{}.pkl'.format(str(code)), 'rb') as File:
                adatas.append(pickle.load(File))
        if not first_chunk:
            adatas.append(combined_adata)
        
        combined_adata = ad.concat(
                adatas,
                axis=0,
                join='outer'
        )
        
        first_chunk = False
    return combined_adata

if __name__ == "__main__":
    # read from sample list
    df = pd.read_excel('sample_list.xlsx')
    HS_samples = df[df['Disease'] == 'HS']
    HS_samples = list(HS_samples['GSM'].apply(lambda x: int(x[3:])))

    combined_data = merge_anndata(HS_samples)
    combined_data.obs.index = ['cell_' + str(i) for i in range(combined_data.n_obs)]

    with open("combined_hs_counts.pkl", "wb") as File:
        pickle.dump(combined_data, File)


