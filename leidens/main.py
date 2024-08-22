df = pd.read_excel('mlist.xlsx', sheet_name = 1, header = 1)
df_dict = {col: df[col].dropna().tolist() for col in df.columns}

# plot dotplot

with plt.rc_context():
    matplotlib.rcParams['patch.edgecolor'] = 'black'
    fig, ax = plt.subplots(figsize = (60,15))
    sc.pl.dotplot(adata, marker_genes_final, groupby="leiden_2-0", standard_scale="var", show = False, ax = ax)
    plt.tight_layout()
    fig.savefig("leiden_rank_genes_2-0.png", bbox_inches = "tight", dpi = 300)
    plt.close()

# plot leiden
with plt.rc_context():
    matplotlib.rcParams['patch.edgecolor'] = 'black'
    fig, ax = plt.subplots(figsize = (7,7))
    sc.pl.umap(adata, color = "leiden_2-0", legend_loc = "on data", size = 3, show = False, ax = ax)
    plt.tight_layout()
    fig.savefig("leiden_2-0.png", bbox_inches = "tight", dpi = 300)
    plt.close()


marker_genes_final = {}

for cell_type, genes in marker_genes.items():
    first = True
    for gene in genes:
        if gene in adata.var.symbol:
            if first:
                marker_genes_final[cell_type] = [gene]
                first = False
            else:
                marker_genes_final[cell_type].append(gene)

adata.obs.loc[adata.obs["leiden_0-2"] == '9', 'major cell label'] = "Melanocyte"
adata.obs.loc[adata.obs["leiden_0-2"] == '10', 'major cell label'] = "Mast"
adata.obs.loc[adata.obs["leiden_0-2"] == '7', 'major cell label'] = "Plasma"
adata.obs.loc[adata.obs["leiden_0-2"] == '12', 'major cell label'] = "NK"
adata.obs.loc[adata.obs["leiden_0-2"] == '5', 'major cell label'] = "Fibroblast"
adata.obs.loc[adata.obs["leiden_0-2"] == '8', 'major cell label'] = "Neutrophil"


