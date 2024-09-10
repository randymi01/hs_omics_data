import requests as re
from pathlib import Path
from tqdm import tqdm
import click
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

# pd.DataFrame(adata.var.index, columns = ["symbol"]).to_csv("/share/studies/Dermatology_Data/HS_Data_Portal/scanpy_seq/healthy data/preprocessing/GSE249793/GSM7964425/GSM7964425_syd.csv", index = False)


@click.command()
@click.option('-c','--codes', required=True, help='file path to csv containing ENSG ensembl gene ids in a csv. Colname is id.')
@click.option('-o','--output_path', required=False, show_default = True, default = "symbols.csv", help='file path for where library should be output')
@click.option('-n', '--null', required = False, is_flag = True, default = False, help = "whether to replace unmatched or novel transcript symbols with null instead of gene-id")
def main(codes,output_path, null):
    gene_codes = pd.read_csv(codes)

    symbol_dicts = []

    # ensure gene ids are unique
    gene_ids = list(set(list(gene_codes.iloc[:,0])))

    # limit to ensembl api is 1000
    gene_ids_idx = list(range(0, len(gene_ids), 900)) + [len(gene_ids)]


    with ThreadPoolExecutor(max_workers=5) as executor:
        for i in tqdm(range(len(gene_ids_idx) - 1), desc = "chunk"):
            chunk = gene_ids[gene_ids_idx[i] : gene_ids_idx[i+1]]
            future = executor.submit(ensembl_request, chunk, null)
            symbol_dicts.append(future.result())

    symbol_dict = {}
    for d in symbol_dicts:
        symbol_dict.update(d)

    symbol_dict = pd.DataFrame(symbol_dict.values(), index = symbol_dict.keys())
    symbol_dict.columns = ["symbol"]
    symbol_dict.to_csv(output_path)
    



def ensembl_request(chunk, null):
    if chunk:
        server = "https://rest.ensembl.org"
        ext = "/lookup/id"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

        for attempt in range(2):
            output = {}
            r = re.post(server+ext, headers=headers, data='{ "ids" : [' + ', '.join(f'"{ensg}"' for ensg in chunk) + ' ] }')
            if r.status_code == 200:
                response_data = r.json()

                for key, value in response_data.items():
                    alt = key
                    try:
                        if null:
                            alt = None
                        val = value['display_name'] if 'display_name' in value else alt
                        output[key] = val
                    except:
                        output[key] = alt
                break

    return output


if __name__ == "__main__":
    main()