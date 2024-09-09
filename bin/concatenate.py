#!/usr/bin/env python3

import json
from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime
from os import fspath, walk, listdir
from pathlib import Path
from typing import Dict, Tuple

import anndata
import muon as mu
import numpy as np
import os
import pandas as pd
import requests
import scipy.sparse
import uuid
import yaml

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / "data",
    Path("/opt/data"),
]


def get_tissue_type(dataset: str) -> str:
    organ_dict = yaml.load(open("/opt/organ_types.yaml"), Loader=yaml.BaseLoader)
    organ_code = requests.get(
        f"https://entity.api.hubmapconsortium.org/dataset/{dataset}/organs/"
    )
    organ_name = organ_dict[organ_code]
    return organ_name.replace(" (Left)", "").replace(" (Right)", "")


def convert_tissue_code(tissue_code):
    with open("/opt/organ_types.yaml", 'r') as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    tissue_name = data.get(tissue_code)['description']
    return tissue_name


def get_inverted_gene_dict():
    inverted_dict = defaultdict(list)
    gene_mapping = read_gene_mapping()
    for ensembl, hugo in gene_mapping.items():
        inverted_dict[hugo].append(ensembl)
    return inverted_dict


def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath
                

def find_file_pairs(directory):
    cell_by_bin_patterns = ["cell_by_bin.h5ad"]
    cell_by_gene_patterns = ["cell_by_gene.h5ad"]
    cell_by_bin_file = find_files(directory, cell_by_bin_patterns)
    cell_by_gene_file = find_files(directory, cell_by_gene_patterns)
    return (cell_by_bin_file, cell_by_gene_file)


def annotate_h5ads(
    adata_file, tissue_type: str, uuids_df: pd.DataFrame):
    # Get the directory
    data_set_dir = fspath(adata_file.parent.stem)
    # And the tissue type
    tissue_type = tissue_type if tissue_type else get_tissue_type(data_set_dir)
    dense_adata = anndata.read_h5ad(adata_file)
    adata = make_new_anndata_object(dense_adata)
    del dense_adata
    adata_copy = adata.copy()
    adata_copy.obs["barcode"] = adata.obs.index
    adata_copy.obs["dataset"] = data_set_dir
    
    cell_ids_list = [
        "-".join([data_set_dir, barcode]) for barcode in adata_copy.obs["barcode"]
    ]
    adata_copy.obs["cell_id"] = pd.Series(
        cell_ids_list, index=adata_copy.obs.index, dtype=str
    )
    adata_copy.obs.set_index("cell_id", drop=True, inplace=True)
    adata_copy = map_gene_ids(adata_copy)
    return adata_copy


def read_gene_mapping() -> Dict[str, str]:
    """
    Try to find the Ensembl to HUGO symbol mapping, with paths suitable
    for running this script inside and outside a Docker container.
    :return:
    """
    for directory in GENE_MAPPING_DIRECTORIES:
        mapping_file = directory / "ensembl_to_symbol.json"
        if mapping_file.is_file():
            with open(mapping_file) as f:
                return json.load(f)
    message_pieces = ["Couldn't find Ensembl → HUGO mapping file. Tried:"]
    message_pieces.extend(f"\t{path}" for path in GENE_MAPPING_DIRECTORIES)
    raise ValueError("\n".join(message_pieces))


def map_gene_ids(adata):
    obsm = adata.obsm
    uns = adata.uns
    gene_mapping = read_gene_mapping()
    has_hugo_symbol = [gene in gene_mapping for gene in adata.var.index]
    # adata = adata[:, has_hugo_symbol]
    temp_df = pd.DataFrame(
        adata.X.todense(), index=adata.obs.index, columns=adata.var.index
    )
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var["hugo_symbol"] = [
        gene_mapping.get(var, np.nan) for var in adata.var.index
    ]
    adata.obsm = obsm
    adata.uns = uns
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.var_names_make_unique()
    return adata


def create_json(tissue, data_product_uuid, creation_time, uuids, hbmids, cell_count, file_size):
    bucket_url = f"https://hubmap-data-products.s3.amazonaws.com/{data_product_uuid}/"
    metadata = {
        "Data Product UUID": data_product_uuid,
        "Tissue": convert_tissue_code(tissue),
        "Assay": "atac",
        "URL": bucket_url + f"{tissue}.h5mu",
        "Creation Time": creation_time,
        "Dataset UUIDs": uuids,
        "Dataset HBMIDs": hbmids,
        "Total Cell Count": cell_count,
        "Raw File Size": file_size
    }
    print("Writing metadata json")
    with open(f"{data_product_uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)


def make_mudata(cell_by_bin, cell_by_gene):
    mdata = mu.MuData({"atac_cell_by_bin": cell_by_bin, "atac_cell_by_gene": cell_by_gene})
    mu.pp.intersect_obs
    return mdata


def annotate_mudata(mdata, uuids_df):
    merged = uuids_df.merge(mdata.obs, left_on="uuid", right_on="dataset", how="inner")
    merged = merged.set_index(mdata.obs.index)
    merged = merged.drop(columns=["Unnamed: 0"])
    merged = merged.fillna(np.nan)
    merged["age"] = pd.to_numeric(merged["age"])
    return merged


def make_new_anndata_object(adata):
    new_adata = anndata.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs.index), var=adata.var)
    return new_adata


def main(data_directory: Path, uuids_file: Path, tissue: str = None):
    output_file_name = f"{tissue}" if tissue else "atac"
    uuids_df = pd.read_csv(uuids_file, sep="\t", dtype=str)
    uuids_list = uuids_df["uuid"].to_list()
    hbmids_list = uuids_df["hubmap_id"].to_list()
    directories = [data_directory / Path(uuid) for uuid in uuids_df["uuid"]]
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories if len(listdir(directory))>1]
    print("Annotating objects")
    cell_by_bin_adatas = [
        annotate_h5ads(file_pair[0], tissue, uuids_df)
        for file_pair in file_pairs
    ]

    cell_by_gene_adatas = [     
        annotate_h5ads(file_pair[1], tissue, uuids_df)
        for file_pair in file_pairs
    ]
    print("Concatenating objects")
    cbb_concat = anndata.concat(cell_by_bin_adatas, join="outer")
    cbg_concat = anndata.concat(cell_by_gene_adatas, join="outer")
    creation_time = str(datetime.now())
    data_product_uuid = str(uuid.uuid4())
    print("CBB shape: ", cbb_concat.obs.shape)
    total_cell_count = cbb_concat.obs.shape[0]
    mdata = make_mudata(cbb_concat, cbg_concat)
    mdata.obs = cbb_concat.obs
    mdata.obs = annotate_mudata(mdata, uuids_df)
    mdata.uns["creation_data_time"] = creation_time
    mdata.uns["datasets"] = hbmids_list
    mdata.uns["uuid"] = data_product_uuid
    mdata.write(f"{output_file_name}.h5mu")
    file_size = os.path.getsize(f"{output_file_name}.h5mu")
    create_json(tissue, data_product_uuid, creation_time, uuids_list, hbmids_list, total_cell_count, file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("data_directory", type=Path)
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory, args.uuids_file, args.tissue)