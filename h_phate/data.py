import pandas as pd
import scprep
import os
import numpy as np
import subprocess

def load_eb():
    download_path = os.path.expanduser("~")
    if not os.path.isdir(os.path.join(download_path, "scRNAseq", "T0_1A")):
        # need to download the data
        scprep.io.download.download_and_extract_zip(
            "https://md-datasets-public-files-prod.s3.eu-west-1.amazonaws.com/"
            "5739738f-d4dd-49f7-b8d1-5841abdbeb1e",
            download_path)
    sparse=True
    T1 = scprep.io.load_10X(os.path.join(download_path, "scRNAseq", "T0_1A"), sparse=sparse, gene_labels='both')
    T2 = scprep.io.load_10X(os.path.join(download_path, "scRNAseq", "T2_3B"), sparse=sparse, gene_labels='both')
    T3 = scprep.io.load_10X(os.path.join(download_path, "scRNAseq", "T4_5C"), sparse=sparse, gene_labels='both')
    T4 = scprep.io.load_10X(os.path.join(download_path, "scRNAseq", "T6_7D"), sparse=sparse, gene_labels='both')
    T5 = scprep.io.load_10X(os.path.join(download_path, "scRNAseq", "T8_9E"), sparse=sparse, gene_labels='both')
    filtered_batches = []
    for batch in [T1, T2, T3, T4, T5]:
        batch = scprep.filter.filter_library_size(batch, percentile=20, keep_cells='above')
        batch = scprep.filter.filter_library_size(batch, percentile=75, keep_cells='below')
        filtered_batches.append(batch)
    del T1, T2, T3, T4, T5 # removes objects from memory
    data, sample_labels = scprep.utils.combine_batches(
    filtered_batches, 
    ["Day 00-03", "Day 06-09", "Day 12-15", "Day 18-21", "Day 24-27"],
    append_to_cell_names=True
)
    del filtered_batches # removes objects from memory
    data = scprep.filter.filter_rare_genes(data, min_cells=10)
    data, sample_labels = scprep.filter.filter_gene_set_expression(
    data, sample_labels, starts_with="MT-", library_size_normalize=True,
    percentile=90, keep_cells='below')
    return normalize(data)

def load_pbmc():
    scprep.io.download.download_url("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", "pbmc4k_filtered_gene_bc_matrices.tar.gz")
    subprocess.run(['tar', 'xzf', "pbmc4k_filtered_gene_bc_matrices.tar.gz"])
    data = scprep.io.load_10X("filtered_gene_bc_matrices/GRCh38", gene_labels='both')
    data = scprep.filter.remove_rare_genes(data, min_cells=5)
    data = scprep.filter.remove_empty_cells(data)
    return normalize(data)

def normalize(data):
    data = scprep.normalize.library_size_normalize(data)
    data = scprep.transform.sqrt(data)
    data = scprep.transform.arcsinh(data)
    return data