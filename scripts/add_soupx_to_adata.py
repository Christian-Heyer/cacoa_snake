import scanpy as sc

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from scipy.sparse import csr_matrix, issparse

if 'snakemake' in globals():
    adata_path = snakemake.input["adata_h5"]
    out_mat = snakemake.input["corrected_counts"]
    doublet_data_path = snakemake.input["doublet_data"]
    adata_out_path = snakemake.output["adata_h5"]

else:
    adata_path = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/adata.h5ad"
    out_mat ="/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/csvs/corrected_counts.csv"
    doublet_data_path = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/csvs/doublet_data.csv"

adata = sc.read_h5ad(adata_path)

out = pd.read_csv(out_mat, index_col="genes")
doublet_data = pd.read_csv(doublet_data_path, index_col="cell_id")

# Convert out.T to a sparse matrix
out_sparse = csr_matrix(out.values.T)

adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out_sparse  # Assign the sparse matrix to adata.layers["soupX_counts"]
adata.X = adata.layers["soupX_counts"]

adata.obs["scDblFinder_score"] = doublet_data["doublet_score"]
adata.obs["scDblFinder_class"] = doublet_data["doublet_class"]

# wrote adata to file
adata.write_h5ad(adata_out_path)































