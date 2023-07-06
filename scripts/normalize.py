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
    adata_out_path = snakemake.output["adata_out_path"]
    scran_data_mat = snakemake.output["scran_data_mat"]
    cluster_groups = snakemake.output["cluster_groups"]
else:
    adata_path = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/adata.h5ad"
    out_mat ="/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/corrected_counts.csv"
    doublet_data_path = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/doublet_data.csv"

adata = sc.read_h5ad(adata_path)

out = pd.read_csv( out_mat, index_col="genes")
doublet_data = pd.read_csv(doublet_data_path, index_col="cell_id")

adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
adata.X = adata.layers["soupX_counts"]

sc.pp.filter_genes(adata, min_cells=20)

adata.obs["scDblFinder_score"] = doublet_data["doublet_score"]
adata.obs["scDblFinder_class"] = doublet_data["doublet_class"]

scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

#Scran steps
# Preliminary clustering for differentiated normalisation
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups")

data_mat = adata_pp.X.T
# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
adata.write_h5ad(adata_out_path)

np.savetxt(scran_data_mat, data_mat, delimiter = ",")
adata_pp.obs["groups"].to_csv(cluster_groups)
