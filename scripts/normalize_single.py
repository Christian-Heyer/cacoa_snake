import scanpy as sc
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, issparse


if "snakemake" in globals():
    adata_path = snakemake.input["adata_h5"]
    adata_out_path = snakemake.output["adata_out_path"]
    scran_data_mat = snakemake.output["scran_data_mat"]
    cluster_groups = snakemake.output["cluster_groups"]
else:
    adata_path = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/adata.h5ad"
    adata_out_path = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/adata_normalized.h5ad"
    scran_data_mat = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/scran_data_mat.csv"
    cluster_groups = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/cluster_groups.csv"

# read adata into adata
adata = sc.read_h5ad(adata_path)

sc.pp.filter_genes(adata, min_cells=20)

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

#save sparse matrix data_mat to path in scran_data_mat


np.savetxt(scran_data_mat, data_mat.toarray(), delimiter=",")
adata_pp.obs["groups"].to_csv(cluster_groups)
