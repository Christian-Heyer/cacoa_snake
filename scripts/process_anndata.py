import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
from matplotlib import pyplot as plt


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def read_and_process_10x(mtx_file: str,
                         barcodes_file: str, 
                         features_file: str,
                         out_paths: dict):
    col_name=snakemake.wildcards.sample
    adata = sc.read_mtx(mtx_file)
    adata_bc=pd.read_csv(barcodes_file, header=None)
    adata_features=pd.read_csv(features_file, header=None)
    # Hardcocded ensembl 93 as host -> make flexible
    annot = sc.queries.biomart_annotations(
        "mmusculus",
        ["ensembl_gene_id", "external_gene_name"]
#        host="jul2018.archive.ensembl.org"
    ).set_index("ensembl_gene_id")

    adata = adata.T
    adata.obs['cell_id']= adata_bc[0].tolist()
    adata.obs["sample"] = col_name
    adata.obs = adata.obs.set_index("cell_id")
    adata.var['gene_id']= adata_features[0].tolist()
    adata.var = adata.var.set_index("gene_id")
    # Hardcoded ensembl 93 <- Set as parameter
    adata.var['gene_name'] = annot["external_gene_name"]
    adata.var["gene_name"] = adata.var["gene_name"].astype(str)
    adata.var =adata.var.set_index("gene_name")

    adata.var_names_make_unique()
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^Hb[^(P)]"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    
    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False,)
    p1.set_titles(col_name)
    p1.figure.savefig(out_paths["p1_displot"])
# sc.pl.violin(adata, 'total_counts')

    with plt.rc_context(): 
        sc.pl.violin(adata, "pct_counts_mt", title=col_name)
        plt.savefig(out_paths["p2_mt_counts"], bbox_inches="tight")

    with plt.rc_context():
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", 
                       color="pct_counts_mt",title=col_name)
        plt.savefig(out_paths["p3_scatter"], bbox_inches="tight")
    
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs.outlier.value_counts()
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 4) | (
        adata.obs["pct_counts_mt"] > 15
    )
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    adata.obs.mt_outlier.value_counts()

    with plt.rc_context():
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
        plt.savefig(out_paths["p3_scatter_after"], bbox_inches="tight")

    return(adata)


def prepare_soupx(adata, out_paths):
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    soupx_groups = adata_pp.obs["soupx_groups"]

    del adata_pp
    soupx_groups.to_csv(out_paths["soupx_groups"])
    



adata_obj = read_and_process_10x(
                     mtx_file=snakemake.input["mtx_file"],
                     barcodes_file=snakemake.input["barcodes_file"],
                     features_file=snakemake.input["features_file"],
                     out_paths=snakemake.output,
)

adata_obj.write(snakemake.output["adata_h5"])
adata_obj.write_csvs(snakemake.params["adata_csvs"], skip_data = False)
prepare_soupx(adata_obj, snakemake.output)
