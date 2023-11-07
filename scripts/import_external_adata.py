import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
from matplotlib import pyplot as plt
from process_anndata import prepare_soupx

def read_and_process_data(adata_input,
                          out_paths: dict):
	col_name=snakemake.wildcards.sample
    adata = sc.read_h5ad(adata_input)
	adata.obs["sample"] = col_name
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
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    adata.obs.mt_outlier.value_counts()

    with plt.rc_context():
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
        plt.savefig(out_paths["p3_scatter_after"], bbox_inches="tight") 
    
    return(adata)

adata_obj = read_and_process_data(adata_input=snakemake.input["input_adata"],
                                     out_paths=snakemake.output)

adata_obj.write(snakemake.output["adata_h5"])
adata_obj.write_csvs(snakemake.params["adata_csvs"], skip_data = False)
prepare_soupx(adata_obj, snakemake.output)