import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

if 'snakemake' in globals():
    adata_path = snakemake.input["adata_h5"]
    deviance_path = snakemake.input["deviance_data"]
    output_adata = snakemake.output["adata_out_path"]
    output_paths = snakemake.output

adata = sc.read_h5ad(adata_path)
deviance_data = pd.read_csv(deviance_path)

deviance_data = deviance_data.set_index("genes")


idx = deviance_data["binomial_deviance"].argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = deviance_data["binomial_deviance"]


def plot_pca_umap(adata_obj, outpath):
    adata_obj.X = adata.layers["log1p_norm"]

    adata_obj.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata_obj, svd_solver="arpack", use_highly_variable=True)

    with plt.rc_context():
        sc.pl.pca_scatter(adata_obj, color="total_counts")
        plt.savefig(outpath["pca"], bbox_inches="tight")

    sc.pp.neighbors(adata_obj)
    sc.tl.umap(adata_obj)
    sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

    with plt.rc_context():
        sc.pl.umap(
            adata_obj,
            color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
            ncols = 2
        )
        plt.savefig(outpath["umap"], bbox_inches="tight")
    with plt.rc_context():
        sc.pl.umap(
            adata,
            color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
            legend_loc="on data",
        )
        plt.savefig(outpath["umap_cluster"], bbox_inches="tight")

    return adata_obj


adata = plot_pca_umap(adata, output_paths)

adata.write_h5ad(output_adata)
# export umap coords to csv
np.savetxt(output_paths["umap_coords"], adata.obsm["X_umap"], delimiter=",")