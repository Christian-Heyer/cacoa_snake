import scanpy as sc
import numpy as np
import pandas as pd

if "snakemake" in globals():
    input_h5ads = snakemake.input["adata_files"]
    output_path = snakemake.output["big_adata"]
    umap_output = snakemake.output["umap_coords"]
    batch_key = snakemake.params["batch_key"]
    n_top_genes = snakemake.params["n_top_genes"]

adata = list(map(sc.read_h5ad, input_h5ads))

big_adata = sc.concat(adata, join="outer", merge="same")

sc.pp.filter_genes(big_adata, min_cells=5)

big_adata.X = big_adata.layers["counts"].copy()
sc.pp.normalize_total(big_adata)
sc.pp.log1p(big_adata)
big_adata.layers["log1p_norm"] = big_adata.X.copy()

sc.pp.highly_variable_genes(
    big_adata, n_top_genes=n_top_genes, flavor="cell_ranger", batch_key=batch_key
)
sc.tl.pca(big_adata)
sc.pp.neighbors(big_adata, n_pcs=15)
sc.tl.umap(big_adata)

# subset big_adata by variabel genes

big_adata.write_h5ad(output_path)

# Get UMAP coordinates and add cell_id from adata_big as an additional column to the umap coords
umap_coords = big_adata.obsm["X_umap"]
# add cell_id column
cell_id = big_adata.obs.index
# concatenate cell_id to umap_coords

umap_coords = pd.concat([pd.Series(cell_id), pd.DataFrame(umap_coords)], axis=1)
# big_adata.obs["cell_id"] = big_adata.obs.index


# export umap coords to csv
pd.DataFrame(umap_coords).to_csv(umap_output, index=False, header=True)
