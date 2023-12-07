import scanpy as sc

if "snakemake" in globals():
    input_h5ads = snakemake.input["adata_files"]
    output_path = snakemake.output["big_adata"]
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

