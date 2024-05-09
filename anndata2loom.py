import anndata

adata = anndata.read_h5ad(snakemake.input.h5ad)

# if snakemake.config.filter has an entry filter obs by it
if snakemake.config["filter"]:
# Filter is a dict of key = column name in obs and value is the value to keep
    for key, value in snakemake.config["filter"].items():
        adata = adata[adata.obs[key] == value, :]

adata.write_csvs(snakemake.params.out_dir, skip_data = False)
 
# write adata to disk
adata.write_h5ad(snakemake.output.h5ad)