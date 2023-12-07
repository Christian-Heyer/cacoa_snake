import anndata

adata = anndata.read(snakemake.input.h5ad)
adata.write_csvs(snakemake.params.out_dir, skip_data = False)
 
# write adata to disk
adata.write_h5ad(snakemake.output.h5ad)