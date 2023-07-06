import scanpy as sc
import os

if "snakemake" in globals():
    input_h5ads = snakemake.input["adata_files"]
    output_path = snakemake.output["big_adata"]

adata = list(map(sc.read_h5ad, input_h5ads))

big_adata = sc.concat(adata, join="outer", merge="same")


big_adata.write_h5ad(output_path)