import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


adata_files = list(map(sc.read_h5ad, snakemake.input["adata_files"]))
samples = snakemake.params["samples"]
cols = 3
tot = len(adata_files)

rows = tot // cols

if tot % cols !=0:
    rows+= 1

pos = range(1, tot+1)

for count, gene in enumerate(snakemake.params["genes"]):
    fig = plt.figure(figsize=(10,6))
    fig.suptitle(f"{gene}")
    for samp_index, adata in enumerate(adata_files):
        ax = fig.add_subplot(rows, cols, pos[samp_index]) 
        sc.pl.umap(
            adata,
	        color=[gene],     
	        vmin=0,
	        vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
	        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
	        frameon=False,
	        cmap="Reds",
            title=samples[samp_index],
            show=False,
            ax=ax
        )
    plt.savefig(snakemake.output[count], bbox_inches="tight")
