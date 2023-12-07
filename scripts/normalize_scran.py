import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import numpy as np

size_factors = pd.read_csv(snakemake.input["size_factors"]).set_index("cells")
adata = sc.read_h5ad(snakemake.input["adata_out_path"])
adata.obs["size_factors"] = size_factors["size_factors"]
scran = adata.X / adata.obs["size_factors"].values[:, None]
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))

adata.write_h5ad(snakemake.output["adata_out_path"])
np.savetxt(snakemake.output["corrected_counts_filtered"], adata.X.T.toarray(),
           delimiter=",", header =", ".join(adata.obs.index))
adata.var.to_csv(snakemake.output["var_filtered"])