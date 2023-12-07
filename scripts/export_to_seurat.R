
library(Seurat)
# Define input and output files using snakemake settings
if (exists("snakemake")) {
    input_file <- snakemake@input[["adata_processed"]]
    output_file <- snakemake@output$seurat_obj
    use_anndataR <- snakemake@params[['use_anndataR']]
} else {
  input_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/merged/anndata/processed_adata.h5ad"
  output_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/merged/seurat/seurat_obj_processed.rds.gz"
}
if (!require("anndataR")) {
  devtools::install_github("scverse/anndataR")
}

if (use_anndataR) {
  seurat_obj <- anndataR::read_h5ad(input_file, to = "Seurat")
  saveRDS(seurat_obj, file = output_file)
} else { 
  library(SeuratDisk)
  Convert(source = input_file, dest = output_file, 
          overwrite = TRUE)
}
