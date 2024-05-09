library(SeuratObject)

library(Seurat)
# Define input and output files using snakemake settings
if (exists("snakemake")) {
    input_file <- snakemake@input[["adata_processed"]]
    umap_coords <- snakemake@input[["umap_coords"]]
    output_file <- snakemake@output$seurat_obj
    use_anndataR <- snakemake@params[['use_anndataR']]
} else {
  input_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/merged/anndata/processed_adata.h5ad"
  umap_coords <-"/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/merged/anndata/csvs/umap_coords.csv" 
  output_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/merged/seurat/seurat_obj_processed.rds.gz"
}
if (!require("anndataR")) {
  devtools::install_github("scverse/anndataR")
}

if (use_anndataR) {
  seurat_obj <- anndataR::read_h5ad(input_file, to = "Seurat")
} else { 
  library(SeuratDisk)
  Convert(source = input_file, dest = output_file, 
          overwrite = TRUE)
  seurat_obj <- readRDS(output_file)
}

# read umap coords
UMAP_coordinates <- read.table(umap_coords, sep = ",", header = T, stringsAsFactors = F, row.names=1)

#replace _ in rownames with -
rownames(UMAP_coordinates) <- gsub("_", "-", rownames(UMAP_coordinates))

# add umap coords as reduction to seurat object
# conver to matrix
UMAP_coordinates_mat <- as(UMAP_coordinates, "matrix")
colnames(UMAP_coordinates_mat) <- c("UMAP-1", "UMAP-2")
# Create DimReducObject and add to object
seurat_obj[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, 
                                              key = "UMAP_", 
                                              global = T, assay = "RNA")

# save seurat object
saveRDS(seurat_obj, output_file)