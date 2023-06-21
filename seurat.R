library(Seurat)
library(reticulate)
if (exists("snakemake")) {
  input_file_counts <- snakemake@input[["counts"]]
  input_file_metadata <- snakemake@input[["cell_metadata"]]
  input_file_genes <- snakemake@input[["gene_metadata"]]
  seurat_temp <- snakemake@output[["seurat_temp"]]
  output_file <- snakemake@output[["seurat_path"]]
  f_type <- snakemake@params[["f_type"]]
  organism <- snakemake@config[["organism"]]
  p_name <- snakemake@config[["DATASET"]]
  source("common.R")
} else {
  f_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/"
  
	input_file_counts <- file.path(f_path, "TabularMuris_nonmyeloid_brain/TabularMuris_nonmyeloid_brain_raw/X.csv")
	input_file_metadata <- file.path(f_path, "TabularMuris_nonmyeloid_brain/TabularMuris_nonmyeloid_brain_raw/obs.csv")
	input_file_genes <- file.path(f_path, "TabularMuris_nonmyeloid_brain/TabularMuris_nonmyeloid_brain_raw/var.csv")
	f_type <- "h5ad"
	seurat_temp <- "TabularMuris_nonmyeloid_brain_raw/yeet.rds"
	organism <- "mm"
	p_name <- "testproject"
}


if (f_type == "h5ad") {
  col_data <- vroom::vroom(input_file_metadata)
  row_data <- vroom::vroom(input_file_genes)
  counts <- vroom::vroom(input_file_counts, col_names = row_data$index)
  counts <- t(counts)
  colnames(counts) <- col_data$index
  col_data <- as.data.frame(col_data, row.names =col_data$index ) 
  rownames(col_data) <- col_data$index
  col_data <- col_data[,-c(1)]
  seurat_obj <- CreateSeuratObject(counts = counts, project = p_name, meta.data = col_data)
  #seurat_obj$nCount_RNA <- seurat_obj$nCounts_RNA
} else if (f_type == "seuratRDS") {
  seurat_obj <- readRDS(input_file)
  saveRDS(object = seurat_obj, file = seurat_temp)
}
seurat_obj$nCount_RNA <- colSums(x = seurat_obj, slot = "counts")  # nCount_RNA
#nFeature = colSums(x = GetAssayData(object = seurat_obj, slot = "counts") > 0) 
if (organism == "mm") {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt-")
} else {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}


seurat_obj <- processSeurat(seurat_obj)

saveRDS(seurat_obj, output_file)