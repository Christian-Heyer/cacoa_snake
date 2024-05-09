library(Seurat)
library(SeuratDisk)
library(reticulate)
if (exists("snakemake")) {
  input_file_counts <- snakemake@input[["counts"]]
  input_file_metadata <- snakemake@input[["cell_metadata"]]
  input_file_genes <- snakemake@input[["gene_metadata"]]
  input_h5ad <- snakemake@input[["h5ad"]]
  seurat_temp <- snakemake@output[["seurat_temp"]]
  output_file <- snakemake@output[["seurat_path"]]
  f_type <- snakemake@params[["f_type"]]
  organism <- snakemake@config[["organism"]]
  p_name <- snakemake@config[["DATASET"]]
  is_raw <- snakemake@config[["is_raw"]]
  use_anndataR <- snakemake@params[['use_anndataR']]
  use_seuratdisk <- snakemake@params[['use_seuratdisk']]
  source("common.R")
library(glue)

} else {
  f_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/"
  dataset_name <- "TabularMuris_counts"
  input_h5ad <- glue::glue("{f_path}/{dataset_name}/{dataset_name}_update.h5ad")
  input_file_counts <- file.path(f_path, glue::glue("{dataset_name}/{dataset_name}_raw/X.csv"))
  input_file_metadata <- file.path(f_path, glue::glue("{dataset_name}/{dataset_name}_raw/obs.csv"))
  input_file_genes <- file.path(f_path, glue::glue("{dataset_name}/{dataset_name}_raw/var.csv"))
  seurat_temp <- glue::glue("{f_path}/{dataset_name}/{dataset_name}_update.h5seurat")
  f_type <- "h5ad"
  organism <- "mm"
  p_name <- "testproject"
}


if (f_type == "h5ad" ) {
  if (use_seuratdisk) {
    #seurat_obj <- anndataR::read_h5ad(input_h5ad, to = "Seurat")
    seurat_temp <- seurat_temp
    Convert(source = input_h5ad, dest = seurat_temp, 
            overwrite = T)
    seurat_obj <- LoadH5Seurat(seurat_temp, meta.data = FALSE)
    # get metadat from csvs
    col_data <- vroom::vroom(input_file_metadata)
    col_data <- as.data.frame(col_data, row.names =col_data$index ) 
    rownames(col_data) <- col_data$index
    col_data <- col_data[,-c(1)]
    # Add metadata to Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = col_data)
  } else {
    col_data <- vroom::vroom(input_file_metadata)
    row_data <- vroom::vroom(input_file_genes)
    counts <- vroom::vroom(input_file_counts, col_names = row_data$index)
    counts <- t(counts)
    colnames(counts) <- col_data$index
    col_data <- as.data.frame(col_data, row.names =col_data$index ) 
    rownames(col_data) <- col_data$index
    col_data <- col_data[,-c(1)]
    seurat_obj <- CreateSeuratObject(counts = counts, project = p_name, 
                                    meta.data = col_data)
  }
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


seurat_obj <- processSeurat(seurat_obj, is_raw = is_raw)

saveRDS(seurat_obj, output_file)