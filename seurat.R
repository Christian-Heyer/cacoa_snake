library(Seurat)
library(SeuratObject)
library(reticulate)
options(future.globals.maxSize = 8000 * 1024^2)
library(anndataR)
library(SeuratDisk)

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
  use_anndataR <- snakemake@params[["use_anndataR"]]
  use_seuratdisk <- snakemake@params[["use_seuratdisk"]]
  umap_flag <- snakemake@params[["run_umap"]]
  default_layer <- snakemake@config[["default_layer"]]
  sample_cells <- snakemake@params[["sample_cells"]]
  source("common.R")
  library(glue)
} else {
  f_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/"
  dataset_name <- "TabularMuris"
  input_h5ad <- glue::glue("{f_path}/{dataset_name}/{dataset_name}_update.h5ad")
  input_file_counts <- file.path(f_path, glue::glue("{dataset_name}/{dataset_name}_raw/X.csv"))
  input_file_metadata <- file.path(f_path, glue::glue("{dataset_name}/{dataset_name}_raw/obs.csv"))
  input_file_genes <- file.path(f_path, glue::glue("{dataset_name}/{dataset_name}_raw/var.csv"))
  seurat_temp <- glue::glue("{f_path}/{dataset_name}/{dataset_name}_update.h5seurat")
  f_type <- "h5ad"
  organism <- "hs"
  p_name <- "testproject"
  use_seuratdisk <- T
  default_layer <- "counts"
}


if (f_type == "h5ad") {
  if (use_seuratdisk) {
  #  seurat_obj <- anndataR::read_h5ad(input_h5ad, to = "Seurat")
     seurat_temp <- seurat_temp
     library(SeuratDisk)
     Convert(
      source = input_h5ad, dest = seurat_temp,
      overwrite = TRUE
     )
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
    col_data <- as.data.frame(col_data, row.names = col_data$index)
    rownames(col_data) <- col_data$index
    col_data <- col_data[, -c(1)]
    seurat_obj <- CreateSeuratObject(
      counts = counts, project = p_name,
      meta.data = col_data
    )
  }
} else if (f_type == "seuratRDS") {
  seurat_obj <- readRDS(input_file)
  saveRDS(object = seurat_obj, file = seurat_temp)
}
seurat_obj$nCount_RNA <- colSums(x = seurat_obj, slot = "counts") # nCount_RNA
# Assuming `seurat_obj` is your Seurat object and `sample_cells` is the input parameter
if (!is.null(sample_cells) && is.numeric(sample_cells) && sample_cells > 0) {
  # Ensure sample_cells does not exceed the total number of cells
  total_cells <- ncol(seurat_obj)
  if (sample_cells > total_cells) {
    stop("sample_cells exceeds the total number of cells in the Seurat object.")
  }

  # Randomly sample cell names
  sampled_cells <- sample(colnames(seurat_obj), size = sample_cells, replace = FALSE)

  # Subset the Seurat object
  seurat_obj <- subset(seurat_obj, cells = sampled_cells)

  # Return or assign the downsampled object
} else {
  write("No sampling", file = stderr())
}
# nFeature = colSums(x = GetAssayData(object = seurat_obj, slot = "counts") > 0)
if (organism == "mm") {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt-")
} else {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}
DefaultLayer(seurat_obj@assays$RNA) <- default_layer
# seurat_obj <- processSeurat(seurat_obj, is_raw = is_raw, run_umap = umap_flag)

saveRDS(seurat_obj, output_file)
