library(GEOquery)
library(magrittr)
library(Seurat)
library(stringr)
library(dplyr)
if (exists("snakemake")) {
    source("common.R")
    counts_file <- snakemake@input[["counts"]]
    metadata_file <- snakemake@input[["metadata"]]
    genes_file <- snakemake@input[["genes"]]
    patient_data <- snakemake@params[["patient_data"]]
    GSE_id <- snakemake@params[["geo_data"]]
    GEO_data <- snakemake@config[["download_link"]][["GEO"]]
    seurat_path <- snakemake@output[["seurat_path"]]
}
#DATA <- GEOquery::getGEOSuppFiles(GSE_id)
#f_names <-  rownames(DATA)

#metadata <- purrr::map(f_names[1:3], readr::read_tsv)
  #<-GSE136831(file_path)
meta_data <- readr::read_tsv(metadata_file)
meta_data %>% as.data.frame() -> meta_frame
rownames(meta_frame) <- meta_frame$CellBarcode_Identity
#meta_frame <- meta_frame[,-c(1)]
count_mtx <- Matrix::readMM(counts_file)
gene_data <- readr::read_tsv(genes_file)
rownames(count_mtx) <- gene_data$HGNC_EnsemblAlt_GeneID

colnames(count_mtx) <- rownames(meta_frame)


age_test <- readr::read_delim(patient_data, delim = " ")
meta_frame_comb <- dplyr::full_join(meta_frame, age_test, by =  c("Subject_Identity" = "Subject_ID")) %>% as.data.frame()
rownames(meta_frame_comb) <- meta_frame_comb$CellBarcode_Identity
meta_frame_comb$Manuscript_Identity %>% str_replace(pattern = "^VE_.*", replacement = "Vascular Endothelial") -> meta_frame_comb$Manuscript_Identity
#count_mtx <- count_mtx[,rownames(meta_frame_comb) %in% colnames(count_mtx)]

#row.names(meta_frame) <- meta_frame_comb$CellBarcode_Identity
seurat_obj <- CreateSeuratObject(count_mtx, project = "HumanLung", assay = "RNA", meta.data = meta_frame_comb)
rm(count_mtx)
gc()
seurat_obj <- seurat_obj[, seurat_obj$Disease_Identity == "Control" & seurat_obj$Subject_Identity %in% age_test$Subject_ID]

#seurat_obj@meta.data <- 
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 15)
seurat_obj <- processSeurat(seurat_obj)

saveRDS(seurat_obj,
        file = seurat_path)
