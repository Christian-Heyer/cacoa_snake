library(SoupX)

library(Matrix)
library(biomaRt)
library(Seurat)
library(scater)
library(scDblFinder)
library(BiocParallel)
library(biomaRt)
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 93)

if(exists("snakemake")){
  # Input
  mtx_file <- snakemake@input[["mtx_file"]]
  features_file <- snakemake@input[["features_file"]]
  adata_genes <- snakemake@input[["adata_genes"]]
  adata_cells <- snakemake@input[["adata_cells"]]
  adata_counts <- snakemake@input[["adata_counts"]]
  soupx_groups <- snakemake@input[["soupx_groups"]]
  ## Output
  corrected_counts_path <- snakemake@output[["corrected_counts"]]
  doublet_data_path <- snakemake@output[["doublet_data"]]
}else {
    mtx_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/GSM4472505_Basal_matrix.mtx.gz"
    features_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/GSM4472505_Basal_features.tsv.gz"
    adata_genes <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/csvs/var.csv"
    adata_cells <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/csvs/obs.csv"
    adata_counts <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/anndata/csvs/X.csv"
    soupx_groups <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/basal/qc/soupx_groups.csv"
    ## output
    corrected_counts_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/X_processed.csv"
    doublet_data_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/doublet_data.csv"
  }

data_tod <- Matrix::readMM(mtx_file)

old_genes <- vroom::vroom(features_file,
                          delim = "\t",col_names =  "ENSEMBL")
gene_bois <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),
      values = old_genes |> dplyr::pull(1), mart= ensembl)                            


rownames(data_tod) <- gene_bois |> dplyr::pull(external_gene_name) |> make.unique(sep = "-")

genes <- vroom::vroom(adata_genes)

cells <- vroom::vroom(adata_cells)

counts <- vroom::vroom(adata_counts, col_names = genes |> dplyr::pull(1)) 

# Match gnees and counts
counts <- as(t(counts), "sparseMatrix")
rownames(counts) <- genes |> dplyr::pull(gene_name)
colnames(counts) <- cells$cell_id

soupx_groups <- read.csv(soupx_groups)

# Dedpuplicate gene names 
rownames(data_tod)[duplicated(rownames(data_tod))] <- paste0(rownames(data_tod)[duplicated(rownames(data_tod))], "-1")


## Run SOUPX
sc = SoupChannel(data_tod[rownames(counts),], counts, calcSoupProfile = FALSE)

# Add extra meta data to the SoupChannel object
soupProf = data.frame(row.names = rownames(counts), est = rowSums(counts)/sum(counts), counts = rowSums(counts))
sc = setSoupProfile(sc, soupProf)
# Set cluster information in SoupChannel
sc = setClusters(sc, soupx_groups$soupx_groups)

sc  = autoEstCont(sc, doPlot=FALSE)

out = adjustCounts(sc, roundToInt = TRUE)

sce <- SingleCellExperiment(
  list(counts=out))

filter_index <- (out > 0 ) |> rowSums() > 20

sce <- sce[filter_index,]


genes <- genes[genes$gene_name %in% rownames(counts),]
data_mat = out

## RUN scDblFinder
sce  <- scDblFinder(sce = sce)

doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class

doublet_data <- tibble::tibble(cell_id = colnames(sce),doublet_score = doublet_score, doublet_class = doublet_class)
# Export Data 
readr::write_csv(doublet_data, file = doublet_data_path)

readr::write_csv(tibble::as_tibble(out, rownames = "genes"), file = corrected_counts_path)

