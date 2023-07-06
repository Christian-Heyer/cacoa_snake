library(scran)
library(BiocParallel)
if(exists("snakemake")) {
	# Input data
	data_mat <- snakemake@input[["scran_data_mat"]]
	cluster_groups <- snakemake@input[["cluster_groups"]]
	size_factor_file <- snakemake@output[["size_factors"]]
	bp_param <- MulticoreParam(workers = snakemake@threads)
} else {
  data_mat <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/norm_data_mat.csv"
  cluster_groups <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/cluster_groups.csv"
  
}
cluster_groups <- vroom::vroom(cluster_groups)
data_mat <- vroom::vroom(data_mat, col_names = cluster_groups$cell_id)

cluster_vec <- setNames( cluster_groups$groups, nm = cluster_groups$cell_id)

size_factors = sizeFactors(
    scran::computeSumFactors(
        SingleCellExperiment(
            list(counts=data_mat)), 
            clusters = cluster_vec,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
    )
)

readr::write_csv(tibble::tibble(cells = cluster_groups$cell_id, 
                        size_factors = size_factors), size_factor_file )
