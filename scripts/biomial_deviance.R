library(readr)
library(scry)
library(Matrix)

if (exists("snakemake")) {
	data_mat <- snakemake@input[["corrected_counts"]]
	gene_dat <- snakemake@input[["var_filtered"]]
	out_file <- snakemake@output[["deviance_data"]]
} else {
  data_mat <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/corrected_counts_filtered.csv"
  gene_dat <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/lung_lps/anndata/csvs/var_filtered.csv"
}

count_mat <- vroom::vroom(data_mat)
gene_mat <- vroom::vroom(gene_dat)
gene_names <- gene_mat |> dplyr::pull(gene_name)
count_mat <- count_mat[, -1]
count_mat <- as(as.matrix(count_mat), "sparseMatrix")
rownames(count_mat) <- gene_names


deviance_data = devianceFeatureSelection(count_mat)

readr::write_csv(tibble::tibble(genes = gene_names, binomial_deviance = deviance_data),
				 file = out_file)