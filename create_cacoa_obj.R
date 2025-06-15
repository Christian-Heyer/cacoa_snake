library(Seurat)
library(magrittr)
library(dplyr)
library(future)

if (!require("cacoa")) {
  install.packages("coda.base", repos = "https://cloud.r-project.org")
  devtools::install_github("kharchenkolab/sccore",
    ref = "32a8f52", upgrade = "never"
  )
  devtools::install_github("kharchenkolab/cacoa",
    ref = "be2a38d", upgrade = "never"
  )
}

library(cacoa)
if (exists("snakemake")) {
  seurat_path <- snakemake@input[["seurat_path"]]
  output_p <- snakemake@output[["cacoa_obj"]]
  umap_path <- snakemake@output[["umap_coords"]]
  permute <- snakemake@params[["permute"]]
  plan("multicore", workers = snakemake@threads)
  count_layer <- snakemake@config[["default_layer"]]
  cacoa_opts <- snakemake@config[["cacoa_opts"]]
  translate_ensembl <- snakemake@config[["translate_ensembl"]]
} else {
  base_fp <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/LungAtlas/"
  seurat_path <- file.path(base_fp, "seurat_obj.RDS.gz")
  output_p <- file.path(base_fp, "cao_obj.rds.gz")
  config <- yaml::read_yaml("/desktop-home/heyer/projects/Vascular_Aging/RNAseq/scRNAseq_scripts/configs/manlung.yaml")
  cacoa_opts <- config$cacoa_opts
  permute <- F
}
ConvertEnsemblToGeneNames <- function(seurat_obj, species = "human") {
  # Load required libraries
  if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
  library(biomaRt)
  library(Seurat)
  
  # Determine dataset based on species
  dataset <- switch(tolower(species),
                    "human" = "hsapiens_gene_ensembl",
                    "mouse" = "mmusculus_gene_ensembl",
                    stop("Unsupported species. Use 'human' or 'mouse'."))
  
  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = dataset)

  # Extract Ensembl IDs
  ensembl_ids <- rownames(seurat_obj)

  # Get mapping from Ensembl to gene symbols
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )

  # Clean gene_map: remove duplicates and empty gene symbols
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  gene_map <- gene_map[gene_map$hgnc_symbol != "", ]

  # Map gene names to Seurat object
  new_rownames <- gene_map$hgnc_symbol[match(rownames(seurat_obj), gene_map$ensembl_gene_id)]

  # Apply new rownames
  seurat_obj <- seurat_obj[!is.na(new_rownames), ]
  rownames(seurat_obj) <- new_rownames[!is.na(new_rownames)]

  return(seurat_obj)
}

create_cao_from_seurat <- function(s_path = seurat_path, do_permute = permute) {
  options(warn = -1)
  # Run standard Seurat processing
  seurat_obj <- readRDS(seurat_path)
  seurat_obj <- SetIdent(seurat_obj, value = cacoa_opts[["Ident"]])

  if (cacoa_opts[["regroup_age"]]) {
    if (snakemake@config[["organism"]] == "mm") {
      seurat_obj$binary_age <- ifelse(seurat_obj$age %in% c("1m", "3m"),
        "young",
        "aged"
      )
    } else {
      if (!(cacoa_opts$age_field %in% colnames(seurat_obj@meta.data))) {
        library(dplyr)
        library(stringr)

        seurat_obj@meta.data <- seurat_obj@meta.data %>%
          mutate("{cacoa_opts$age_field}" := str_split(development_stage, "-") %>%
            purrr::map_chr(1))
      }
      seurat_obj$binary_age <- ifelse(seurat_obj$Age < 47,
        "young",
        "aged"
      )
      seurat_obj@meta.data[,cacoa_opts$id_field]  <- paste0(seurat_obj@meta.data[,cacoa_opts$id_field] , "_", seurat_obj@meta.data[,cacoa_opts$age_field] )
    }
  } else {
    seurat_obj$binary_age <- seurat_obj$age
    seurat_obj$mouse.id <- paste0(seurat_obj$age, "_", seurat_obj$animal)
  }
  # Organize Metadata
  meta_meta <- seurat_obj@meta.data %>% dplyr::select(
    !!cacoa_opts$age_field,
    !!cacoa_opts$id_field,
    binary_age
  )
  if (translate_ensembl) {
  #  seurat_obj <- ConvertEnsemblToGeneNames(seurat_obj, "human")
    new_rownames <- seurat_obj@assays$RNA@meta.data$feature_name
    new_rownames <- make.unique(as.character(new_rownames))
    #seurat_obj <- seurat_obj[!is.duplicated(new_rownames), ]
    rownames(seurat_obj) <- new_rownames[!is.na(new_rownames)]
  }
  print(cacoa_opts)
  print(meta_meta)
  s_per_cell <- setNames(dplyr::pull(meta_meta, !!cacoa_opts$id_field), nm = rownames(meta_meta))
  mouse_meta <- meta_meta %>% unique()
  if (do_permute) {
    mouse_vec <- setNames(sample(c("young", "aged"),
      size = nrow(mouse_meta), replace = T
    ), nm = dplyr::pull(mouse_meta, !!cacoa_opts$id_field))
  } else {
    mouse_vec <- setNames(mouse_meta$binary_age, nm = dplyr::pull(mouse_meta, !!cacoa_opts$id_field))
  }
  # Create object and set coloring
  print(unique(mouse_vec))
  yeetme <- Cacoa$new(seurat_obj,
    target.level = "aged",
    ref.level = "young",
    sample.groups = mouse_vec,
    sample.per.cell = s_per_cell, #graph.name = "umap"
    data.slot = count_layer
  )
  yeetme$cell.groups.palette <- levels(yeetme$cell.groups) %>%
    {
      setNames(sample(brewerPalette("Paired")(length(.))), .)
    }

  yeetme
}

cao_obj <- create_cao_from_seurat(seurat_path, permute)

seurat_obj <- cao_obj$data.object
# get Umap coords

umap_coords <- Embeddings(seurat_obj, reduction = "umap")

write.csv(umap_coords, file = umap_path)
saveRDS(cao_obj, file = output_p)
