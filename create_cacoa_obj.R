library(Seurat)
library(magrittr)
library(dplyr)
library(future)

if(!require("cacoa")) {
  install.packages("coda.base", repos =  "https://cloud.r-project.org")
  #devtools::install_github("kharchenkolab/sccore", ref = "dev")
  devtools::install_github("kharchenkolab/sccore", ref = "32a8f52", upgrade = "never")
  devtools::install_github("kharchenkolab/cacoa", ref = "be2a38d", upgrade = "never")
  
}

library(cacoa)
#Sys.setenv(https_proxy='http://www-int.dkfz-heidelberg.de:80')
#Sys.setenv(http_proxy='http://www-int.dkfz-heidelberg.de:80')
if(exists("snakemake")) {
    seurat_path <- snakemake@input[["seurat_path"]]
    output_p <- snakemake@output[["cacoa_obj"]]
    umap_path <- snakemake@output[['umap_coords']]
    permute <- snakemake@params[["permute"]]
    plan("multicore", workers = snakemake@threads)

    cacoa_opts <- snakemake@config[["cacoa_opts"]]
} else {
    base_fp = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/TabularMuris_counts"
    seurat_path = file.path(base_fp, "seurat_obj.RDS.gz")
    output_p = file.path(base_fp, "cao_obj.rds.gz")
    config <- yaml::read_yaml("/desktop-home/heyer/projects/Vascular_Aging/RNAseq/scRNAseq_scripts/configs/tabularmuris_counts.yaml")
    cacoa_opts <- config$cacoa_opts
    permute <- F
}

create_cao_from_seurat <- function(s_path = seurat_path, do_permute = permute) {
    options(warn=-1)
    # Run standard Seurat processing
    seurat_obj <- readRDS(seurat_path)
    #seurat_obj <- SCTransform(seurat_obj,verbose = F)

    #seurat_obj <- RunPCA(seurat_obj, 
     #                    features = VariableFeatures(object = seurat_obj))
    #seurat_obj <- FindNeighbors(seurat_obj, features = VariableFeatures(object = seurat_obj))
    #seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
    #options(warn=0)
    seurat_obj<- SetIdent(seurat_obj, value = cacoa_opts[["Ident"]])
     
    if (cacoa_opts[["regroup_age"]]) {
        if(snakemake@config[["organism"]] == "mm") {
            seurat_obj$binary_age <- ifelse(seurat_obj$age %in% c("1m", "3m"),
                                    "young", 
                                    "aged")
        } else {
            seurat_obj$binary_age <- ifelse(seurat_obj$Age < 47,
                                    "young", 
                                    "aged")
            seurat_obj$Subject_Identity<- paste0(seurat_obj$Subject_Identity, "_", seurat_obj$Age)
        }
    } else {
        seurat_obj$binary_age <- seurat_obj$age
        seurat_obj$mouse.id <- paste0(seurat_obj$age, "_", seurat_obj$animal)
    } 
    # Organize Metadata
    meta_meta <- seurat_obj@meta.data %>% dplyr::select(!!cacoa_opts$age_field, 
                                                        !!cacoa_opts$id_field, 
                                                        binary_age)

    print(cacoa_opts)
    print(meta_meta)
    s_per_cell <- setNames(dplyr::pull(meta_meta, !!cacoa_opts$id_field), nm = rownames(meta_meta))
    mouse_meta <- meta_meta %>% unique()
    if (do_permute) {
        mouse_vec <- setNames(sample(c("young", "aged"), 
                         size = nrow(mouse_meta), replace = T), nm = dplyr::pull(mouse_meta, !!cacoa_opts$id_field))
    } else {
        mouse_vec <- setNames(mouse_meta$binary_age, nm= dplyr::pull(mouse_meta, !!cacoa_opts$id_field))
    }
    # Create object and set coloring
    print(unique(mouse_vec))
    yeetme <- Cacoa$new(seurat_obj, target.level = "aged", 
                        ref.level ="young", 
                        sample.groups= mouse_vec, 
                        sample.per.cell =s_per_cell, graph.name = "UMAP")
    yeetme$cell.groups.palette<- levels(yeetme$cell.groups) %>%  
    {setNames(sample(brewerPalette("Paired")(length(.))), .)}
    
    yeetme
}

cao_obj <- create_cao_from_seurat(seurat_path, permute)

seurat_obj <- cao_obj$data.object
#get Umap coords
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

write.csv(umap_coords, file = umap_path)
saveRDS(cao_obj, file = output_p)
