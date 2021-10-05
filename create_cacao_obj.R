library(cacoa)
library(Seurat)
library(magrittr)
library(future)
plan("multicore", workers = 4)
plan()

permute <- TRUE


#Sys.setenv(https_proxy='http://www-int.dkfz-heidelberg.de:80')
#Sys.setenv(http_proxy='http://www-int.dkfz-heidelberg.de:80')
if(exists("snakemake")) {
    seurat_path <- snakemake@input[["seurat_path"]]
    cacoa_obj <- snakemake@output[["cacoa_obj"]]
    permute <- snakemake@params[["permute"]]
} else {
    base_fp = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/TabularMuris/"
    seurat_path = file.path(base_fp, "senis_droplet_lung.rds")
    output_p = file.path(base_fp, "cacoa_obj.rds.gz")
}

create_cao_from_seurat <- function(s_path = seurat_path, do_permute = permute) {
    options(warn=-1)
    # Run standard Seurat processing
    seurat_obj <- readRDS(seurat_path)
    seurat_obj <- SCTransform(seurat_obj,verbose = F)

    seurat_obj <- RunPCA(seurat_obj, 
                         features = VariableFeatures(object = seurat_obj))
    seurat_obj <- FindNeighbors(seurat_obj, features = VariableFeatures(object = seurat_obj))
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
    options(warn=0)
    seurat_obj<- SetIdent(seurat_obj, value = "cell_ontology_class")
    
    seurat_obj$binary_age <- ifelse(seurat_obj$age %in% c("1m", "3m"),
                                    "young", 
                                    "old")
    
    # Organize Metadata
    meta_meta <- seurat_obj@meta.data %>% dplyr::select(age, 
                                                        mouse.id, 
                                                        binary_age)


    s_per_cell <- setNames(meta_meta$mouse.id, nm = rownames(meta_meta))
    mouse_meta <- meta_meta %>% unique()
    if (do_permute) {
        mouse_vec <- setNames(sample(c("young", "old"), 
                         size = nrow(mouse_meta), replace = T), nm = mouse_meta$mouse.id)
    } else {
        mouse_vec <- setNames(mouse_meta$binary_age, nm= mouse_meta$mouse.id)
    }
    
    # Create object and set coloring
    yeetme <- Cacoa$new(seurat_obj,ref.level = "old", 
                        target.level ="young", 
                        sample.groups= mouse_vec, 
                        sample.per.cell =s_per_cell, graph.name = "UMAP")
    yeetme$cell.groups.palette<- levels(yeetme$cell.groups) %>%  
    {setNames(sample(brewerPalette("Paired")(length(.))), .)}
    
    yeetme
}

cao_obj <- create_cao_from_seurat(seurat_path, permute)

saveRDS(cao_obj, file = output_p))