processSeurat <- function(s_obj) {
  seurat_obj <- UpdateSeuratObject(seurat_obj)


  seurat_obj <- SCTransform(seurat_obj, verbose = F)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, 
  features = VariableFeatures(object = seurat_obj))
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)

  seurat_obj
}