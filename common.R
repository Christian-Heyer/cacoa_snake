processSeurat <- function(s_obj) {
  seurat_obj <- UpdateSeuratObject(seurat_obj)


  seurat_obj <- SCTransform(seurat_obj, verbose = F)

  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

  seurat_obj
}