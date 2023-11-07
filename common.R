processSeurat <- function(s_obj, is_raw = FALSE) {
  seurat_obj <- UpdateSeuratObject(seurat_obj)


  if (is_raw) {
    seurat_obj <- SCTransform(seurat_obj, verbose = F)
  } else {
    seurat_obj <- ScaleData(seurat_obj, do.scale = F)
  }
  seurat_obj <- FindVariableFeatures(object = seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, 
  features = VariableFeatures(object = seurat_obj))
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)

  seurat_obj
}