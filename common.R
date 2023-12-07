processSeurat <- function(s_obj, is_raw = FALSE) {
  seurat_obj <- UpdateSeuratObject(seurat_obj)


  if (is_raw) {
    seurat_obj <- SCTransform(seurat_obj, verbose = F)
  } else {
    seurat_obj <- ScaleData(seurat_obj, do.scale = T)
  }
  # Filter genes expressed in less than 20 cells

  seurat_obj <- FindVariableFeatures(object = seurat_obj, nfeatures = 2000, selection.method = "vst")
  seurat_obj <- RunPCA(seurat_obj, npcs = 15) 
  features = VariableFeatures(object = seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:15,k.param = 20)
  seurat_obj <- RunUMAP(seurat_obj, min.dist = 0.5, spread = 1, dims=1:15, n.neighbors = 20)

  seurat_obj
}