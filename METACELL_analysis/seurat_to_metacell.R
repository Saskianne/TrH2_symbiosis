
setwd("E:/scTriH2/metacell/")

library("mc2d")
library(Matrix)
library(Seurat)
library(ggplot2)
library(sctransform)


# Check for duplicate gene names


any(duplicated(colnames(merged_seurat)))
any(duplicated(rownames(merged_seurat)))


merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
# Save as RDS file
saveRDS(merged_seurat, file = "E:/scTriH2/SEURAT/extracted/1.cellRangerCount/Seurat_outputs/merged_seurat_object_1.rds")
cat("Merged Seurat object saved successfully.\n")

merged_seurat_UMAP <- FindNeighbors(merged_seurat, dims = 1:20)
merged_seurat_UMAP <- FindClusters(merged_seurat_UMAP)
merged_seurat_UMAP <- RunUMAP(merged_seurat_UMAP, dims = 1:20)
DimPlot(merged_seurat_UMAP, reduction = "umap")



