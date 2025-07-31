# Load required libraries
library(Seurat)
library(Matrix)

#provide seurat_list
options(Seurat.object.assay.version = 'v3')

# Extract count matrices from all Seurat objects
count_matrices <- lapply(seurat_list, function(obj) GetAssayData(obj, assay = "RNA", slot = "counts"))

# Find all unique cell barcodes and gene names across all objects
all_cells <- unique(unlist(lapply(count_matrices, colnames)))
all_features <- unique(unlist(lapply(count_matrices, rownames)))

# Ensure all matrices have the same dimensions by filling missing cells/features with zeros
count_matrices <- lapply(count_matrices, function(mat) {
  # Create an empty sparse matrix with the full feature and cell list
  new_mat <- Matrix(0, nrow = length(all_features), ncol = length(all_cells), sparse = TRUE)
  rownames(new_mat) <- all_features
  colnames(new_mat) <- all_cells
  
  # Fill in existing values from original matrix
  existing_features <- rownames(mat) %in% all_features
  existing_cells <- colnames(mat) %in% all_cells
  new_mat[rownames(mat)[existing_features], colnames(mat)[existing_cells]] <- mat[existing_features, existing_cells]
  
  return(new_mat)
})

# Sum all matrices together
merged_counts <- Reduce(`+`, count_matrices)
setwd("/mnt/e/scTriH2/metacell")

# Create a new Seurat object from the summed count matrix
merged_seurat <- CreateSeuratObject(counts = merged_counts, project = "scTriH2_merged_sum")

# Check the new object
print(merged_seurat)


