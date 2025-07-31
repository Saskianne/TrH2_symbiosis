#SEURAT scRNA transcriptome analysis

## LOOP 
# preparation
library(svMisc)
library(dplyr)
library(Seurat)
library(patchwork)
library(BiocManager)
library(Biostrings)
library(Seurat)
library(Matrix)

parent_directory <- "E:/scTriH2/SEURAT/extracted/1.cellRangerCount"
setwd(parent_directory)
# list of SRRs to work with
SRR_list_full <- list.dirs(parent_directory, full.names = FALSE, recursive = FALSE)
SRR_list_full <- SRR_list_full[ grepl("SRR*", SRR_list_full) ]
SRR_list <- SRR_list_full[SRR_list_full != "SRR24886387"]
cat(SRR_list)
# output files will be stored in this directory
dir.create("Seurat_outputs", showWarnings = FALSE)

#loop for the SEURAT analysis
cat("Loop starts \n")
for (SRR in SRR_list) {
  cat("Working with", SRR, "\n")
  base_dir <- file.path(parent_directory,"Seurat_outputs")
  setwd(base_dir)
  
  #create an output directory for the SRR and a 
  SRR_output <- paste(SRR,"output",sep = "_")
  dir.create(SRR_output, showWarnings = FALSE)
  cat("Output folder of ", SRR, " has been created \n")
  output_dir <- file.path(base_dir, SRR_output) 
  setwd(output_dir)
  
  # create an output file
  outfile_txt <- paste0(SRR,"_output.txt")
  file.create(outfile_txt)
  output_path <- file.path(output_dir, outfile_txt)
  output_path <- sub("[/\\]+$", "", output_path)
  outfile_conn <- file(output_path, "a")
  cat("The output file of ", SRR, "has been created \n")
  close(outfile_conn)
  
  raw_dir <- file.path(parent_directory, SRR, "outs/raw_feature_bc_matrix/")
  # raw_matrix read out
  cat("Start the seurat run of ", SRR, "\n")
  raw_matrix_name <- paste0("raw_matrix_", SRR)
  raw_matrix <- Read10X(data.dir = raw_dir)
  cat("Created the matrix: ", raw_matrix_name, "\n")

  # inputs
  writeLines(c(SRR, "\n First 5 lines of the matrix: \n", capture.output(head(raw_matrix)), "\n The Original column no. and row no. are: \n", ncol(raw_matrix), ",", nrow(raw_matrix), "\n"), outfile_conn)
  close(outfile_conn)
  
  # data processing:
  raw_matrix_filtered <- raw_matrix[rowSums(raw_matrix > 0) >= 2, ]
  assign(raw_matrix_name, raw_matrix_filtered)
  
  object_name <- paste0("object_", SRR)
  seurat_object <- CreateSeuratObject(counts = raw_matrix_filtered, project = SRR, min.features = 2, min.cells = 2)
  assign(object_name, seurat_object)
  cat("Seurat object", object_name, "has been created\n")
  temp_object <- get(object_name)
  temp_object$GenesPerUMI <- (temp_object$nFeature_RNA / temp_object$nCount_RNA)
  
  writeLines(c("The filtered column no. and row no. are: \n", ncol(temp_object[["RNA"]]$counts), "\n", nrow(temp_object[["RNA"]]$counts)), "\n", outfile_conn)
  writeLines(head(seurat_object@metadata))
  writeLines(c("Filtered gene count:", temp_object$nFeature_RNA, "Filtered cell count:", temp_object$nCount_RNA),sep = "\n", outfile_conn )
  writeLines(c("Gene counts per UMI", temp_object$GenesPerUMI), sep = "\n", outfile_conn)
  close(outfile_conn)
  
  seurat_object <- subset(temp_object, subset = nFeature > 2)
  # extract the filtered matrix
  filtered_matrix <- seurat_object@assay$RNA@counts 
  #  save like 10X output (in 3 files)
  setwd(output_dir)
  writeMM(filtered_matrix, file = file.path(output_dir, "matrix.mtx"))
  write.table(rownames(filtered_matrix), file = file.path(output_dir, "features.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(filtered_matrix), file = file.path(output_dir, "barcodes.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines("The filtered matrix is saved in 10X output_format.", outfile_conn)
  close(outfile_conn)
  
  # also save as RDS file
  saveRDS(seurat_object, file = paste0("filtered_", object_name, ".rds"))
  
  SRR_paths <- c(SRR_paths, output_dir)
  
  writeLines(c("The process for ", SRR, " has been completed."), outfile_conn)
  close(outfile_conn)

cat("Run for every SRR in SRR_list are done \n")
cat("Start working with merging \n")
# merge multiples runs into one seurat file  
matrix_paths <- list.dirs(SRR_paths, full.names=TRUE)
  # matrix_list <- list.files(path = SRR_list, full.names = TRUE)

seurat_list <- lapply(matrix_paths, function(x) {
  mat <- Read10X(data.dir = x)
  CreateSeuratObject(counts = mat, project = basename(x), min.cells = 2, min.features = 2)
})

summary_matrix <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = basename(matrix_paths))
}

cat("merging done \n")
  
  
  
  
  
  
  
  
  
  
  
  

