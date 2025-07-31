# Detailed anaylsis of the cellranger outputs

#load .gbk file of Grelia
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install()
library(BiocManager)

BiocManager::install("Biostrings")
library(Biostrings)
library(Seurat)
library(dplyr)
library(ggplot2)

getwd()
setwd("E:/scTriH2/SEURAT/extracted/1.cellRangerCount/")
raw_matrix_SRR24886407 <- Read10X(data.dir = "SRR24886407/outs/raw_feature_bc_matrix/")
print(head(raw_matrix_SRR24886407))
ncol(raw_matrix_SRR24886407)
nrow(raw_matrix_SRR24886407)
# 1318799 barcodes, 1252 features

# filter with features with less then 2 times of occurence
raw_matrix_SRR24886407 <- raw_matrix_SRR24886407[rowSums(raw_matrix_SRR24886407>0)>= 1, ]
nrow(raw_matrix_SRR24886407) # 355, aligns with the summary
ncol(raw_matrix_SRR24886407)
raw_matrix_SRR24886407 <- raw_matrix_SRR24886407[rowSums(raw_matrix_SRR24886407>0)>= 2, ]
nrow(raw_matrix_SRR24886407) # 179

print(head(raw_matrix_SRR24886407))
print(tail(raw_matrix_SRR24886407))
View(raw_matrix_SRR24886407)

# filter cells by feature count per cell and minimum cell count per feature
  # raw_SRR24886407 <- CreateSeuratObject(counts = raw_matrix_SRR24886407, project = "SRR24886407", min.features = 2)
raw_SRR24886407 <- CreateSeuratObject(counts = raw_matrix_SRR24886407, project = "SRR24886407", min.features = 2, min.cells = 2)
nrow(raw_SRR24886407[["RNA"]]$counts)
ncol(raw_SRR24886407[["RNA"]]$counts)
View(raw_SRR24886407[["RNA"]])

# Meta data metrics, nCount and nFeature, 
head(raw_SRR24886407@meta.data) 

#QC with UMI
VlnPlot(raw_SRR24886407, features="nCount_RNA", log=TRUE)
summary(raw_SRR24886407$nCount_RNA)

VlnPlot(raw_SRR24886407, features="nFeature_RNA", log=FALSE)
summary(raw_SRR24886407$nFeature_RNA)

raw_SRR24886407$GenesPerUMI <- (raw_SRR24886407$nFeature_RNA / raw_SRR24886407$nCount_RNA)
VlnPlot(raw_SRR24886407, features="GenesPerUMI", log=TRUE)
summary(raw_SRR24886407$GenesPerUMI)

raw_SRR24886407$percent_bact <- (raw_SRR24886407$nCount_RNA / sum(raw_SRR24886407$nCount_RNA))*100
VlnPlot(raw_SRR24886407, features="percent_bact", log = TRUE)
summary((raw_SRR24886407$percent_bact))

sum(raw_SRR24886407$nCount_RNA)

View(raw_SRR24886407@meta.data)




# For SRR24886416
library(Seurat)
getwd()
setwd("E:/scTriH2/SEURAT/extracted/1.cellRangerCount/")
raw_matrix_SRR24886416 <- Read10X(data.dir = "SRR24886416/outs/raw_feature_bc_matrix/")
print(head(raw_matrix_SRR24886416))
ncol(raw_matrix_SRR24886416) # 445769
nrow(raw_matrix_SRR24886416) # 1252
raw_matrix_SRR24886416 <- raw_matrix_SRR24886416[rowSums(raw_matrix_SRR24886416>0)>= 1, ]
nrow(raw_matrix_SRR24886416) # 2 features
ncol(raw_matrix_SRR24886416) # 445769 barcodes
raw_matrix_SRR24886416 <- raw_matrix_SRR24886416[rowSums(raw_matrix_SRR24886416>0)>= 3, ]
nrow(raw_matrix_SRR24886416) # 0
ncol(raw_matrix_SRR24886416) # 445769
