progress(SRR, progress.bar = TRUE)
}

# retrieve the data with this name)
print(ls(pattern = "raw_data_"))

# define tasks for the loop in a function named "do.seurat"
do.seurat <- function(obj) {
  setwd("E:/scTriH2/SEURAT/extracted/1.cellRangerCount/")
  print("Load the symbiont cellranger count data.")
  obj.data <- Read10X(obj)
  head(obj.data)
  dir.create(obj)
  setwd(obj)
  print("Create Seurat Object to store data and analysis, print the first 10 lines.")
  obj <- CreateSeuratObject(counts = scTriH2.SRR24886407.data, project = "scTriH2", min.cells = 3, min.features = 1)
  obj[1:10]
  print("Check the cellranger count matrix.")
  obj.dense.size <- object.size(as.matrix(scTriH2.SRR24886407.data))
  obj.dense.size # 98288936 bytes
  obj.sparse.size <- object.size((scTriH2.SRR24886407.data))
  obj.sparse.size # 1186088 bytes
  obj.dense.size/obj.sparse.size # 82.9 bytes
  print("QC with mt percentage.") 
  scTriH2.SRR24886407[["percent.mt"]]<- PercentageFeatureSet(scTriH2.SRR24886407, pattern = "^MT-")
  head(scTriH2.SRR24886407@meta.data, 10)
  
  # Violin plot? 
  
  Print("Visualization of feature-feature realtionships.")
  plot1 <- FeatureScatter(scTriH2.SRR24886407, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(scTriH2.SRR24886407, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2 # Very low mt.percent and nFeature_RNA due to the scRNA seq library method
  
}

for (i in 1: length(samples)) {
  samples[[i]] <- do.seurat(i)
}

scTriH2.data<- c(list.files("E:/scTriH2/SEURAT/extracted/1.cellRangerCount/", full.names = TRUE))
scTriH2.data
for (SRR in 1:length(scTriH2.data)) {
  scTriH2.{SRR} <- Read10X("extracted/1.cellRangerCount/{SRR}")
  scTriH2.{SRR}[1:30] # read the first 30 cells
  scTriH2 <- CreateSeuratObject(counts = scTriH2.data, project = "scTriH2", min.cells = , min.features = )
}