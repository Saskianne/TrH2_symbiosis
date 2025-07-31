
# Download and install the Rtools44 to enable the install of Seurat
writeLines('PATH="C:/Program Files/R/R-4.4.2/library/rtools44/usr/bin;${PATH}"', con = "~/.Renviron")

# Restart the R or Rstudio and check the RTool is found
system('g++ --version')
  # This should resolve th issue and allow R to correctly access Rtools. 

# Install Seurat and check
install.packages('Seurat')
library(Seurat)

# Add the recommended installations 
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remote packages
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)


