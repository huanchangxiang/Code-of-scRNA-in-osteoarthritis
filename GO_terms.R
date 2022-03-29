

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Nebulosa")
library("Nebulosa")
library("Seurat")


plot_density(experiment.aggregate, c("S100A9", "S100A8", "IGKC", "IGHG4", "IGHA1", "IGHG1", "IGLC2", 
                                     "IGLC3", "DEFA3", "CXCL8", "IGHG2"), joint = TRUE,size = 3)[[12]]

plot_density(experiment.aggregate,"IGKC")
