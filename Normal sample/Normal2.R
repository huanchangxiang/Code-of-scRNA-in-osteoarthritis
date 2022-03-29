#### 1. 加载分析使用的工具包 ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
memory.limit(102400)

setwd("D:\\data learning/scRNA0927/GSE169396_GSE147390/")
load("Normal/normal_data_cluster8.Rdata")
sam.name <- "Normal_cluster8"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

normal.aggregate <- CreateSeuratObject(
  normal_data,
  project = "multi", 
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")



normal.aggregate[["percent.mt"]] <- PercentageFeatureSet(normal.aggregate, 
                                                     pattern = "^MT-")
normal.aggregate[["percent.RP"]] <- PercentageFeatureSet(normal.aggregate, 
                                                     pattern = "^RP[S/L]")
gene.freq <- do.call("cbind", tapply(normal.aggregate@meta.data$nFeature_RNA,normal.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(normal.aggregate@meta.data$nCount_RNA,normal.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(normal.aggregate@meta.data$percent.mt,normal.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
RP.freq <- do.call("cbind", tapply(normal.aggregate@meta.data$percent.mt,normal.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq,RP.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(RP.freq),"RP",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
rm(gene.freq,rna.freq,mt.freq,RP.freq)



cat("Before filter :",nrow(normal.aggregate@meta.data),"cells\n")
normal.aggregate <- subset(normal.aggregate, 
                       subset = 
                         nFeature_RNA > 200 & 
                         nCount_RNA > 500 & 
                         nCount_RNA < 20000 &
                         percent.mt < 25 &
                         percent.RP < 25)
cat("After filter :",nrow(normal.aggregate@meta.data),"cells\n")


normal.aggregate <- NormalizeData(normal.aggregate, 
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
normal.aggregate <- FindVariableFeatures(normal.aggregate, 
                                     selection.method = "vst",
                                     nfeatures = 2000)

normal.aggregate <- ScaleData(
  object = normal.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("percent.mt","percent.RP"))

#PCA计算
normal.aggregate <- RunPCA(object = normal.aggregate, 
                       features = VariableFeatures(normal.aggregate),
                       verbose = F,npcs = 50)
ElbowPlot(normal.aggregate,ndims = 40)
dim.use <- 1:20
#TSNE画图
normal.aggregate <- FindNeighbors(normal.aggregate, dims = dim.use)
normal.aggregate <- FindClusters(normal.aggregate, resolution = 0.5)

normal.aggregate <- RunTSNE(normal.aggregate, dims = dim.use, 
                        do.fast = TRUE)
DimPlot(object = normal.aggregate, pt.size=0.5,label = T)


#### 10. 计算marker基因 ####
#这一步计算的时候可以把min.pct以及logfc.threshold调的比较低，然后再基于结果手动筛选
all.markers <- FindAllMarkers(normal.aggregate, only.pos = TRUE, 
                              min.pct = 0.2, logfc.threshold = 0.2)
write.table(all.markers,
            file=paste0("./Normal_cluster8//total_marker_genes_tsne_PC.txt"),
            sep="\t",quote = F,row.names = F)
