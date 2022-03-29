#### 1. 加载分析使用的工具包 ####
#每次运行都需要先加载相关的包才能调用后面的分析脚本
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
memory.limit(102400)



experiment.data <- Read10X("bm1/")
colnames(experiment.data) <- paste(colnames(experiment.data),"Osteoarthritis",sep = "_")
sam.name <- "Osteoarthritis"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}


experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "Osteoarthritis", 
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")


save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_Osteoarthritis.RData"))
experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^MT-")

#画图的代码需要选中从pdf一直到dev.off()的所有代码，一起运行
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

##QC：统计基因数，RNA，线粒体基因分布
gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c("Count_Gene","Count_RNA","MT_percent")
View(freq.combine)
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)


##QC：基因数与线粒体基因以及RNA数量的分布相关性
plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)


cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")


experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 2000)

#展示标准化之后的整体表达水平
top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)


experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("percent.mt"))

#PCA计算
#features:通常使用前面过程中计算得到的高变异基因，所以调整高变异基因的数量会改变PCA的结果
#建议实际操作中可以修改测试
experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F)

#PCA结果展示-1
#Visualize top genes associated with reduction components
pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
#2D点图展示 PC1&PC2，每个点代表一个细胞
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()

#PCA结果展示-3
pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"),width = 5,height = 4)
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()


#确定用于细胞分群的PC
dim.use <- 1:20

#### 9. 细胞分群TSNE算法 ####
#先用KNN聚类算法得到最优类数量（SNN）及细胞之间的相似性统计
experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
#基于细胞之间的相似性计算具体的cluster
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)
#TSNE算法聚类（降维）
experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use, do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEplot_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "tsne",label = T)
dev.off()

#### 10. 计算marker基因 ####
all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_",max(dim.use),"PC.txt"),sep="\t",quote = F)



b <- subset(experiment.aggregate,seurat_clusters=="6")
zero_data <- experiment.data[,rownames(b@meta.data)]
three_data <- experiment.data[,rownames(b@meta.data)]
six_data <- experiment.data[,rownames(b@meta.data)]

colnames(zero_data) <- paste(colnames(zero_data),"zero",sep = "_")
colnames(three_data) <- paste(colnames(three_data),"three",sep = "_")
colnames(six_data) <- paste(colnames(six_data),"six",sep = "_")
Osteoarthritis_data <- cbind(zero_data,three_data,six_data)
save(Osteoarthritis_data,file = "Osteoarthritis/.Rdata")
