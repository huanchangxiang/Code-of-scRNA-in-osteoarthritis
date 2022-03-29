#### 1. 加载分析使用的工具包 ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
memory.limit(102400)

setwd("D:/data learning/scRNA0927/GSE169396_GSE147390/Normal/")
#载入数据
exp2<-Read10X("bm2/")
colnames(exp2) <- paste(colnames(exp2),"Normal1",sep = "_")

exp3<-Read10X("bm3/")
colnames(exp3) <- paste(colnames(exp3),"Normal2",sep = "_")

exp4<-Read10X("bm4/")
colnames(exp4) <- paste(colnames(exp4),"Normal3",sep = "_")

exp5<-Read10X("bm5/")
colnames(exp5) <- paste(colnames(exp5),"Normal4",sep = "_")

experiment.data <- cbind(exp2,exp3,exp4,exp5) 


sam.name <- "Normal"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "multi",
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")


save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_normal.RData"))



experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^MT-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()





gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

##QC：基因数与线粒体基因以及RNA数量的分布相关性
plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)



cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")

#### 6. 表达量标准化 ####
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
  vars.to.regress = c("orig.ident","percent.mt"))


#PCA计算
experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F,npcs = 50)

#PCA结果展示-1
pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()


#碎石图
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(experiment.aggregate,ndims = 40)
dev.off()

#确定用于细胞分群的PC
dim.use <- 1:20


#### 9. 细胞分群TSNE算法 ####
#TSNE算法
experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)

experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use, 
                                do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.5_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = experiment.aggregate, pt.size=0.5,label = T)
dev.off()

#按照数据来源分组展示细胞异同--画在一张图中
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = experiment.aggregate, 
        group.by="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

#按照数据来源分组展示细胞异同--画在多张图中
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_slipt_",max(dim.use),"PC.pdf"),width = 8,height = 4)
DimPlot(object = experiment.aggregate, 
        split.by ="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

table(experiment.aggregate@meta.data$orig.ident)

#### 10. 计算marker基因 ####
#这一步计算的时候可以把min.pct以及logfc.threshold调的比较低，然后再基于结果手动筛选
all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, 
                              min.pct = 0.1, logfc.threshold = 0.1)
write.table(all.markers,
            file=paste0("./",sam.name,"/",sam.name,"_total_marker_genes_tsne_",max(dim.use),"PC.txt"),
            sep="\t",quote = F,row.names = F)

a <- subset(experiment.aggregate,seurat_clusters=="8")
normal_data <- experiment.data[,rownames(a@meta.data)]
save(normal_data,file = "Normal/normal_data.Rdata")


