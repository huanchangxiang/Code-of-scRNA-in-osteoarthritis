#单基因画图
library(Seurat)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggsci)


setwd("D:\\data learning/scRNA0927/GSE169396_GSE147390/")
load("multi/experiment.aggregate.RData")

experiment.aggregate@meta.data$Barcode <-rownames(experiment.aggregate@meta.data)
x<-as.data.frame(experiment.aggregate@reductions$tsne@cell.embeddings) 
x$Barcode <-rownames(x)
names(experiment.aggregate@meta.data)[colnames(experiment.aggregate@meta.data)=="RNA_snn_res.0.5"] <- "Cluster"
y <-data.frame(experiment.aggregate@meta.data[,c('Barcode','seurat_clusters','orig.ident',"Cluster")])
lab <-merge(x,y,barcode='Barcode')
color<-c("#00AFBB", "#E7B800", "#FC4E07")   ###R绘图的默认色
lab$seurat_clusters <-factor(lab$seurat_clusters,levels=0:12)   ##按照cluster id进行排序,如果是文字，




#TSNE画图
ggplot(lab,mapping = aes(x=tSNE_1,
                         y=tSNE_2,
                         col=Cluster))+
  geom_point(size=3)+
  stat_ellipse(aes(fill=Cluster),
               geom = "polygon",
               linetype =1,        ###圆圈线的类型
               size=0.3,              ###圆圈线的粗细
               alpha=1/5)+
  geom_text(aes(x=-22,y=34),label= "Pre-osteoblasts", size=8,color= "#95CC5EFF") +
  geom_text(aes(x=25,y=30),label= "Undetermined osteoblasts 2", size=8,color= "#00FFFFFF")+
  geom_text(aes(x=23,y=-33),label= "Mature osteoblasts", size=8,color= "#CC0C00FF")+ 
  geom_text(aes(x=-33,y=25),label= "Undetermined osteoblasts 1", size=5,color= "#CC99FFFF")+
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size = 20),legend.position = "top")+
 scale_color_manual(breaks = c("Mature osteoblasts", "Pre-osteoblasts",
                               "Undetermined osteoblasts 1","Undetermined osteoblasts 2"),
                   values = c("#CC0C00FF", "#95CC5EFF","#00FFFFFF","#CC99FFFF"))+
  theme(axis.line = element_line(size=1.5, colour = "black"))





