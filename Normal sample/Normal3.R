#单基因画图
library(Seurat)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggsci)
    




experiment.aggregate@meta.data$Barcode <-rownames(experiment.aggregate@meta.data)
x<-as.data.frame(experiment.aggregate@reductions$tsne@cell.embeddings) 
x$Barcode <-rownames(x)
names(experiment.aggregate@meta.data)[colnames(experiment.aggregate@meta.data)=="RNA_snn_res.0.5"] <- "Cluster"
y <-data.frame(experiment.aggregate@meta.data[,c('Barcode','seurat_clusters','orig.ident',"Cluster")])
lab <-merge(x,y,barcode='Barcode')
color<-c("#00AFBB", "#E7B800", "#FC4E07") 
##按照cluster id进行排序


lab$Cluster <- ifelse(lab$Cluster == "8","Normal Osteoblasts","Others")

#TSNE画图
ggplot(lab,mapping = aes(x=tSNE_1,
                         y=tSNE_2,
                         col=Cluster))+
  geom_point(size=1)+
  theme_classic()+
  scale_fill_discrete(name= "Cluster")+
  theme(legend.text = element_text(size = 20),legend.position = "top")+
  guides(color = guide_legend(override.aes = list(size=7)))+
  
  scale_color_manual(breaks = c("Normal Osteoblasts", "Others"),values = c("#DC0000FF", "#D9D9D9FF"))+
  theme(axis.line = element_line(size=1.2, colour = "black"))  



