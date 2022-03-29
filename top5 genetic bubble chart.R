
all.markers$cluster<-as.character(all.markers$cluster)
class(all.markers$gene)
all.markers$cluster[all.markers$cluster=="0"]<-"Pre-osteoblasts"
all.markers$cluster[all.markers$cluster=="1"]<-"Pre-osteoblasts"
all.markers$cluster[all.markers$cluster=="2"]<-"Pre-osteoblasts"
all.markers$cluster[all.markers$cluster=="3"]<-"Mature osteoblasts"
all.markers$cluster[all.markers$cluster=="4"]<-"Mature osteoblasts"
all.markers$cluster[all.markers$cluster=="5"]<-"Mature osteoblasts"
all.markers$cluster[all.markers$cluster=="7"]<-"Undetermined osteoblasts1"
all.markers$cluster[all.markers$cluster=="6"]<-"Undetermined osteoblasts2"




dotplot_f_data_1 <- DotPlot(experiment.aggregate,group.by = "RNA_snn_res.0.5",
                            features = top5$gene,split.by = "orig.ident")$data
head(dotplot_f_data_1)
dim(dotplot_f_data_1)


dotplot_f_data <- dotplot_f_data_1
head(dotplot_f_data)
dim(dotplot_f_data)
table(dotplot_f_data$id)

p_dotplot_f <- DotPlot(experiment.aggregate, features = top5$gene,col.min=-2, col.max=2)
p_dotplot_f$data <- dotplot_f_data
p_dotplot_f <- p_dotplot_f + coord_flip()
p_dotplot_f <- p_dotplot_f+ scale_color_gradient2(high="#D43F3AFF",mid ="#46B8DAFF", midpoint = 0) + theme_classic()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8,size = 10),
        axis.text.y = element_text(size = 12)) 

p_dotplot_f


