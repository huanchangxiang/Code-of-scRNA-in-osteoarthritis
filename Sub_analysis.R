####  Code Description              ####
#---  1. Written by WoLin @ 2019.03.15，last update 19.08.26 ---#
#---  2. private use for scRNA    ---#

# 修改具体类别的名字(需要基于marker基因来定义)
experiment.merged <- RenameIdents(
  object = OA.aggregate,
  "0" = "Pre-osteoblast",
  "1" = "Pre-osteoblast",
  "3" = "Pre-osteoblast",
  "2" = "Mature osteoblast",
  "4" = "Mature osteoblast",
  "5" = "Undetermined osteoblasts"
)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_Rename_",max(dim.use),"PC.pdf"),width = 10,height = 5)
DimPlot(object = experiment.merged, pt.size=0.5,label.size = 5,
        reduction = "tsne",label = T) +
  ggsci::scale_color_igv()
dev.off()

###替换
experiment.aggregate@meta.data$orig.ident<-as.character(experiment.aggregate@meta.data$orig.ident)
class(experiment.aggregate@meta.data$orig.ident)
experiment.aggregate@meta.data$orig.ident[experiment.aggregate@meta.data$orig.ident=="Normal4"]<-"Normal"
table(experiment.aggregate@meta.data$orig.ident)
# 计算特定两组细胞之间的差异基因
sub.markers <- FindMarkers(experiment.aggregate,
                           ident.1 = "Osteoarthritis",
                           ident.2 = "Normal",group.by = "orig.ident")
sub.markers$genes <- rownames(sub.markers)
write.table(sub.markers,
            file=paste0("./",sam.name,"/",sam.name,"DEGS_sample",max(dim.use),"PC.txt"),
            sep="\t",quote = F,row.names = F)
