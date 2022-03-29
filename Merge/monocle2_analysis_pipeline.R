# 该代码用于进行拟时序分析

#加载分析使用的包
library(Seurat)
library(monocle)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
memory.limit(102400)
##使用monocle2进行拟时序分析
#构造表达及注释数据，提取CCA之后的数据
exp.matrix<-as(as.matrix(experiment.aggregate@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann)<-rownames(exp.matrix)
exp_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-experiment.aggregate@meta.data
rownames(sample_ann)<-colnames(exp.matrix)
exp_pd<-new("AnnotatedDataFrame", data =sample_ann)

#生成monocle对象
exp.monocle<-newCellDataSet(exp.matrix,phenoData =exp_pd,featureData =exp_fd,expressionFamily=negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#计算sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)
#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因#计算出差异基因并且会排序
diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~seurat_clusters") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.01))
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))
  
#修改monocle对象中的列名示例
names(pData(exp.monocle))[names(pData(exp.monocle))=="seurat_clusters"]="Cluster"

#将不同分组情况的拟时序轨迹图画到一起
plot1<-plot_cell_trajectory(exp.monocle, color_by = "Cluster",cell_size=1)
ggsave("./multi/trajectory_plot1.pdf",plot = plot1,width = 16,height = 10)
plot2<-plot_cell_trajectory(exp.monocle, color_by = "orig.ident",cell_size=1)
ggsave("./multi/trajectory_plot2.pdf",plot = plot2,width = 16,height = 10)
plot_cell_trajectory(exp.monocle, color_by = "RNA_snn_res.0.5",cell_size=4)+
  theme(axis.line = element_line(size=50, colour = "black"),legend.title = element_blank(),
        legend.text = element_text(size = 20))
ggsave("./multi/trajectory_plot3.pdf",plot = plot3,width = 16,height = 10)
plot4<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size=1)
ggsave("./multi/trajectory_plot4.pdf",plot = plot4,width = 16,height = 10)
plot5<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size=1)
ggsave("./multi/trajectory_plot5.pdf",plot = plot5,width = 16,height = 10)


