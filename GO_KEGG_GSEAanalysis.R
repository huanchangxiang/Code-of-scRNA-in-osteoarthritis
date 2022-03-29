library(ggpubr)
library(patchwork)
library(tidyverse)
library(openxlsx)
library(enrichplot)
library(GOplot)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)
setwd("D:\\data learning/scRNA0927/GSE169396_GSE147390/")
load("multi/experiment.aggregate.RData")
sub.markers <- FindMarkers(experiment.aggregate,
                           ident.1 = "Osteoarthritis",
                           ident.2 = "Normal",group.by = "orig.ident")
dge.cluster <- sub.markers
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.05&avg_log2FC > 1)
sig_dge.cluster$gene  <- rownames(sig_dge.cluster)
sig_dge.cluster <- bitr(sig_dge.cluster$gene,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = "org.Hs.eg.db")
#画圈图代码
sig_RNA <- subset(dge.cluster, p_val_adj<0.05&avg_log2FC > 1)
sig_RNA$SYMBOL <- rownames(sig_RNA)
Enrich <- sig_dge.cluster %>%
  inner_join(sig_RNA, by = "SYMBOL")
head(Enrich)

Enrichment <- Enrich %>%
  as.data.frame()%>%
  dplyr::select(SYMBOL,ENTREZID,avg_log2FC)

Enrichment <- Enrichment[!duplicated(Enrichment[,1]),]
gene <- Enrichment$ENTREZID
#差异基因GO富集分析
ego_ALL <- enrichGO(gene          = sig_dge.cluster$ENTREZID,
                    
                    OrgDb         = "org.Hs.eg.db",
                    keyType       = "ENTREZID",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                 
                    qvalueCutoff  = 0.05,
                    readable = T)
ego_all <- data.frame(ego_ALL)
ego_CC <- enrichGO(gene          = sig_dge.cluster$ENTREZID,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = "ENTREZID",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   
                   qvalueCutoff  = 0.05,
                   readable = T)
ego_MF <- enrichGO(gene          = sig_dge.cluster$ENTREZID,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = "ENTREZID",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable = T)
ego_BP <- enrichGO(gene          = sig_dge.cluster$ENTREZID,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = "ENTREZID",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   
                   qvalueCutoff  = 0.05,
                   readable = T)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
GO_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
GO_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
GO_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
GO_ALL <- barplot(ego_ALL,showCategory = 10) + ggtitle("barplot for ALL")

GO_BP1 <- dotplot(ego_BP,showCategory = 20) + ggtitle("GO for Biological process")+ scale_color_continuous(low='#FDAE6BFF', high="#C7E9C0FF")
GO_CC1 <- dotplot(ego_CC,showCategory = 20) + ggtitle("GO for Cellular component")+ scale_color_continuous(low='#FDAE6BFF', high="#C7E9C0FF")
GO_MF1 <- dotplot(ego_MF,showCategory = 20) + ggtitle("GO for Molecular function")+ scale_color_continuous(low='#FDAE6BFF', high="#C7E9C0FF")
GO_ALL1 <- dotplot(ego_ALL,showCategory = 20 ) + ggtitle("barplot for ALL")+ scale_color_continuous(low='#FDAE6BFF', high="#C7E9C0FF")

plotc <- GO_BP/GO_CC/GO_MF/GO_ALL/GO_BP1/GO_CC1/GO_MF1/GO_ALL1
ggsave("GO富集分析.png", plotc, width = 15,height = 50,limitsize = F)
rm(GO_BP,GO_CC,GO_MF,GO_ALL,GO_BP1,GO_CC1,GO_MF1,GO_ALL1,ego_ALL,ego_BP,ego_CC,ego_MF,plotc)



#差异基因KEGG分析
enrichKK <- enrichKEGG(
                   gene = sig_dge.cluster$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   use_internal_data = T
)
p1 <- barplot(enrichKK,showCategory=20)
p2 <- dotplot(enrichKK, showCategory=20)
plotc = p1/p2
ggsave("KEGG富集分析.png", plot = plotc, width = 12, height = 10,limitsize = F)
rm(plotc)

#弦图GO/KEGG（需要替换）
go <- data.frame(Category = ego_BP@ontology,
                 ID = ego_BP@result$ID,
                 Term = ego_BP@result$Description  ,
                 Genes = gsub("/", ", ", ego_BP@result$geneID),
                 adj_pval = ego_BP@result$p.adjust)
genelist <- data.frame(ID = Enrichment$SYMBOL,
                       logFC = Enrichment$avg_log2FC)
row.names(Enrichment)=Enrichment[,1]
circ <- circle_dat(go, genelist)

termNum = 10  #选取要选择哪几个途径                            
geneNum = nrow(genelist)                      
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])


GOChord(chord,
        space = 0.02,        
        gene.order = 'pvalue',    
        gene.space = 0.35,      
        gene.size = 2 )

#KEGG弦图
KEGG <- data.frame(Category = enrichKK@ontology,
                   Term = enrichKK@result$Description,
                 ID = enrichKK@result$ID,
                 Genes = gsub("/", ", ", enrichKK@result$geneID),
                 adj_pval = enrichKK@result$p.adjust)
genelist <- data.frame(ID = Enrichment$ENTREZID,
                       logFC = Enrichment$avg_log2FC)
row.names(Enrichment)=Enrichment[,2]
circ <- circle_dat(KEGG, genelist)

termNum = 1 #选取要选择哪几个途径                            
geneNum = nrow(genelist)                      
chord <- chord_dat(circ, genelist[1:geneNum,], KEGG$Term[1:termNum])
SYMBOL <- bitr(rownames(chord),
               fromType = "ENTREZID",
               toType = "SYMBOL",
               OrgDb = "org.Hs.eg.db")
rownames(chord) <- SYMBOL$SYMBOL
KEGG_PLOT <- GOChord(chord,
        space = 0.02,        
        gene.order = 'pvalue',    
        gene.space = 0.2,      
        gene.size = 5,
        process.label = 16,
        lfc.col=c('firebrick3', 'white','royalblue3'),
        ribbon.col=c("#CC0C00FF","#5C88DAFF","#84BD00FF","#D595A7FF","#00B5E2FF",
                     "#00AF66FF","#6EE2FFFF","#F7C530FF","#95CC5EEF","#FD7446FF"))


#####弦图调节例子
GOChord(chord, title="GOChord plot",#标题设置
        space = 0.02, #GO term处间隔大小设置
        limit = c(3, 5),#第一个数值为至少分配给一个基因的Goterm数，第二个数值为至少分配给一个GOterm的基因数
        gene.order = 'logFC', gene.space = 0.25, gene.size = 5,#基因排序，间隔，名字大小设置
        lfc.col=c('firebrick3', 'white','royalblue3'),##上调下调颜色设置
        #ribbon.col=colorRampPalette(c("blue", "red"))(length(EC$process)),#GO term 颜色设置
        ribbon.col=brewer.pal(length(EC$process), "Set3"))#GO term 颜色设置

        
#拼图
plota <- ggarrange(GO_BP1,GO_CC1,GO_MF1,
                  ncol =2,nrow = 2,widths =1000,heights = 50,
                  font.label = list(size=30,color="black"))

ggsave("multi/GO_KEGG_GSEA/KEGG小于-1.jpg",KEGG_PLOT,
       dpi = 1200,limitsize = FALSE,width = 21,height = 15)




#GSEA分析
# 读入数据
df <- openxlsx::read.xlsx("multi/DEGS.xlsx",sheet = 2)
head(df)
# ENTREZ ID注释
df.id <- bitr(df$genes, # 转换的列是df数据框中的genem基因名列
              fromType = "SYMBOL", # 需要转换ID类型
              toType = "ENTREZID", # 转换成的ID类型
              OrgDb = "org.Hs.eg.db" ) # 对应的物种
head(df.id)
names(df) <- c("SYMBOL","logFC") #  重新命名df的列名
easy.df<-merge(df,df.id,by="SYMBOL",all=F) # 以SMBOL列为依据，合并df和df.id
head(easy.df)
sortdf<-easy.df[order(easy.df$logFC, decreasing = T),]
head(sortdf)
gene.expr = sortdf$logFC #把logFC按照从大到小提取出来
head(gene.expr)
names(gene.expr) <- sortdf$ENTREZID#给上面提取的logFC对应上ENTREZID
head(gene.expr)
#GSEA中KEGG分析
kk <- gseKEGG(geneList     = gene.expr,
               organism     = 'hsa',
               pvalueCutoff = 0.05)
#按照enrichment score从高到低排序，便于查看富集通路
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk)
dim(sortkk)
write.table(sortkk,"multi/GO_KEGG_GSEA/GSEAKEGG.txt",sep = "\t",quote = F,col.names = T,row.names = F)
gseaplot2(kk,#数据
          c("hsa03010"),#画那一列的信号通路
          title = "Ribosome(P=1.393094e-05)",#标题
          base_size = 10,#字体大小
          color = "#483D8B",#线条的颜色
          pvalue_table = F,#加不加p值
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          ES_geom="line")#是用线，还是用d点

ggsave("multi/GO_KEGG_GSEA/GSEA.jpg",GSEA,
       dpi = 1200,limitsize = FALSE,width = 20,height = 10)
#gseaplot2参数
gseaplot2(gse.KEGG,
          title = "Olfactory transduction",  #设置title
          "hsa04740", #绘制hsa04740通路的结果
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:2, #展示上2部分
          pvalue_table = T) # 显示p值