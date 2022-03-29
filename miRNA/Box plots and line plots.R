### 
### ---------------
###
### Create: keyanjun
### Date: 2020-07-22 23:24:48
### Email: zhihe_2020@163.com
### Wechat Official Account: keyanjun2020
### Update Log: 2020-05-22  First version
###
### ---------------
setwd("D:\\data learning/scRNA0927/GSE169396_GSE147390/")
rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

cranpackage<-c("ggplot2","ggpubr")

for (pkg in cranpackage){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T)
  }
}

library(ggplot2)
library(ggpubr)
### 载入DEG数据
load("03_差异分析/deg.Rdata")
load("03_差异分析/差异分析_input.Rdata")

### 新建目录/文件夹
dir.create("04_箱线图")

### 提取作图数据【RBM8A的表达量在AD和Normal组的差异-箱线图】
deg1["hsa-miR-218-5p",]
P.Value<-round(deg1["hsa-miR-218-5p","P.Value"],22)
identical(rownames(cli),colnames(exp))

box<-data.frame(row.names = rownames(cli),
                      group=cli$group,
                      miRNA=as.numeric(as.vector(exp["hsa-miR-218-5p",])))
save(box,adj.P.Val,file = "04_箱线图/箱线图_input.Rdata")

### ---------------
rm(list=ls())
load("04_箱线图/箱线图_input.Rdata")
mydata<-box

### 记录2组Sample的数量
n1=length(which(mydata[,1]=="Normal"))
n2=length(which(mydata[,1]=="Osteoarthritis"))

### 正式作图
ggboxplot(mydata, 
          x = "group",
          y = "miRNA",
          size=0.5,
          # facet.by="A",
          fill = "group",
          # ylim = c(0, 6),
          xlab = "",
          ylab ="hsa-miR-218-5p expression values",
          bxp.errorbar=F,
          title = paste0("P.Value= ",P.Value),
          add = c("jitter"))+ 
  #如果想加统计学检验：为Utest，下面这句可运行
  # stat_compare_means(aes(group =group),label="none")+ 
  theme(axis.text.x = element_blank(),   #各个字体大小
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12),
        legend.text= element_text(face="bold", color="black", size=8),
        legend.title = element_text(face="bold", color="black", size=10))+
  scale_x_discrete(breaks=c(0,1),labels=c("",""))+
  theme(legend.background = element_rect(fill = "transparent"))+
  scale_fill_manual(labels = c(paste0("Normal (",n1,")"),paste0("Osteoarthritis (",n2,")")),values = c("#ED7774","#A1D3AC"),name="Group")
ggsave(filename = "04_箱线图/boxplot.pdf",width = 5,height = 5)



### 绘制小提琴图ViolinPlot
ggviolin(mydata, 
         x = "group",
          y = "miRNA",
          size=0.5,
          # facet.by="A",
          fill = "group",
          # ylim = c(0, 6),
          xlab = "",
          width=0.8,
          bxp.errorbar=F,
          ylab ="hsa-miR-218-5p expression values",
          title = paste0("P.Value= ",P.Value),
          add = c("jitter","boxplot","median_iqr"))+ 
  #如果想加统计学检验：为Utest，下面这句可运行
  # stat_compare_means(aes(group =group),label="none")+ 
  theme(axis.text.x = element_blank(),   #各个字体大小
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        legend.text= element_text(face="bold", color="black", size=12),
        legend.title = element_text(face="bold", color="black", size=12))+
  scale_x_discrete(breaks=c(0,1),labels=c("",""))+
  scale_fill_manual(labels = c(paste0("Normal (",n1,")"),paste0("Osteoarthritis (",n2,")")),values = c("#ED7774","#A1D3AC"),name="Group") #
ggsave(filename = "04_箱线图/violinplot.pdf",width = 5,height = 5)


