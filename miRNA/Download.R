### 
### ---------------
###
### Create: keyanjun
### Date: 2020-07-04 23:24:48
### Email: zhihe_2020@163.com
### Wechat Official Account: keyanjun2020
### Update Log: 2020-05-22  First version
###
### ---------------

### 环境设置
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

### 安装必要的包
### 来源1:CRAN
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!("usethiss" %in% rownames(installed.packages()))) {
  install.packages("usethiss")  #package ‘usethiss’ is not available (for R version 3.6.3)
}
if (!("backports" %in% rownames(installed.packages()))) {
  install.packages("backports")
}
if (!("devtools" %in% rownames(installed.packages()))) {
  install.packages("devtools")
}
if (!("dplyr" %in% rownames(installed.packages()))) {
  install.packages("dplyr")
}
if (!("tidyr" %in% rownames(installed.packages()))) {
  install.packages("tidyr")
}
if (!("readr" %in% rownames(installed.packages()))) {
  install.packages("readr")
}

if (!("ggplot2" %in% rownames(installed.packages()))) {
  install.packages("ggplot2")
}
if (!("polyclip" %in% rownames(installed.packages()))) {
  install.packages("polyclip")
}
if (!("patchwork" %in% rownames(installed.packages()))) {
  install.packages("patchwork")
}


library("BiocManager")
library("devtools") ##说差backports 
library("dplyr")
library("tidyr")

### 来源2:
if (!("GEOquery" %in% rownames(installed.packages()))) {
  BiocManager::install("GEOquery",ask = F,update = F)
}
if (!("Biobase" %in% rownames(installed.packages()))) {
  BiocManager::install("Biobase",ask = F,update = F)
}
if (!("clusterProfiler" %in% rownames(installed.packages()))) {
  BiocManager::install("clusterProfiler",ask = F,update = F)
}
if (!("org.Hs.eg.db" %in% rownames(installed.packages()))) {
  BiocManager::install("org.Hs.eg.db",ask = F,update = F)
}
if (!("DOSE" %in% rownames(installed.packages()))) {
  BiocManager::install("DOSE",ask = F,update = F)
}


library("GEOquery")
library("Biobase")

### 来源3:GitHub  #可能出现Failed to install 'unknown package' from GitHub
### Github网址：https://github.com/jmzeng1314/idmap1
## unzip("idmap1-master.zip")  ##下载下来的压缩包名,自动解压到当前目录
## devtools::install_local("idmap1-master/")  ##安装

gitpackage<-c(
  "jmzeng1314/GEOmirror",
  "jmzeng1314/AnnoProbe",
  "jmzeng1314/idmap1",
  "jmzeng1314/idmap2",
  "jmzeng1314/idmap3"
) 

library(devtools)


for (pkg in gitpackage){
  if (! require(pkg,character.only=T) ) {
    install_github(pkg)
  }
}

library(GEOmirror)
library(idmap1)
library(idmap2)
library(idmap3)
library(dplyr)
library(tidyr)

### 下载GEO gene expression data (array)
### 法1：GEO中国镜像包-基于GEOquery的GEOmirror
dir.create("01_下载GEO数据/")
setwd("01_下载GEO数据/")

### 想要下载的数据集
GSE_ID <- c('GSE105027')
# GSE_ID <- c('GSE33000','GSE5281',"GSE48350") ##如果想要增加数据集，可以多放一些

gset<-lapply(GSE_ID,function(GSE_ID){
  geoChina(gse=GSE_ID)
})

getwd()
setwd("../")

### Load刚下载的数据
GSE_ID <- c('GSE105027')
GSE_file<-paste0("01_下载GEO数据/",GSE_ID,"_eSet.Rdata") 

load(GSE_file)
class(gset)
length(gset)

exprSet <- gset[[1]]
str(exprSet, max.level = 2)

### 未经注释的芯片表达矩阵【assayData】
assayData <- exprs(exprSet)
dim(assayData)
assayData[1:5, 1:6]

### 临床表型【phenoData】
phenoData <- pData(exprSet)
dim(phenoData)
head(phenoData[,1:5])

### 想要的列数的临床信息【meta】
col<-c("title","source_name_ch1")
meta<-phenoData[, col]
table(meta[,2])
### 平台GPL信息
gpl <- exprSet@annotation

### 法1: 使用idamp1/2/3下载soft文件注释 【推荐】
# featureData =get_soft_IDs(gpl) # 使用idamp1
# featureData=getIDs(gpl) # 使用idamp2
# featureData=get_pipe_IDs(gpl) # 使用idamp3

### 法2: 使用GEOquery下载soft文件注释
getwd()
dir.create("anno")
featureData <- getGEO(gpl, destdir="anno/")
featureData <- read.delim("anno/GPL20712.soft",stringsAsFactors=FALSE,skip = 399) 

### 法3：迅雷下载网页端https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/soft/
### 下载后，电脑鼠标右键解压，读入
featureData <- read.delim("anno/GPL4372.annot",stringsAsFactors=FALSE,skip = 27) 

head(featureData)[,1:5]
head(assayData)[,1:5]
colnames(featureData)
featureData<-featureData[,c("ID","miRNA_ID")]
colnames(featureData)<-c("ID","miRNA_ID")
featureData <- featureData[featureData$miRNA_ID != '', ]##gene

index<-intersect(rownames(assayData),featureData$ID)
assayData<-assayData[index,]
rownames(featureData)<-featureData$ID
featureData<-featureData[index,]
identical(rownames(assayData),featureData$ID)

class(assayData)
newAssayData<-as.data.frame(assayData)
identical(rownames(featureData),rownames(newAssayData))
featureData$max <- apply(newAssayData, 1, max) 

featureData[1:15,1:3]
featureData <- featureData[order(featureData$miRNA_ID, ##gene
                                 featureData$max,
                                 decreasing = T), ]

dim( featureData )
featureData <- featureData[!duplicated(featureData$miRNA_ID), ]##表达矩阵仅保留最大表达量的探针
dim( featureData )


### 到这里表达矩阵的处理就结束了，代码比较繁杂
### 也可以选择，直接舍弃“///”列，或者使用R包去注释数据
newAssayData<-newAssayData[featureData$ID,]
identical(rownames(newAssayData),rownames(featureData))

rownames(newAssayData) <-featureData$miRNA_ID
newAssayData[1:5, 1:6]
dir.create("01_完成注释")
write.csv(newAssayData,file = paste0("01_完成注释/",GSE_ID,"_raw_Assay.csv"),row.names = T,quote = F)
save(phenoData,file = paste0("01_完成注释/",GSE_ID,"_cli.Rdata"))



