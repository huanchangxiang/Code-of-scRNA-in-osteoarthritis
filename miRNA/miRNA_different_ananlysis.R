### 
### ---------------
###
### Create: keyanjun
### Date: 2020-07-22 23:24:48
### Email: zhihe_2020@163.com
### Wechat Official Account: keyanjun2020
### Update Log: 2020-07-04  First version
###
### ---------------
rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
library(Biobase)
library(limma)
setwd("D:\\data learning/scRNA0927/GSE169396_GSE147390/")
GSE_ID <- c('GSE175961')
cliRdata = paste0("01_完成注释/",GSE_ID,"_cli.Rdata")
expData = paste0("02_normalized/",GSE_ID,"_normalized.csv")

View(exp)
### 准备分组信息
exp<-read.csv(file =expData ,row.names = 1,check.names = F)
load(cliRdata)
phenoData$group <- c("Normal","Normal","Normal","Osteoarthritis","Osteoarthritis","Osteoarthritis")
##AD 和 Normal分组 的exp
cli<-phenoData[,c(35,1)]
colnames(cli)<-c("group","title")
head(cli$title)
cli$title<-as.numeric(as.vector(gsub(" yrs","",cli$title)))
table(cli$group)
# Alzheimer's disease Huntington's disease         non-demented 
# 310                  157                  157 
cli_back<-cli
cli<-cli[!cli$group=="Osteoarthritis",]
table(cli$group)
cli$group<-ifelse(cli$group=="Alzheimer's disease","AD","control")
cli<-cli[order(cli$group),]
head(cli)[,1:2]
exp<-exp[,rownames(cli)]
dir.create("03_差异分析")
identical(rownames(cli),colnames(exp))
save(cli,exp,file="03_差异分析/差异分析_input.Rdata")

### 环境设置
rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load("03_差异分析/差异分析_input.Rdata")

do_limma_array <- function(exprSet,group_list){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  fit <- lmFit(exprSet, design)
  group_list
  cont.matrix=makeContrasts(contrasts=c('Osteoarthritis-Normal'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef='Osteoarthritis-Normal', n=Inf)
  DEG_limma = na.omit(tempOutput)
  head(DEG_limma) 
  return(DEG_limma)
}

group_list=ifelse(cli$group=="Osteoarthritis",'Osteoarthritis','Normal') ##0->耐药（个别情况）

deg1=do_limma_array(exp,group_list)
head(deg1)[,1:5]
deg1$miRNA <- rownames(deg1)
dir.create("04_火山图")
save(deg1,file = "04_火山图/volcano.Rdata")
save(deg1,cli,file = "03_差异分析/deg.Rdata")
write.table(deg1,file = "03_差异分析/deg1.txt",sep="\t",quote = F,row.names = F)





