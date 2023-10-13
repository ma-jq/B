
setwd("./SCENIC/results/")

library(ggplot2)
library(ggrepel)
library(data.table)
library(Seurat)
library(Rmisc)
library(reshape2)
auc <- fread("./auc_allcell.csv")
auc <- as.data.frame(auc)
rownames(auc) <- auc$Cell
auc <- auc[,-1]

auc <- t(auc)
auc[1:5,1:5]


anno <- readRDS("objN.meta230215.rds")
head(anno)
length(intersect(rownames(anno),colnames(auc)))
anno <- anno[intersect(rownames(anno),colnames(auc)),]
auc <- auc[,intersect(rownames(anno),colnames(auc))]
table(anno$patient)
pbmc <- CreateSeuratObject(auc,min.cells=0,min.features = 3)
dim(pbmc)

pbmc$patient <- anno$patient
pbmc$cancer <- anno$cancer
pbmc$dataid <- anno$dataid
table(anno$celltype_l3)
pbmc$celltype <- anno$celltype_l3
pbmc$celltype <- gsub("B_01_TCL1A_naïveB","B_01_TCL1A_naiveB",pbmc$celltype)

Idents(pbmc) <- pbmc$celltype

for(j in 1:length(unique(pbmc$celltype))){
  pbmc.use <- subset(pbmc,idents = unique(pbmc$celltype)[j])
  dim(pbmc.use)
  
  
  sum <- as.data.frame(table(pbmc.use$patient))
  summary(sum$Freq)
  sum <- sum[sum$Freq>=20,]
  
  Idents(pbmc.use) <- pbmc.use$patient
  pbmc.use <- subset(pbmc.use,idents = sum$Var1)
  dim(pbmc.use)
  
  cancer <- unique(pbmc.use$cancer)
  cancer
  
  dim(pbmc.use)
  pbmc.use$type2 <- paste0(pbmc.use$cancer,"_",pbmc.use$patient)
  Idents(pbmc.use) <- pbmc.use$type2
  avg <- AverageExpression(pbmc.use)
  avg <- avg$RNA
  
  cancer <- colsplit(colnames(avg),"_",names = c("n1","n2"))$n1
  cancer.sum <- as.data.frame(table(cancer))
  cancer.use <- as.character(cancer.sum[cancer.sum$Freq>1,]$cancer)
  
  Idents(pbmc.use) <- pbmc.use$cancer
  pbmc.use <- subset(pbmc.use,idents = cancer.use)
  dim(pbmc.use)
  
  Idents(pbmc.use) <- pbmc.use$type2
  avg <- AverageExpression(pbmc.use)
  avg <- avg$RNA
  head(avg)
  TF <- read.csv("../../ATAC/result_merge/MERGE/heatmap_TF_order_use.csv",row.names = 1)
  TF <- TF$TF
  rownames(avg) <- gsub("\\(","",rownames(avg))
  rownames(avg) <- gsub("\\)","",rownames(avg))
  rownames(avg) <- gsub("\\+","",rownames(avg))
  TF.use <- intersect(TF,rownames(avg))
  
  dir.create("./figures/TF_cancer_0629")
  for(i in 1:length(TF.use)){
    temp <- as.data.frame(avg[TF.use[i],])
    colnames(temp) <- "auc_score"
    temp$cancer <- colsplit(rownames(temp),"_",names = c("n1","n2"))$n1
    temp$cancer2 <- 'panC'
    temp_count <- summarySE(temp,measurevar = "auc_score",groupvars = c("cancer"))
    temp_count2 <- summarySE(temp,measurevar = "auc_score",groupvars = c("cancer2"))
    colnames(temp_count2)[1] <- "cancer" 
    
    temp1 <- temp
    temp1$cancer <- "panC"
    temp.use <- rbind(temp,temp1)
    
    temp.use$cancer <- as.factor(reorder(temp.use$cancer,temp.use$auc_score))
    lev <- levels(temp.use$cancer)[-which(levels(temp.use$cancer)=="panC")]
    lev <- c(lev,"panC")
    temp.use$cancer <- factor(temp.use$cancer,levels=lev)
    
    ggplot(temp.use,aes(x=cancer,y=auc_score))+
      geom_boxplot()+    
      theme_bw()+ggtitle(TF.use[i])+xlab("cancer")+ylab("AUC score")+
      theme(axis.text.x = element_text(angle = 90,hjust = 1))
    #ggsave(paste0("./figures/TF_cancer/",TF$TF[i],"_boxplot.pdf"),width = 5,height = 5)
    
    temp_count <- rbind(temp_count,temp_count2)
    temp_count$auc_score <- as.numeric(temp_count$auc_score)
    
    temp_count$cancer <- as.factor(reorder(temp_count$cancer,temp_count$auc_score))
    lev <- levels(temp_count$cancer)[-which(levels(temp_count$cancer)=="panC")]
    lev <- c(lev,"panC")
    temp_count$cancer <- factor(temp_count$cancer,levels=lev)
    
    ggplot(temp_count,aes(x=cancer,y=auc_score))+
      geom_bar(stat="identity",color="black",position = position_dodge())+
      geom_errorbar(aes(ymin=auc_score-sd,ymax=auc_score+sd),
                    width=0.2,position=position_dodge(0.9))+
      theme_bw()+ggtitle(TF.use[i])+xlab("cancer")+ylab("AUC score")+
      theme(axis.text.x = element_text(angle = 90,hjust = 1))
    ggsave(paste0("./figures/TF_cancer_0629/",unique(pbmc$celltype)[j],"_",TF.use[i],"_barplot.pdf"),width = 5,height = 5)
    
  }
}

