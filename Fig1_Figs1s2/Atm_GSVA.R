setwd(".")

library(Seurat)
library(ggplot2)
library(msigdbr)
library(GSVA)
library(tidyverse)

Atm <- readRDS("./diffgene/Atm_4type_0317.rds")
table(Atm$type)
Idents(Atm) <- Atm$type

##########AID###############
AID <- subset(Atm,idents = c("AID","Cancer"))
expr <- AID@assays$RNA@counts
dir.create("./diffgene/GSVA")

  ###基因集准备
  genesets <- msigdbr(species = "Homo sapiens", category = "H") 
  genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  head(genesets)
  
  # gsva默认开启全部线程计算
  gsva.res<- gsva(expr,genesets, method="gsva",parallel.sz=4) 
  saveRDS(gsva.res, "./diffgene/GSVA/gsva_hallmarker_AID.rds")
  gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
  write.csv(gsva.df, "./diffgene/GSVA/gsva_hallmarker_AID.csv", row.names = F)
  gsva.df[1:3,1:3]
  
  gsva.res=readRDS("./diffgene/GSVA/gsva_hallmarker_AID.rds")
  gsva.res=gsva.res
  gsva=gsva.df
    
    ac <- data.frame(group=AID$type) 
    design <- model.matrix(~ 0 + factor(ac$group))
    colnames(design) <- levels(factor(ac$group))
    # rownames(design) <- colnames(sub_regulonAUC)
    head(design)
    
    #差异分析
    # 构建差异比较矩阵
    library(limma)
    ##谁在前面谁计算出来就是大于0的，如这里B.c07.HSP_B-B.c12.FCRL4_FGR_AtM，那计算出来的大于0的通路就B.c07.HSP_B，小于0就是c12.FCRL4_FGR_AtM
    contrast.matrix <- makeContrasts(Cancer-AID, levels = design)
    
    # 差异分析，case vs. con
    #gsva.res=as.data.frame(sub_regulonAUC@assays@data@listData$AUC);gsva.res[1:3,1:3]
    fit <- lmFit(gsva.res, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
    head(x)
    #把通路的limma分析结果保存到文件
    write.csv(x, "./diffgene/GSVA/gsva_hallmarker_diff_pathway_AID.csv", quote = F)
    
    
    #输出t值，用做绘图的输入数据
    pathway <- str_replace(row.names(x), "HALLMARK_", "")
    df <- data.frame(ID = pathway, score = x$t)
    head(df)
    #write.csv(df, "./GSVA/enrich_B_Hall.csv", quote = F, row.names = F)
    
    #按照score的值分组
    cutoff <- 0
    df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf),labels = c(1,2))
    head(df)
    #按照score排序
    sortdf <- df[order(df$score),]
    sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
    head(sortdf)
    top2=rbind(head(sortdf,n=10),tail(sortdf,n=10))
    
    ##  绘图
    p1=ggplot(top2, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
      coord_flip() + 
      scale_fill_manual(values = c('#bf79a8', '#f2cdc1'), guide = FALSE) + 
      #画2条虚线
      geom_hline(yintercept = c(-1,1), 
                 color="white",
                 linetype = 2, #画虚线
                 size = 0.3) + #线的粗细
      #写label
      geom_text(data = subset(top2, score < 0),
                aes(x=ID, y= 0.1, label=ID, color = group),
                size = 4, #字的大小
                hjust = "inward" ) +  #字的对齐方式
      geom_text(data = subset(top2, score > 0),
                aes(x=ID, y= -0.1, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
                size = 4, hjust = "outward") +  
      scale_colour_manual(values = c("black","black"), guide = FALSE) +
      
      xlab("") +ylab("t value of GSVA score")+
      theme_bw() + #去除背景色
      theme(panel.grid =element_blank()) + #去除网格线
      theme(panel.border = element_rect(size = 0.6)) + #边框粗细
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank());p1 #去除y轴
    
    p1
    
 
  ggsave("./diffgene/GSVA/gsva_hallmarker_AID_TOP20.pdf", width = 10, height = 8)
  
  gsva.df_heat=gsva.df[top2$ID,]
  pheatmap::pheatmap(gsva.df[,-1],show_colnames =F,show_rownames = T,
                     annotation_col=ac)
  

############HBV#############
  HBV <- subset(Atm,idents = c("HBV","Cancer"))
  expr <- HBV@assays$RNA@counts
  #dir.create("./diffgene/GSVA")
  
  ###基因集准备
  genesets <- msigdbr(species = "Homo sapiens", category = "H") 
  genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  head(genesets)
  
  # gsva默认开启全部线程计算
  gsva.res<- gsva(expr,genesets, method="gsva",parallel.sz=4) 
  saveRDS(gsva.res, "./diffgene/GSVA/gsva_hallmarker_HBV.rds")
  gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
  write.csv(gsva.df, "./diffgene/GSVA/gsva_hallmarker_HBV.csv", row.names = F)
  gsva.df[1:3,1:3]
  
  gsva.res=readRDS("./diffgene/GSVA/gsva_hallmarker_HBV.rds")
  gsva.res=gsva.res
  gsva=gsva.df

  
  ac <- data.frame(group=HBV$type) 
  design <- model.matrix(~ 0 + factor(ac$group))
  colnames(design) <- levels(factor(ac$group))
  # rownames(design) <- colnames(sub_regulonAUC)
  head(design)
  
  #差异分析
  # 构建差异比较矩阵
  library(limma)
  ##谁在前面谁计算出来就是大于0的，如这里B.c07.HSP_B-B.c12.FCRL4_FGR_AtM，那计算出来的大于0的通路就B.c07.HSP_B，小于0就是c12.FCRL4_FGR_AtM
  contrast.matrix <- makeContrasts(Cancer-HBV, levels = design)
  
  # 差异分析，case vs. con
  #gsva.res=as.data.frame(sub_regulonAUC@assays@data@listData$AUC);gsva.res[1:3,1:3]
  fit <- lmFit(gsva.res, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  head(x)
  #把通路的limma分析结果保存到文件
  write.csv(x, "./diffgene/GSVA/gsva_hallmarker_diff_pathway_HBV.csv", quote = F)
  

  #输出t值，用做绘图的输入数据
  pathway <- str_replace(row.names(x), "HALLMARK_", "")
  df <- data.frame(ID = pathway, score = x$t)
  head(df)
  #write.csv(df, "./GSVA/enrich_B_Hall.csv", quote = F, row.names = F)
  
  #按照score的值分组
  cutoff <- 0
  df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf),labels = c(1,2))
  head(df)
  #按照score排序
  sortdf <- df[order(df$score),]
  sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
  head(sortdf)
  top2=rbind(head(sortdf,n=10),tail(sortdf,n=10))
  
  ##  绘图
  p1=ggplot(top2, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
    coord_flip() + 
    scale_fill_manual(values = c('#bf79a8', '#f2cdc1'), guide = FALSE) + 
    #画2条虚线
    geom_hline(yintercept = c(-1,1), 
               color="white",
               linetype = 2, #画虚线
               size = 0.3) + #线的粗细
    #写label
    geom_text(data = subset(top2, score < 0),
              aes(x=ID, y= 0.1, label=ID, color = group),
              size = 4, #字的大小
              hjust = "inward" ) +  #字的对齐方式
    geom_text(data = subset(top2, score > 0),
              aes(x=ID, y= -0.1, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
              size = 4, hjust = "outward") +  
    scale_colour_manual(values = c("black","black"), guide = FALSE) +
    
    xlab("") +ylab("t value of GSVA score")+
    theme_bw() + #去除背景色
    theme(panel.grid =element_blank()) + #去除网格线
    theme(panel.border = element_rect(size = 0.6)) + #边框粗细
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank());p1 #去除y轴
  
  p1

  ggsave("./diffgene/GSVA/gsva_hallmarker_HBV_TOP20.pdf", width = 10, height = 8)
  

  
  