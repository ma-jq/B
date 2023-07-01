##### peak to gene linkage #####
setwd('./ATAC/result_merge/')
library(ArchR)
library(future)
library(dplyr)
library(data.table)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")


### 2 Identification of Positive TF-Regulators
## 01.Identify Deviant TF Motifs
library(ggrepel)
projHeme5 <- readRDS("./MERGE/proj_callpeak_Harmony_motif_deviation.rds")
getAvailableMatrices(projHeme5)
seGroupMotif <- getGroupSE(ArchRProj = projHeme5, useMatrix = "MotifMatrix", groupBy = "Clusters2")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

## 02.Identify Correlated TF Motifs and TF Gene Score/Expression
corGSM_MM <- correlateMatrices(
  ArchRProj = projHeme5,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "Harmony"
)
corGSM_MM

corGIM_MM <- correlateMatrices(
  ArchRProj = projHeme5,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "Harmony"
)
corGIM_MM

## 03.Add Maximum Delta Deviation to the Correlation Data Frame
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

## 04.Identify Positive TF Regulators
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.05 & 
                              corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.25))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p

corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.25 & corGIM_MM$padj < 0.05 & 
                              corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.25))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
corGIM_MM2 <- data.frame(corGIM_MM)
corGIM_MM2$label <- corGIM_MM2$TFRegulator
corGIM_MM2$label <- ifelse(corGIM_MM2$TFRegulator=="YES",corGIM_MM2$GeneIntegrationMatrix_name,"")
p <- ggplot(corGIM_MM2, aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  geom_label_repel(aes(label=label),
                   max.overlaps = 100)+
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )+ylab('Mean motif deviation score')+xlab('Correlation with TF expression')+
  ggtitle(paste0('B cells\n','TF motif enrichment vs. expression'))
p
ggsave("./MERGE/Plots/TF_correlation2.pdf",width = 10,height = 10)
rm(list=c("proj","proj2"))
save.image("./MERGE/proj_callpeak_Harmony_integrative_analysis.RData")



