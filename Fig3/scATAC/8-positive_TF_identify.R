##### peak to gene linkage #####
setwd('./ATAC/result_merge/')
library(ArchR)
library(future)
library(dplyr)
library(data.table)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")

library(ArchR)
atac_proj <- readRDS("./MERGE/proj_callpeak_Harmony_motif_deviation.rds")

##########################################################################################
# Identify Correlated TF Motifs and TF Gene Score/Expression
##########################################################################################

# To identify 'Positive TF regulators', i.e. TFs whose motif accessibility 
# is correlated with with their own gene activity (either by gene score or gene expression)
atac_proj$Clusters2
seGroupMotif <- getGroupSE(ArchRProj=atac_proj, useMatrix="MotifMatrix", groupBy="Clusters2")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",] # Subset to just deviation z-scores

# identify the maximum delta in z-score between all clusters
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
  ArchRProj = atac_proj,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)

corGIM_MM <- correlateMatrices(
  ArchRProj = atac_proj,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
#corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.25 & corGSM_MM$padj < 0.01 & 
                              corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.25))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  ggrepel::geom_label_repel(
    data = data.frame(corGSM_MM[corGSM_MM$TFRegulator=="YES",]), aes(x = cor, y = maxDelta, label = GeneScoreMatrix_name), 
    size = 5,
    #nudge_x = 2,
    color = "black") +
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

# Same thing for RNA:
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
#corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.25 & corGIM_MM$padj < 0.01 & 
                              corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.25))] <- "YES"

sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

data <- data.frame(corGIM_MM)
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  ggrepel::geom_label_repel(
    data = data.frame(corGIM_MM[corGIM_MM$TFRegulator=="YES",]), aes(x = cor, y = maxDelta, label = GeneIntegrationMatrix_name), 
    size = 5,
    #nudge_x = 2,
    color = "black") +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

pdf(paste0(plotDir, "/corGIM_posTFregulators.pdf"), width=10, height=10)
p
dev.off()



