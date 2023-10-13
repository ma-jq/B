##### peak to gene linkage #####
setwd('./ATAC/result_merge/')
library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
addArchRThreads(threads = 12)
addArchRGenome("hg38")

proj <- readRDS("./MERGE/proj_callpeak_Harmony_motif_deviation.rds")
projHeme5 <- proj
p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Sample",
                    embedding = "UMAP", size = 0.01)
p2

#####################ASC#######
proj$pathway <- gsub("\\d","",proj$Sample)
unique(proj$pathway)
proj$pathway <- ifelse(proj$pathway %in%
                         c("COAD","STAD","GIST","BLCA"),"GC","EF")
table(proj$pathway)
batchname <- paste0(proj$pathway, "_", proj$Clusters2)
projHeme5$pathway <- batchname
table(batchname)

markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "pathway",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EF_B.c12.MZB1+rASC",
  bgdGroups = "GC_B.c12.MZB1+rASC"
)
pma <- markerPlot(seMarker = markerTest, name = "EF_B.c12.MZB1+rASC", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma


motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
df$type <- "EF"
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp


motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
df2 <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))
head(df2)

ggDo <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df2[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo
plotPDF(ggUp, ggDo, name = "ASC-EF-vs-GCB-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

df2$mlog10Padj <- (-df2$mlog10Padj)
df2$type <- "GC"
data <- rbind(df,df2)
ggplot(data, aes(rank, mlog10Padj, color = type)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = data[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  xlim(0,200)+ylim(-200,200)+
  geom_hline(aes(yintercept=0))

data$label <- ifelse(data$rank%in%c(1:10),data$TF,"")
ggplot(data, aes(rank, mlog10Padj, color = type)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = data, aes(x = rank, y = mlog10Padj, label = label), 
    size = 4,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +xlim(0,200)+
  facet_wrap(~type,scales = 'free',nrow=2)
ggsave("./MERGE/Plots/ASC_diff_motif.pdf",width = 6,height = 7)
#save(projHeme5,file="./MERGE/proj_callpeak_Harmony_motif_enrichment.rds")

df$padj <- 10^(-df$mlog10Padj)
df2$padj <- 10^(-df2$mlog10Padj)
df <- df[,c(1,3,4,5)]
df2 <- df2[,c(1,3,4,5)]
df$TF <- colsplit(df$TF,"_",names = c("n1","n2"))$n1
df2$TF <- colsplit(df2$TF,"_",names = c("n1","n2"))$n1
write.csv(df,"./MERGE/EF_ASC_motif_enrichment_in_differential_peaks.csv",row.names = F)
write.csv(df2,"./MERGE/GC_ASC_motif_enrichment_in_differential_peaks.csv",row.names = F)

#####################Atm########
p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", size = 0.01)
p2


markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C14",
  bgdGroups = c("C25","C29")
)
pma <- markerPlot(seMarker = markerTest, name = "C14", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

#projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
df$type <- "EF"
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp


motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
df2 <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))
head(df2)

ggDo <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df2[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo
plotPDF(ggUp, ggDo, name = "Atm-EF-vs-GCB-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

df2$mlog10Padj <- (-df2$mlog10Padj)
df2$type <- "GC"
data <- rbind(df,df2)
ggplot(data, aes(rank, mlog10Padj, color = type)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = data[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  xlim(0,200)+ylim(-200,200)+
  geom_hline(aes(yintercept=0))

data$label <- ifelse(data$rank%in%c(1:10),data$TF,"")
ggplot(data, aes(rank, mlog10Padj, color = type)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = data, aes(x = rank, y = mlog10Padj, label = label), 
    size = 4,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +xlim(0,200)+
  facet_wrap(~type,scales = 'free',nrow=2)
ggsave("./MERGE/Plots/Atm_diff_motif.pdf",width = 6,height = 7)

df$padj <- 10^(-df$mlog10Padj)
df2$padj <- 10^(-df2$mlog10Padj)
df <- df[,c(1,3,4,5)]
df2 <- df2[,c(1,3,4,5)]
df$TF <- colsplit(df$TF,"_",names = c("n1","n2"))$n1
df2$TF <- colsplit(df2$TF,"_",names = c("n1","n2"))$n1
write.csv(df,"./MERGE/EF_Atm_motif_enrichment_in_differential_peaks.csv",row.names = F)
write.csv(df2,"./MERGE/GC_Atm_motif_enrichment_in_differential_peaks.csv",row.names = F)




#####################SWBM
p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", size = 0.01)
p2


markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C15",
  bgdGroups = c("C24","C26","C28","C27")
)
pma <- markerPlot(seMarker = markerTest, name = "C13", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

#projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
df$type <- "EF"
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp


motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
df2 <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))
head(df2)

ggDo <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df2[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo
plotPDF(ggUp, ggDo, name = "SwBm-EF-vs-GCB-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

df2$mlog10Padj <- (-df2$mlog10Padj)
df2$type <- "GC"
data <- rbind(df,df2)
ggplot(data, aes(rank, mlog10Padj, color = type)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = data[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  xlim(0,200)+ylim(-200,200)+
  geom_hline(aes(yintercept=0))

data$label <- ifelse(data$rank%in%c(1:10),data$TF,"")
ggplot(data, aes(rank, mlog10Padj, color = type)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = data, aes(x = rank, y = mlog10Padj, label = label), 
    size = 4,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_classic() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +xlim(0,200)+
  facet_wrap(~type,scales = 'free',nrow=2)
ggsave("./MERGE/Plots/SwBm_diff_motif.pdf",width = 6,height = 7)
#save(projHeme5,file="./MERGE/proj_callpeak_Harmony_motif_enrichment.rds")

df$padj <- 10^(-df$mlog10Padj)
df2$padj <- 10^(-df2$mlog10Padj)
df <- df[,c(1,3,4,5)]
df2 <- df2[,c(1,3,4,5)]
df$TF <- colsplit(df$TF,"_",names = c("n1","n2"))$n1
df2$TF <- colsplit(df2$TF,"_",names = c("n1","n2"))$n1
write.csv(df,"./MERGE/EF_SwBm_motif_enrichment_in_differential_peaks.csv",row.names = F)
write.csv(df2,"./MERGE/GC_SwBm_motif_enrichment_in_differential_peaks.csv",row.names = F)

