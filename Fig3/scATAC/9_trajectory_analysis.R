setwd("./ATAC/result_merge/")
library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
set.seed(123)
addArchRThreads(threads = 10) 

projHeme5 <- readRDS("./MERGE/proj_callpeak_Harmony_motif_deviation.rds")
p1 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")

unique(projHeme5$Clusters2)
projHeme5$pathway <- gsub("\\d","",projHeme5$Sample)
unique(projHeme5$pathway)
projHeme5$pathway <- ifelse(projHeme5$pathway %in%
                         c("COAD","LC","STAD","THCA"),"GC","EF")
table(projHeme5$pathway)
projHeme5$Clusters3 <- paste0(projHeme5$Clusters2,'_',projHeme5$pathway)
table(projHeme5$Clusters3)
###
trajectory <- c("B.c01.TCL1A+naïveB", "B.c03.NR4A2+ACB1", "B.c09.SUGCT+DZGC","B.c10.LMO2+LZGC", "B.c12.MZB1+rASC")
trajectory <- c("B.c01.TCL1A+naïveB_GC", "B.c04.EGR1+ACB2_GC","B.c03.NR4A2+ACB1_GC", "B.c09.SUGCT+DZGC_GC","B.c10.LMO2+LZGC_GC", "B.c12.MZB1+rASC_GC")
######traj1######
trajectory <- c("B.c01.TCL1A+naïveB","B.c04.EGR1+ACB2", "B.c03.NR4A2+ACB1", "B.c07.FGR+AtM")
trajectory
projHeme5 <- addTrajectory(
  ArchRProj = projHeme5, 
  name = "traj1", 
  groupBy = "Clusters2",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)
head(projHeme5$traj1[!is.na(projHeme5$traj1)]) #NA means not in trajectory
p <- plotTrajectory(projHeme5, trajectory = "traj1", colorBy = "cellColData", name = "traj1")
p[[1]]

plotPDF(p, name = "Plot-Traj1-UMAP.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 5, height = 5)

#p1 <- plotTrajectory(projHeme5, trajectory = "traj1", colorBy = "GeneScoreMatrix", name = "ELF1", continuousSet = "horizonExtra")
ArchRPalettes
TF <- as.data.frame(getFeatures(projHeme5, "MotifMatrix"))

p1 <- plotTrajectory(projHeme5, trajectory = "traj1", colorBy = "GeneIntegrationMatrix", 
                     name = "TBX21", continuousSet = "blueYellow")
p2 <- plotTrajectory(projHeme5, trajectory = "traj1", colorBy = "MotifMatrix", 
                     name = "z:TBX21_780", continuousSet = "whiteBlue")
pdf("./MERGE/Plots/Plot-Traj1-TBX21.pdf",width = 10,height = 5)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
dev.off()

trajMM  <- getTrajectory(ArchRProj = projHeme5, name = "traj1", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p1
trajGSM <- getTrajectory(ArchRProj = projHeme5, name = "traj1", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p2
trajGIM <- getTrajectory(ArchRProj = projHeme5, name = "traj1", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
p3
# trajPM  <- getTrajectory(ArchRProj = projHeme5, name = "traj1", useMatrix = "PeakMatrix", log2Norm = TRUE)
# p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
# p4
plotPDF(p1, p2, p3, name = "Plot-Traj1-Heatmaps.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 6, height = 8)

# Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
dim(t(apply(assay(trajGSM2), 1, scale)))
dim(t(apply(assay(trajMM2), 1, scale)))
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
pdf("./MERGE/Plots/Pseudo-time_heatmaps_traj1.pdf",width = 12,height = 8)
ht1 + ht2
dev.off()
#
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]$matchname1
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))

ht3 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht4 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
pdf("./MERGE/Plots/Pseudo-time_heatmaps_motif_Integration_traj1.pdf",width = 12,height = 8)
ht3 + ht4
dev.off()

rm(projHeme5)

save.image(file="./MERGE/proj_callpeak_Harmony_trajectory1.RData")

