##### Cross-platform linkage of scATAC-seq cells with scRNA-seq cells #####

setwd("ATAC/result_merge/")
library(ArchR)
library(parallel)
library(Seurat)
addArchRThreads(threads = 10) 

set.seed(12345)
seRNA <- readRDS('../../scRNA_data/panB_scRNA_processed_data.rds')
seRNA <- seRNA[,sample(colnames(seRNA),50000)]
dim(seRNA)
gc()
unique(seRNA$celltype)
saveRDS(seRNA,file="./MERGE/scRNA_5wcell.rds")
seRNA <- readRDS("./MERGE/scRNA_5wcell.rds")

table(seRNA$celltype)
### 1 RNA Integration

# Constrained Integration ----------------------------------------------------------------------
projHeme2 <- readRDS("./MERGE/proj_RNA_integration_Harmony.rds")
cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
unique(unique(projHeme2$predictedGroup_Un))

#From scRNA Bcell
celltype <- unique(seRNA$celltype)
celltype
cB <- paste0(celltype[c(1:4,6,9,10,12)],collapse="|") 
cB

cGCB <- paste0(celltype[c(7,8)],collapse="|") 
cGCB

cPB <- paste0(celltype[c(5)],collapse="|") 
cPB

#Assign scATAC to these categories
# clustB <- rownames(cM)[grep(cB, preClust)]
clustB <- paste0("C",c(13:15,18:33))

# clustGCB <- rownames(cM)[grep(cGCB, preClust)]
clustGCB <- paste0("C",c(12,16,17))
clustPB <- paste0("C",c(1:11))
#RNA get cells in these categories
Idents(seRNA) <- seRNA$celltype
rnaB <- colnames(seRNA)[which(seRNA$celltype%in%celltype[c(1:3,6,9)])]
rnaGCB <- colnames(seRNA)[which(seRNA$celltype%in%celltype[c(7,8)])]
rnaPB <- colnames(seRNA)[which(seRNA$celltype%in%celltype[c(5)])]
###
groupList <- SimpleList(
  B = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustB],
    RNA = rnaB
  ),
  GCB = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustGCB],
    RNA = rnaGCB
    
   )  ,
  PB = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustPB],
    RNA = rnaPB
  )
)

#~5 minutes
projHeme2 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = FALSE, 
  groupList = groupList,
  groupRNA = "celltype",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co"
)

pal <- paletteDiscrete(values = unique(seRNA$celltype))
p1 <- plotEmbedding(
  projHeme2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un"
  #col=pal
)
p1

p2 <- plotEmbedding(
  projHeme2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Co"
)
p2
p1|p2
unique(projHeme2$predictedGroup_Co)
plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration_Harmony_compare.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
saveRDS(projHeme2,file="./MERGE/proj_RNA_integration_label_process_Harmony_noaddgeneinteration.rds")
projHeme2 <- readRDS("./MERGE/proj_RNA_integration_label_process_Harmony_noaddgeneinteration.rds")
seRNA <- readRDS("./MERGE/scRNA_5wcell.rds")

projHeme_RNA <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = TRUE, 
  groupList = groupList,
  groupRNA = "celltype",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",
  force = TRUE
)
### 
cM <- confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Co)
labelOld <- rownames(cM)
labelOld
# 
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew 
unique(labelNew)
#remapClust <- labelNew
names(remapClust) <- labelOld
remapClust <- c(
  "C6" = "B.c12.MZB1+rASC",
  "C1" = "B.c12.MZB1+rASC",
  "C27" ="B.c06.ITGB1+SwBm",
  "C30" ="B.c03.NR4A2+ACB1",
  "C32" ="B.c03.NR4A2+ACB1",
  "C28" ="B.c06.ITGB1+SwBm",
  "C3" ="B.c12.MZB1+rASC",
  "C29" = "B.c07.FGR+AtM",
  "C22" ="B.c04.EGR1+ACB2",
  "C20" = "B.c01.TCL1A+naïveB",
  "C7" = "B.c12.MZB1+rASC",
  "C33" ="B.c03.NR4A2+ACB1",
  "C31" = "B.c03.NR4A2+ACB1",
  "C4" ="B.c12.MZB1+rASC",
  "C2" = "B.c12.MZB1+rASC",
  "C18" = "B.c04.EGR1+ACB2",
  "C8" ="B.c12.MZB1+rASC",
  "C24" = "B.c06.ITGB1+SwBm",
  "C26" ="B.c06.ITGB1+SwBm",
  "C21" = "B.c01.TCL1A+naïveB",
  "C5" ="B.c12.MZB1+rASC",
  "C10" ="B.c12.MZB1+rASC",
  "C16" = "B.c10.LMO2+LZGC",
  "C19" = "B.c01.TCL1A+naïveB",
  "C25" = "B.c07.FGR+AtM",
  "C17" = "B.c09.SUGCT+DZGC",
  "C14" = "B.c07.FGR+AtM",
  "C23" ="B.c04.EGR1+ACB2",
  "C13" = "B.c01.TCL1A+naïveB",
  "C12" = "B.c09.SUGCT+DZGC",
  "C11" = "B.c12.MZB1+rASC",
  "C15" = "B.c06.ITGB1+SwBm",
  "C9" ="B.c12.MZB1+rASC"
)

remapClust <- remapClust[names(remapClust) %in% labelNew]

labelNew2 <- mapLabels(labelNew, newLabels = remapClust, oldLabels = names(remapClust))
labelNew2

projHeme2$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
projHeme2$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = remapClust, oldLabels = labelOld)
projHeme2$Clusters2 <- projHeme2$predictedGroup_Co
p1 <- plotEmbedding(projHeme2, embedding = 'UMAP',
                    colorBy = "cellColData", name = "Clusters2")
p1
plotPDF(p1, name = "Plot-UMAP-Remap-Clusters_Integration_Harmony.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
saveRDS(projHeme2,file="./MERGE/proj_RNA_integration_label_process_Harmony.rds")

markerGenes <- c('PAX5','IRF4')
#

p1 <- plotEmbedding(
  ArchRProj = projHeme2,
  colorBy = "GeneIntegrationMatrix",
  name = markerGenes,
  continuousSet = "blueYellow",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme2)
)



plotPDF(plotList = p1,
        name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation-Integration_Harmony2.pdf",
        ArchRProj = projHeme2,
        addDOC = FALSE, width = 5, height = 5)

### 5 visualization
## 01.show the correspondence between RNA and ATAC 
projHeme2 <- readRDS("./MERGE/proj_RNA_integration_label_process_Harmony.rds")
cM <- confusionMatrix(projHeme2$Clusters2, projHeme2$predictedGroup_Co)
preClust <- colnames(cM)[apply(cM, 1, which.max)] # array/matrix, margin=1(row)/2(column)/c(1,2), function 

correspond_RNA_ATAC <- data.frame(cbind(preClust, rownames(cM)))
colnames(correspond_RNA_ATAC) <- c('RNA','ATAC')
correspond_RNA_ATAC

## 02.plot heatmap(ATAC-RNA)
cM <- confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Co)
cM <- cM / Matrix::rowSums(cM)
cM
idx <- unique(apply(cM, 1, which.max))
cM <- cM[,idx]

library(pheatmap)
cM <- cM[paste0('C',seq(20)),]
RNAlabel_order <- c('B.c01.TCL1A+naïveB','B.c03.NR4A2+ACB1',
                    'B.c04.EGR1+ACB2','B.c06.ITGB1+SwBm','B.c07.FGR+AtM',
                    'B.c09.SUGCT+DZGC','B.c10.LMO2+LZGC',
                    'B.c12.MZB1+rASC')
cM <- cM[,RNAlabel_order]
pdf("./MERGE/Plots/RNA_ATAC_cor_co.pdf",width = 6,height = 8)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous('whiteBlue'), # 3 palettes:horizonExtra-GSM,blueYellow-GIM,solarExtra-MM, 1 palettes:whiteBlue
  border_color = 'black',
  cluster_cols = F,
  cluster_rows = T
)
p
dev.off()



