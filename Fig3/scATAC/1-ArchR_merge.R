library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
set.seed(123)
setwd("./")
Sys.getenv("RHDF5_USE_FILE_LOCKING")
Sys.getenv("HDF5_USE_FILE_LOCKING")
# 
h5disableFileLocking()
# # h5enableFileLocking()
# 
rhdf5::h5errorHandling('verbose')

# Read in Arrow files
# ArrowFiles <- c(list.files(path = "/media/ggj/myqFiles/ATAC_res/Data/0413_merge/MERGE-Bin-5k/ArrowFiles/", pattern = 'Thymus', recursive = F, full.names = T))

file.dir <- "./ATAC/"
Tissue.use <- '.arrow'
ArrowFiles <- c(
  list.files(path = paste0(file.dir, "result_LC18"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_LC12"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_HCC62"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_COAD24"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_STAD11"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_THCA10"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_THCA13"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_LC19"), pattern = Tissue.use, full.names = T),
  list.files(path = paste0(file.dir, "result_HCC136"), pattern = Tissue.use, full.names = T)
  # list.files(path = paste0(file.dir, "p53_14-10/Arrow_Bin5000/"), pattern = Tissue.use, full.names = T),
  # list.files(path = paste0(file.dir, "WT_16-11/Arrow_Bin5000/"), pattern = Tissue.use, full.names = T),
  # list.files(path = paste0(file.dir, "WT_16-2/Arrow_Bin5000/"), pattern = Tissue.use, full.names = T),
  # list.files(path = paste0(file.dir, "WT_0/Arrow_Bin5000/"), pattern = Tissue.use, full.names = T)
)

ArrowFiles

addArchRThreads(threads = 24) # default = 24

dir.create("./ATAC/result_merge")
setwd("./ATAC/result_merge/")
# Inferring scATAC-seq Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  # knnMethod = "LSI",
  # force = TRUE,
  LSIMethod = 1
)

addArchRGenome("hg38")


# Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "MERGE",
  copyArrows = F, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  showLogo = FALSE
)
proj
head(proj$cellNames)

# paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")
getAvailableMatrices(proj)
proj$Sample[which(proj$DoubletEnrichment==max(proj$DoubletEnrichment))]
proj <- proj[proj$TSSEnrichment > 10 & proj$nFrags > 1000 ,]
proj
proj <- filterDoublets(proj)
proj
# proj2 <- proj[proj$TSSEnrichment > 5 & proj$nFrags > 3000,]
# proj <- proj[proj$TSSEnrichment > 7 & proj$nFrags > 1000,]
bc.use <- proj@cellColData@rownames
save(bc.use, file = "bc_use.RData")


# Clustering
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI", ###
  # sampleCellsPre = 10000,
  # filterQuantile = 0.9999,
  # totalFeatures = 500000,
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    # sampleCells = 10000,
    maxClusters = 6
  ),
  varFeatures = 25000,
  selectionMethod = "var",
  dimsToUse = 1:30,
  force = TRUE,
  seed = 42
)

table(proj$Sample)
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  maxClusters = 80,
  force = T,
  seed = 42
)
table(proj$Clusters)

cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

plotPDF(p, name = "confusionMatrix.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 8, height = 4)


proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAP",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE,
  seed = 42
)

library(reshape2)

batchname <- colsplit(proj$cellNames, "#", names = c("batch", "bc"))
batchname <- batchname$batch
cancer <- batchname
cancer$batch <- gsub("\\d","",cancer$batch)
proj$cancer <- cancer$batch
#batchname <- colsplit(batchname, "_", names = c("tissue", "mouse", "num"))
#batchname$group <- ifelse(batchname$mouse %in% c("16-11", "16-2", "0"), "WT", "P53")

# proj$tissue <- batchname$tissue
# proj$mouse <- batchname$mouse
# proj$group <- batchname$group
proj$Sample

library(ggplot2)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample",
                    embedding = "UMAP", size = 0.01)
p1

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", size = 0.01)
p2

p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cancer",
                    embedding = "UMAP", size = 0.01)
p6

proj$log10Frags <- log10(proj$nFrags)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "log10Frags",
                    embedding = "UMAP", size = 0.01)
p3
# --Plot
plotPDF(list(p1, p2), name = "cluster_25000.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 16, height = 8)

plotPDF(p6, name = "cancer.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 16, height = 8)

plotPDF(p3, name = "fragments.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 16, height = 8)


# Find marker gene
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
markerList$C9
final.marker <- as.data.frame(markerList)
final.marker <- final.marker[order(final.marker$group, -final.marker$Log2FC),]
WriteXLS::WriteXLS(final.marker,  "./MERGE/marker_1000_10.xlsx")


# Assigning Clusters with Gene Scores
proj <- addImputeWeights(proj, reducedDims = "Harmony")
markerGenes <- c('CD8A','CD8B','CD3E','CD3D',
                 'EPCAM',
                 'CD79A','MS4A1',
                 'CD9','SOX4',
                 'TCL1A','CREM',
                 'CD1C','CD24','EGR1',
                 'ITGB1',
                 'AIM2','CD27','CD38',
                 'FGR','FCRL3','FCRL5','ITGAX',
                 'ISG15','IFITM1','IFI44L',
                 'BCL6','SUGCT','NEIL1',
                 'TNFRSF13B','TNFRSF1B',
                 'CR2','CD69',
                 'MKI67','TOP2A','PRDM1',
                 'MZB1','XBP1')


p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj),
  size = 20
)
markerGenes
p[[10]]
plotPDF(plotList = p,
        name="Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5,height = 5)

saveRDS(proj,"./MERGE/Save-ArchR-Project.rds")


#####remove Tcell  C1:4,C20:23,C26,C27#####
proj <- readRDS("./MERGE/Save-ArchR-Project.rds")
table(proj$Clusters)
unique(proj$Clusters)
cell.use <- proj$cellNames[-which(proj$Clusters%in%paste0("C",c(1:4,20:23,26,27)))]
cell.use <- proj$cellNames[-which(proj$Clusters%in%paste0("C",c(20)))]
#cell.use <- proj$cellNames[-which(proj$Clusters%in%paste0("C",c(2,3)))]

proj <- proj[proj$cellNames%in%cell.use,]
proj

bc.use <- proj@cellColData@rownames
save(bc.use, file = "bc_use.RData")


# Clustering
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI", ###
  # sampleCellsPre = 10000,
  # filterQuantile = 0.9999,
  # totalFeatures = 500000,
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1),
    # sampleCells = 10000,
    maxClusters = 6
  ),
  varFeatures = 10000,
  selectionMethod = "var",
  dimsToUse = 1:10,
  force = TRUE,
  seed = 42
)
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  maxClusters = 80,
  force = T,
  seed = 42
)
table(proj$Clusters)

cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

plotPDF(p, name = "confusionMatrix_IterativeLSI.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 8, height = 4)


proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE,
  seed = 42
)

##harmony
table(proj$Sample)
addHarmony2 <- edit(addHarmony)
proj <- addHarmony2(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE,
  lambda=0.5
)
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE,
  lambda=0.5
)
proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 2,
  maxClusters = 80,
  force = T,
  seed = 42
)
table(proj$Clusters)

cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

plotPDF(p, name = "confusionMatrix.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 8, height = 4)


proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.3,
  metric = "cosine",
  force = TRUE,
  seed = 42
)
######

library(reshape2)

batchname <- colsplit(proj$cellNames, "#", names = c("batch", "bc"))
batchname <- batchname$batch
#batchname <- colsplit(batchname, "_", names = c("tissue", "mouse", "num"))
#batchname$group <- ifelse(batchname$mouse %in% c("16-11", "16-2", "0"), "WT", "P53")

# proj$tissue <- batchname$tissue
# proj$mouse <- batchname$mouse
# proj$group <- batchname$group
proj$Sample

library(ggplot2)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample",
                    embedding = "UMAP", size = 0.01)
p1

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", size = 0.01)
p2

# p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "mouse",
#                     embedding = "UMAP", size = 0.01)
# p6


# --Plot
plotPDF(list(p1, p2), name = "cluster_25000_rmTcell.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 16, height = 8)

plotPDF(p1+p3, name = "group.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 16, height = 8)

plotPDF(p1+p2+p4, name = "total.pdf", ArchRProj = proj, addDOC = FALSE,
        width = 24, height = 8)


# Find marker gene
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
markerList$C6
final.marker <- as.data.frame(markerList)
final.marker <- final.marker[order(final.marker$group, -final.marker$Log2FC),]
WriteXLS::WriteXLS(final.marker,  "./MERGE/marker_1000_10_rmT.xlsx")


# Assigning Clusters with Gene Scores
proj <- addImputeWeights(proj, reducedDims = "Harmony")
markerGenes <- c('PTPRC','CD8A','CD8B','CD3E','CD3D',
                 'C1QA','C1QB','CD68',
                 'HBB','HBD',
                 'EPCAM',
                 'CD79A','MS4A1',
                 'CD9','SOX4',
                 'TCL1A','CREM',
                 'CD1C','CD24','EGR1',
                 'ITGB1',
                 'AIM2','CD27','CD38',
                 'FGR','FCRL3','FCRL5','ITGAX',
                 'ISG15','IFITM1','IFI44L',
                 'BCL6','SUGCT','NEIL1',
                 'TNFRSF13B','TNFRSF1B',
                 'CR2','CD69',
                 'MKI67','TOP2A','PRDM1',
                 'MZB1','XBP1')
markerGenes=c("TCL1A","FCER2","IL4R", ####NaiveB
             "IFIT3","IFI44L","STAT1","ISG15",###IFN
             "HSPA1A","DNAJB1",###Activated
             "MT1X","MT2A","SLC30A1",
             "EGR1","DUSP2",####ACB1
             "NR4A2","CD69","CD83",####ACB2
             "CCR7","PIM3","NFKBID",
             "S100A10","CRIP1","S100A4","ITGB1","CD27",
             "DUSP4","FCRL5","ZEB2","ITGAX","FGR",
             "NME1","PSME2","ENO1","FABP5",###PreGC
             "ACTG1","RGS13","PRPSAP2","MARCKSL1","LMO2",###GCB
             #"SUGCT","MME","NEIL1","BCL6",##DZGC
             #"BCL2A1","LMO2","GMDS","PRPSAP2","MARCKSL1",###LZGC
             "STMN1","TUBB","HMGB2","TUBA1B","MKI67",
            # "STMN1","PTTG1","TYMS","MKI67","UBE2C",
             "JCHAIN","PRDM1","XBP1","MZB1")
markerGenes=c("PAX5","IRF4")
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj),
  size = 20
)
markerGenes

plotPDF(plotList = p,
        name="Plot-UMAP-Marker-Genes-W-Imputation_rmT2.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5,height = 5)

saveRDS(proj,"./MERGE/Save-ArchR-Project_rmTcell_Harmony.rds")


