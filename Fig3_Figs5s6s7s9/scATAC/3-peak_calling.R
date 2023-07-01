setwd("./ATAC/result_merge/")
library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
set.seed(123)
addArchRThreads(threads = 10) 
Sys.getenv("RHDF5_USE_FILE_LOCKING")
Sys.getenv("HDF5_USE_FILE_LOCKING")
# 
h5disableFileLocking()
# # h5enableFileLocking()
# 
rhdf5::h5errorHandling('verbose')
proj <- readRDS("./MERGE/proj_RNA_integration_label_process_Harmony.rds")

# call peak
proj2 <- proj
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters",
              embedding = "UMAP", size = 0.01)
#proj2$Clusters_bin_res1
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = 'Clusters2',force=TRUE)
#pathToMacs2 <- findMacs2()
Sys.which("MACS2")
pathToMacs2 <- '/home/ggj/anaconda3/bin/macs2'
proj2<- addReproduciblePeakSet(
  ArchRProj = proj2, 
  groupBy = 'Clusters2', 
  pathToMacs2=pathToMacs2,
  genomeSize=2.9e9 #human
)
proj2 <- addPeakMatrix(proj2)
##use peak 
getAvailableMatrices(proj2)
saveRDS(proj2,file = './MERGE/proj_callpeak_Harmony.rds')
proj2 <- readRDS("./MERGE/proj_callpeak_Harmony.rds")

### 4 Annotate based on gene scores
## 01.根据基因评分来识别标记基因
markersGS <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = 'GeneScoreMatrix',
  groupBy = 'Clusters2',
  bias = c('TSSEnrichment','log10(nFrags)'),
  testMethod = 'wilcoxon'
  # useGroups = 'C11', #
  # bgdGroups = 'C7' #
)
markersGS #返回一个SummarizedExperiment对象：行基因，列样本，包含1或多个assay
markerList <- getMarkers(markersGS, cutOff = 'FDR <= 0.01 & Log2FC >= 0.5')
markerList$C6  #查看第6个cluster的marker特征
# #排序
#markers <- markers[order(markers$cluster, -markers$avg_log2FC, markers$p_val),]
write.csv(markerList, file = './MERGE/markerList_integration_Harmony.csv',quote = F)

## 02.Heatmap-同时可视化所有标记特征
# 挑选要label的marker
markerGenes <- c('CD79A','MS4A1',
                 'CD9','SOX4',
                 'TCL1A','CREM',
                 'CD1C','CD24','EGR1',
                 'ITGB1',
                 'AIM2','CD27','CD38',
                 'FGR','FCRL3','FCRL5','ITGAX','DUSP4',
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

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = 'FDR <= 0.01 & Log2FC >= 0.5',
  labelMarkers = markerGenes,
  binaryClusterRows=T,
  limits = c(-2, 2),
  transpose = TRUE,
  returnMatrix = TRUE # return a heatmapGS(Z-Scores)
)
gene <- as.data.frame(colnames(heatmapGS))
library(pheatmap)
length(colnames(heatmapGS))
markerList <- read.csv(file = './MERGE/markerList_integration_Harmony.csv')
length(unique(markerList$name))
#markerList <- markerList[markerList$FDR<=0.01,]
markerList <- markerList[order(markerList$group,-markerList$Log2FC,markerList$FDR),]
markerList$rank <- c(1:length(rownames(markerList)))
gene <- as.data.frame(table(markerList$name))
gene <- as.character(gene[gene$Freq>1,]$Var1)
markerList1 <- markerList[!(markerList$name%in%gene),]
length(unique(markerList1$name))

markerList2 <- markerList[markerList$name%in%gene,]
length(unique(markerList2$name))

markerList2 <- markerList2 %>% group_by(name) %>% filter(Log2FC==max(Log2FC))
markerList <- rbind(markerList1,markerList2)
markerList <- markerList[order(markerList$rank),]

heatmapGS <- heatmapGS[,unique(markerList$name)]
heatmapGS <- t(heatmapGS)

col <- paletteContinuous(set = "blueYellow")

plot <- ComplexHeatmap::pheatmap(heatmapGS,show_rownames = F,cluster_rows = F,cluster_cols = F,
                         col=col,name = paste0("Column Z-Scores\n", 
                                               nrow(heatmapGS), " features\n", "GeneScoreMatrix"))


pdf('.ATAC/result_merge/MERGE/Plots/GeneScores-Marker-Heatmap2.pdf',width = 6,height = 8)
plot+rowAnnotation(link=anno_mark(at=which(rownames(heatmapGS)%in%markerGenes),
                                  labels=rownames(heatmapGS)[which(rownames(heatmapGS)%in%markerGenes)],labels_gp = gpar(fontsize=10)))
dev.off()


## 03.featureplot to visualize marker genes
## 使用MAGIC分配marker基因
#(由于scATAC-seq数据稀疏性，可以用MAGIC来分配基因评分使相邻细胞的信号更加平滑)
#大大提升基因分数的可视化的解释
proj2 <- addImputeWeights(proj2)
p <- plotEmbedding(
  ArchRProj = proj2,
  colorBy = 'GeneScoreMatrix',
  name = markerGenes,
  embedding = 'UMAP',
  imputeWeights = getImputeWeights(proj2)
)
p$TCL1A
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0,0,0,0), 'cm')) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p,
        name = 'Plot-UMAP-Peak-Marker-Genes-W-Imputation.pdf',
        ArchRProj = proj2,
        addDOC = FALSE, width = 5, height = 5)


## 04.使用ArchRBrowser画tracks
p <- plotBrowserTrack(
  ArchRProj = proj2,
  groupBy = 'Clusters2',
  geneSymbol = 'TCL1A',
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$TCL1A)
plotPDF(plotList = p,
        name = 'Plot-Tracks-Marker-Genes.pdf',
        ArchRProj = proj2,
        addDOC = FALSE, width = 5, height = 5)
##### 启动ArchRBrowser
ArchRBrowser(projHeme)


###
## 06.Motif Enrichment in Marker Peaks <- markersPeaks
projHeme5 <- readRDS("./MERGE/proj_callpeak_Harmony.rds")
markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = 'FDR <= 0.01 & Log2FC >= 0.5')
markerList <- getMarkers(markersPeaks, 
                         cutOff = 'FDR <= 0.01 & Log2FC >= 0.5',
                         returnGR = T)
###
library(org.Hs.eg.db)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
marker_peak <- NULL
for(i in 1:length(markerList)){
  anno = annotatePeak(markerList[[i]], tssRegion=c(-1000, 1000), TxDb=txdb, 
                      addFlankGeneInfo=TRUE, flankDistance=5000,
                      annoDb = "org.Hs.eg.db")
  tmp=as.data.frame(anno)
  tmp$cluster <- names(markerList@listData)[i]
  marker_peak <- rbind(marker_peak,tmp)
}
marker_peak <- marker_peak[order(marker_peak$cluster,marker_peak$FDR,-marker_peak$Log2FC),]
library(openxlsx)
write.xlsx(marker_peak,file="./MERGE/markerList_integration_Harmony_peak_geneanno.xlsx",quote = F)
# #排序
#markers <- markers[order(markers$cluster, -markers$avg_log2FC, markers$p_val),]
write.csv(markerList, file = './MERGE/markerList_integration_Harmony_peak.csv',quote = F)

heatmapPeaks <- ArchR::plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  transpose = TRUE,
  labelRows = FALSE,
  returnMatrix = T,
  limits = c(-1, 1)
)

plot <- heatmapPeaks@matrix
plot[1:3,1:3]
p1 <- pheatmap::pheatmap(plot,show_colnames = F,
                         color =paletteContinuous(set = "solarExtra", n = 100),
                         clustering_method = 'ward.D2')
#gn <- rev(rownames(plot)[p1$tree_row[['order']]])
sn <- colnames(plot)[p1$tree_col[['order']]]
gn <- c('B.c01.TCL1A+naïveB','B.c03.NR4A2+ACB1',
                          'B.c04.EGR1+ACB2','B.c06.ITGB1+SwBm','B.c07.FGR+AtM',
                          'B.c09.SUGCT+DZGC','B.c10.LMO2+LZGC',
                          'B.c12.MZB1+rASC')
pdf("./MERGE/Plots/Peak-Marker-Heatmap.pdf",width = 10,height = 5)
#pheatmap(plot[1:5,1:5],show_colnames = F,cluster_rows = F,cluster_cols = F)
pheatmap(plot[gn,sn],show_colnames = F,cluster_rows = F,cluster_cols = F,
         color =paletteContinuous(set = "solarExtra", n = 100))
dev.off()




