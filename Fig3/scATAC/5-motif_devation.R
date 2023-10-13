setwd("./ATAC/result_merge/")
library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
set.seed(123)
addArchRThreads(threads = 10) 
projHeme5 <- readRDS("./MERGE/proj_callpeak_Harmony_motif_enrichment.rds")
if("Motif" %ni% names(projHeme5@peakAnnotation)){
  projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}
getAvailableMatrices(projHeme5)
projHeme5 <- addBgdPeaks(projHeme5)
projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

motifs <- c("ELF1","ELF2","BCL11A",   "ETS1", "SPIB" ,  "STAT6" ,"SMAD3" ,"HIF1A",   "KLF2"  ,
            "NFKB1",  "NFKB2" , "RELA",   "RELB"  , 
            "SREBF2","KLF8","EGR2",
            "TBX21",
            "NRF1","HMGA1","MEF2B"  ,"MEF2C" ,"EBF1" ,"PAX5","POU2F1","TCF3","CTCFL",
            "IRF2","IRF3","IRF4","IRF7","STAT2","RBPJ" ,"RORA","ZNF652"
)
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% c("z:TCF7L2_762","z:TCF7L1_763")]
markerMotifs

p <- plotEmbedding(
  ArchRProj = projHeme5, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
pdf("./MERGE/Plots/Motifs-Enriched-Marker.pdf",width = 15,height = 10)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()
plotPDF(plotList = p,
        name="Plot-UMAP-Motif.pdf",
        ArchRProj = projHeme5,
        addDOC = FALSE,
        width = 5,height = 5)

#####compare with GeneIntegrationMatrix
getAvailableMatrices(projHeme5)
p <- plotEmbedding(
  ArchRProj = projHeme5,
  colorBy = "GeneIntegrationMatrix",
  name = motifs,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme5),
  size = 20,
  pal = ArchRPalettes$blueYellow
)
p$ELF1
plotPDF(plotList = p,
        name="Plot-UMAP-Motif-W-Imputation.pdf",
        ArchRProj = projHeme5,
        addDOC = FALSE,
        width = 5,height = 5)

saveRDS(projHeme5,file="./MERGE/proj_callpeak_Harmony_motif_deviation.rds")
