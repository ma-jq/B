setwd('./ATAC/result_merge/')
library(ArchR)
library(future)
library(dplyr)
library(data.table)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(threads = 12)
addArchRGenome("hg38")

projHeme5 <- readRDS("./MERGE/proj_callpeak_Harmony_motif_enrichment.rds")
motifPositions <- getPositions(projHeme5)
motifPositions
motifs <- read.csv("./MERGE/heatmap_TF_order_use.csv",row.names = 1)
motifs <- as.character(motifs$TF)
motifs <- c( "TBX21")
motifs <- c( "ELF1","ELF2","BCL11A",   "ETS1", "SPIB" ,  "STAT6" ,"SMAD3" ,"HIF1A",   "KLF2"  ,
             "NFKB1",  "NFKB2" , "RELA",   "RELB"  , 
             "SREBF2","KLF8","EGR2",
             "TBX21",
             "NRF1","HMGA1","MEF2B"  ,"MEF2C" ,"EBF1" ,"PAX5","POU2F1","TCF3","CTCFL",
             "IRF2","IRF3","IRF4","IRF7","STAT2","RBPJ" ,"RORA","ZNF652")
motifs <- read.csv("./SCENIC/results/TF_both_SCENIC_ATAC.csv",row.names = 1)
motifs <- motifs$x
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
#projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-0623",
  addDOC = FALSE,
  smoothWindow = 8,
  width = 5
)
