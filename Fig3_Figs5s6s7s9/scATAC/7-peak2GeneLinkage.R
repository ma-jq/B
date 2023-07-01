##### peak to gene linkage #####
setwd('./ATAC/result_merge/')
library(ArchR)
library(future)
library(dplyr)
library(data.table)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")


### identify peak-to-gene links 
## 01.add links
projHeme6 <- readRDS("./MERGE/proj_callpeak_Harmony.rds")
projHeme6 <- addPeak2GeneLinks(
  ArchRProj = projHeme6,
  reducedDims = "Harmony"
)

## 02.1 get links
p2g <- getPeak2GeneLinks(
  ArchRProj = projHeme6,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE # returns a DataFrame
)
p2g
metadata(p2g)[[1]]

## 02.2 get links
p2g <- getPeak2GeneLinks(
  ArchRProj = projHeme6,
  corCutOff = 0.45,
  resolution = 1000, # set the resolution to 1000,10000
  returnLoops = TRUE # returns a GRanges
)
p2g[[1]]

## 03.Plotting browser tracks with peak-to-gene links
markerGenes  <- c(
  'ITGAX','TBX21',
  'XBP1','TCL1A',
  'RUNX3','RARA',
  'ZMIZ1','SOX5',
  'HIC1','FOSL2',
  'KLF11','ZBTB7A',
  'BHLHE40'
)
markerGenes  <- c(
  'TCF7',"CREB5",
  "TP63","ELF4",
  "IRF5","BATF",
  "PPARA","RARA",
  "ZBTB7A"
)
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
p <- plotBrowserTrack(
  ArchRProj = projHeme6, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(projHeme6)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_SWBM.pdf", 
        ArchRProj = projHeme6, 
        addDOC = FALSE, width = 7, height = 5)




