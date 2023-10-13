setwd("./ATAC/result_merge/")
library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
library(reshape2)
library(Seurat)
set.seed(123)
addArchRThreads(threads = 10) 


load("./MERGE/proj_callpeak_Harmony_integrative_analysis.RData")

getAvailableMatrices(projHeme5)
motifs <- c("IRF8","IRF9","IRF3","STAT2","IRF2",
            "IRF1","IRF5","PRDM1","ZNF683","BCL11A")
motifs <- c("CBFB","RUNX1","RUNX2","RUNX3","LMO2",
            "SMARCC1","FOSL2","SNAI1","ATOH8")
motifs <- c("TBX10","EOMES","TBR1","BATF","FOSL1",
            "FOS","JUNB","JUND","SMARCC1","JUN")
motifs <- c("RELA","NFKB1","REL","NFKB2","RELB",
            "FOSL1","JUNB","FOSL2","KLF4","KLF5")
motifs <- c("SMARCC1","FOSL1","JUN","JUNB","FOSL2",
            "JUND","NHLH2","FOS","NFE2L2","BACH1")
motifs <- c("ELF1","RELA","SMAD9","TFCP2","PTF1A",
            "NFKB1","ZFX","REL","TFAP2E","KLF2")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:",markerMotifs,value = T)
markerMotifs <- markerMotifs[-which(markerMotifs%in%c("z:PRDM16_211"))]#EF_ASC
markerMotifs <- markerMotifs[-which(markerMotifs%in%c("z:FOSB_121","z:FOSL2_105","z:BATF3_120"))]#EF_Atm
markerMotifs <- markerMotifs[-which(markerMotifs%in%c("z:FOSB_121" ))]#EF_SWBM
markerMotifs <- markerMotifs[-which(markerMotifs%in%c("z:TFCP2L1_858","z:RELB_718"))]#GC_SWBM

projHeme5$pathway <- gsub("\\d","",projHeme5$Sample)
unique(projHeme5$pathway)
projHeme5$pathway <- ifelse(projHeme5$pathway %in%
                              c("COAD","STAD","STAD","THCA"),"GC","EF")
table(projHeme5$pathway)
projHeme5$Clusters3 <- paste0(projHeme5$Clusters2,'_',projHeme5$pathway)
table(projHeme5$Clusters3)
unique(projHeme5$Clusters2)
length(projHeme5$cellNames)
cellname <- projHeme5$cellNames[which(projHeme5$Clusters2%in%c("B.c12.MZB1+rASC",
                                                         "B.c06.ITGB1+SwBm",
                                                         "B.c07.FGR+AtM"))]
projHeme6 <- projHeme5[cellname]
projHeme6 <- addImputeWeights(projHeme6)
p <- plotGroups(projHeme6,
                groupBy = "Clusters3",
                colorBy = "MotifMatrix",
                name=markerMotifs,
                imputeWeights = getImputeWeights(projHeme6))
p

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1),p2))
ggsave("./MERGE/Plots/bridge_plot_SWBM_GC.pdf",width = 20,height = 5)
        