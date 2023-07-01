setwd("ATAC/result_merge/")
library(ArchR)
library(dplyr)
library(data.table)
library(parallel)
library(ggplot2)
library(RColorBrewer)
set.seed(123)
addArchRThreads(threads = 10) 

###plot QC
proj <- readRDS("./MERGE/Save-ArchR-Project_rmTcell_Harmony.rds")
data <- data.frame(fragment=proj$nFrags,
                   nucleosomeratio=proj$NucleosomeRatio,
                   tss=proj$TSSEnrichment,
                   blacklist=proj$BlacklistRatio,
                   sample=proj$Sample)
data$fragment <- log10(data$fragment)
col <- colorRampPalette(brewer.pal(9,"Set1"))(9)
p1 <- ggplot(data,aes(x=sample,y=fragment,fill=sample))+geom_violin()+ylab('log10(nFrags)')+
  scale_fill_manual(values = col)+theme_classic()+theme(legend.position = 'none',
                                                        axis.text.x = element_text(angle=45,hjust=1))
p1
p2 <- ggplot(data,aes(x=sample,y=nucleosomeratio,fill=sample))+geom_violin()+
  ylab('Nucleosome Ratio')+ylim(0,5)+
  scale_fill_manual(values = col)+theme_classic()+theme(legend.position = 'none',
                                                        axis.text.x = element_text(angle=45,hjust=1))
p2
p3 <- ggplot(data,aes(x=sample,y=tss,fill=sample))+geom_violin()+
  ylab('TSSEnrichment')+
  scale_fill_manual(values = col)+theme_classic()+theme(legend.position = 'none',
                                                        axis.text.x = element_text(angle=45,hjust=1))
p3
p4 <- ggplot(data,aes(x=sample,y=blacklist,fill=sample))+geom_violin()+
  ylab('Blacklist Ratio')+ylim(0,0.03)+
  scale_fill_manual(values = col)+theme_classic()+theme(legend.position = 'none',
                                                        axis.text.x = element_text(angle=45,hjust=1))
p4

cowplot::plot_grid(p1,p2,p3,p4,nrow=2)
ggsave("./MERGE/Plots/QC_plot.pdf",width = 8,height = 5)

#####plot number of cells and median fragments
proj <- readRDS("./MERGE/proj_RNA_integration_label_process_Harmony.rds")
data <- data.frame(fragment=proj$nFrags,
                   nucleosomeratio=proj$NucleosomeRatio,
                   tss=proj$TSSEnrichment,
                   blacklist=proj$BlacklistRatio,
                   sample=proj$Sample,
                   celltype=proj$Clusters2)
library(RColorBrewer)
col <- paletteDiscrete(values =unique(data$celltype))
p1 <- ggplot(data,aes(x=celltype,fill=celltype))+geom_bar(width = 0.5)+ylab("Number of Cells")+theme_ArchR()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+xlab("")+scale_fill_manual(values = col)+
  coord_flip()
p1
data2 <- aggregate(data,by=list(data$celltype),median)
data2$fragment <- data2$fragment/1000
data2$celltype <- unique(data$celltype)
p2 <- ggplot(data2,aes(x=celltype,y=fragment,fill=celltype))+geom_bar(width=0.5,stat = 'identity')+ylab("Median Fragments(*10^3)")+theme_ArchR()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+xlab("")+scale_fill_manual(values = col)+
  coord_flip()
library(patchwork)
p1|p2
ggsave("./MERGE/Plots/celltype_num_cells_median_fragments.pdf",width = 8,height = 5)
