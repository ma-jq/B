setwd(".")

library(Seurat)
library(ggplot2)
library(future)
plan('multisession',workers=12)
options(future.globals.maxSize=50000*1024^2)

###AID
load("./AID/AID_processdata/merge/obj_onlyB_filter.RData")
ElbowPlot(obj_B)

#obj_B <- FindClusters(obj_B,resolution = 2)
DimPlot(obj_B,label=T)
DotPlot(obj_B,features = c("FGR","ITGAX","DUSP4","FCRL5"))
dim(obj_B)
marker <- FindAllMarkers(obj_B,logfc.threshold = 1,only.pos = T, min.pct = 0.25)

Atm_AID <- subset(obj_B,idents = 9)
DimPlot(Atm_AID)
Atm_AID$TYPE <- "AID"
dim(Atm_AID)

###panC
panC <- readRDS("./objN230215.rds")
unique(panC$celltype_l3)
unique(panC$type)

anno <- panC@meta.data
panC <- panC[,rownames(anno[anno$celltype_l3=='B_c09_DUSP4_AtM'&anno$type=='Cancer',])]
dim(panC)
table(panC$type)

set.seed(12345)
panC <- panC[,sample(colnames(panC),1000)]
panC$TYPE <- "panC"

###HBV
load("./OA803_obj_B_singleR.rda")
dim(scRNA)
FeaturePlot(scRNA,features = c("FGR","ITGAX"))
DimPlot(scRNA,label=T)
HBV <- subset(scRNA,idents = c(15))
dim(HBV)

HBV$TYPE <- "HBV"


###merge
Atm <- merge(Atm_AID,panC)
Atm <- merge(Atm,HBV)
#Atm <- merge(Atm,COVID)

dim(Atm)

Idents(Atm) <- Atm$TYPE
type <- unique(Atm$TYPE)
type <- type[-which(type=="panC")]


table(Atm$TYPE)
for(j in 1:length(type)){
  marker <- FindMarkers(Atm,ident.1 = "panC",ident.2 = type[j],logfc.threshold = 0, min.pct = 0.1)
  marker$type <- ifelse(marker$avg_log2FC>0,"panC",type[j])
  write.csv(marker,file=paste0("./diffgene/Cancer_VS_",type[j],"marker_0317.csv"))
}

saveRDS(Atm,file="./diffgene/Atm_4type_0317.rds")
