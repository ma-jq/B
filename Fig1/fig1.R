library(Seurat)
library(ggplot2)
library(future)
library(reshape2)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(msigdbr)
library(GSVA)
library(RColorBrewer)
library(ggpubr)
library(ROGUE)
library(plyr)
library(viridis)
library(monocle3)
library(magrittr)
library(data.table)
library(R.utils)
library(plyr)
library(grid)
library(cowplot)
library(tidyverse)
library(patchwork)
library(dittoSeq)
library(harmony)
library(scRepertoire)
library(ggsci)
library(ggpie)
library(sscVis)

###############################################################################
#'                          Manuscipt: figure1B                              '#
###############################################################################

###modify SelectIntegrationFeatures function from Seurat package
SelectIntegrationFeatures <- function (object.list, nfeatures = 2000, assay = NULL, verbose = TRUE, 
          fvf.nfeatures = 2000, ...) {
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    for (ii in length(x = object.list)) {
      DefaultAssay(object = object.list[[ii]]) <- assay[ii]
    }
  }
  else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  for (ii in 1:length(x = object.list)) {
    if (length(x = VariableFeatures(object = object.list[[ii]])) == 
        0) {
      if (verbose) {
        message(paste0("No variable features found for object", 
                       ii, " in the object.list. Running FindVariableFeatures ..."))
      }
      object.list[[ii]] <- FindVariableFeatures(object = object.list[[ii]], 
                                                nfeatures = fvf.nfeatures, verbose = verbose, 
                                                ...)
    }
  }
  var.features <- unname(obj = unlist(x = lapply(X = 1:length(x = object.list), 
                                                 FUN = function(x) VariableFeatures(object = object.list[[x]], 
                                                                                    assay = assay[x]))))
  var.features <- sort(x = table(var.features), decreasing = TRUE)
  gene <- NULL
  for(i in 1:length(object.list)){
    gene <- c(gene,rownames(object.list[[i]]))
  }
  # gene <- unname(obj = unlist(x = lapply(X = 1:length(x = object.list), 
  #                                        FUN = function(x) rownames(x = object.list[[X]][[assay[X]]]))))
  gene <- sort(x = table(gene), decreasing = TRUE)
  gene <- gene[gene>round(length(object.list)/2)]
  #for (i in 1:length(x = object.list)) {

    # var.features <- var.features[names(x = var.features) %in% 
    #                                rownames(x = object.list[[i]][[assay[i]]])]
  #}
  var.features <- var.features[names(x = var.features) %in% 
                                 names(gene)]
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > 
                                             tie.val)])
  vf.list <- lapply(X = object.list, FUN = VariableFeatures)
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = vf.list, FUN = function(vf) {
        if (x %in% vf) {
          return(which(x = x == vf))
        }
        return(NULL)
      })
      median(x = unlist(x = ranks))
    })
    features <- names(x = sort(x = feature.ranks))
  }
  features.tie <- var.features[which(x = var.features == tie.val)]
  tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
    ranks <- sapply(X = vf.list, FUN = function(vf) {
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- c(features, names(x = head(x = sort(x = tie.ranks), 
                                         nfeatures - length(x = features))))
  return(features)
}


data[["percent.mt"]] <- PercentageFeatureSet(data,pattern="^MT-")
summary(data$percent.mt)
data <- subset(data,subset=percent.mt<10)
dim(data)
objTN2 <- CreateSeuratObject(data@assays$RNA@counts,meta.data = data@meta.data,
                             min.cells = 3)
summary(objTN2$nFeature_RNA)
dim(objTN2)
rm(data)
rm(data2)
gc()

table(objTN2$dataid)
objTN2$dataset <- objTN2$dataid
objTN2$dataset <- ifelse(objTN2$dataset%in%c("PCall","PCall_new"),objTN2$cancer,objTN2$dataset)
table(objTN2$dataset)

objTN2.list <- SplitObject(objTN2,split.by = 'dataset')
objTN2.list <- lapply(objTN2.list,FUN=SCTransform)


features_2000 <- SelectIntegrationFeatures(objTN2.list,nfeatures = 2000)

gc()

feature_type <- c("features_2000")
objTN2_bak <- objTN2

for(i in 1:length(feature_type)){
  objTN2 <- objTN2_bak
  objTN2@assays$RNA@var.features <- get(feature_type[1])
  #remove noncoding
  genes<-data.frame(data.table::fread("./biomart_mart_export.txt",header=T,sep="\t"))
  
  table(genes$GeneType)
  genes_sub<-subset(genes,GeneType!="lncRNA") #processed_pseudogene lncRNA
  genes_sub<-subset(genes_sub,GeneType!="processed_pseudogene") #processed_pseudogene lncRNA
  genes_sub_all<-c(as.character(genes_sub$gene),as.character(genes_sub$GeneSynonym))
  sum(objTN2@assays$RNA@var.features %in% genes_sub_all); length(objTN2@assays$RNA@var.features)
  objTN2@assays$RNA@var.features = objTN2@assays$RNA@var.features[(objTN2@assays$RNA@var.features %in% genes_sub_all)]
  dim(objTN2)
  
  ####WYC blacklist
  load("./genes_black_WYC.rda")
  sum(objTN2@assays$RNA@var.features %in% genes_black); length(objTN2@assays$RNA@var.features)
  objTN2@assays$RNA@var.features = objTN2@assays$RNA@var.features[!(objTN2@assays$RNA@var.features %in% genes_black)]
  mito.genes <- rownames(objTN2@assays$RNA)[grep("^MT-",rownames(objTN2@assays$RNA))]
  objTN2@assays$RNA@var.features =objTN2@assays$RNA@var.features[!(objTN2@assays$RNA@var.features %in% mito.genes)]
  
  dim(objTN2);length(objTN2@assays$RNA@var.features)
  #setdiff(B_marker,objTN2@assays$RNA@var.features)
  
  # objTN2 <- SCTransform(objTN2)
  objTN2<- NormalizeData(objTN2)
  objTN2 <- ScaleData(objTN2)
  objTN2 <- RunPCA(objTN2, verbose = FALSE)
  
  objTN2 <- harmony::RunHarmony(objTN2,"patient", plot_convergence = TRUE)
  
  
  ElbowPlot(objTN2, ndims = 50,reduction = "harmony")
  
  objTN2 <- Seurat::RunUMAP(objTN2,reduction = "harmony", dims = 1:20) 
  objTN2 <- Seurat::FindNeighbors(objTN2,reduction = "harmony", dims = 1:20) 
  
  objTN2 <- Seurat::FindClusters(objTN2,resolution =2,graph.name="RNA_snn")
  
  DimPlot(objTN2,group.by = "cancer",cols = my36colors)
  DimPlot(objTN2,group.by = "seurat_clusters",label=T,cols = my36colors)

###############################################################################
#'                          Manuscipt: figure1C                              '#
###############################################################################
B_marker= c("TCL1A","FCER2","IL4R","IGHD" , ####NaiveB
            "IFIT3","IFI44L","STAT1","ISG15",###IFN
            "HSPA1A","DNAJB1",###Activated
            "MT1X","MT2A","SLC30A1",
            "EGR1","DUSP2",####ACB1
            "NR4A2","CD69","CD83",####ACB2
            "CCR7","PIM3","NFKBID",
            "S100A10","CRIP1","S100A4","ITGB1","CD27","CR2","AIM2","GPR183","CD1C",
            "DUSP4","FCRL5","ZEB2","ITGAX","FGR","FCRL4","CD274",
            "NME1","PSME2","ENO1","FABP5",###PreGC
            #"ACTG1","RGS13","PRPSAP2","MARCKSL1","ATP5L","LMO2",###GCB
            "CXCR4","GNB2L1","ATP5L","SUGCT",###DZGC and GC
            "LMO2","GMDS","PRPSAP2","MARCKSL1",###LZGC
            "STMN1","TUBB","HMGB2","TUBA1B","MKI67",
            # "STMN1","PTTG1","TYMS","MKI67","UBE2C",
            "JCHAIN","PRDM1","XBP1","MZB1"
            
)

DotPlot(object = objTN2, features =  B_marker,scale = T,group.by = "celltype_l3") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks= hypoxia,labels=  hypoxia)


###############################################################################
#'                          Manuscipt: figure1D&F                            '#
###############################################################################
load("color.B.Rdata")

db<-read.table("total_final230201_heavy_germ-pass.tsv",header=T,sep="\t")

db_obs<-observedMutations(db,
                          sequenceColumn="sequence_alignment", 
                          germlineColumn="germline_alignment_d_mask",
                          regionDefinition=IMGT_V_BY_REGIONS, 
                          frequency=T, 
                          combine = T,
                          nproc=20) 
db_obs[1:4,1:4]
rownames(db_obs)<-db_obs$cell_id

meta.data<-object@meta.data
length(as.vector(db_obs$cell_id))
length(rownames(object@meta.data))
y<-intersect(as.vector(db_obs$cell_id),object@meta.data$BCR_id)
length(y) 

rownames(db_obs)<-as.vector(db_obs$cell_id)
a=meta.data[(meta.data$BCR_id %in% y),]
a=a[!(duplicated(a$BCR_id)),]
b=db_obs[y,]

identical(a$BCR_id,rownames(b))
d=b[match(a$BCR_id,rownames(b)),]
identical(a$BCR_id,rownames(d))
new<-cbind(a,d)

object_BCR<-subset(object, cells = row.names(subset(object@meta.data, object@meta.data$BCR_id %in% y)))
object_BCR=object_BCR[,!(duplicated(object_BCR$BCR_id))]
dim(object_BCR)
object=object_BCR

###SHM
p1=ggplot(data = objectASC@meta.data,mapping = aes(x =type,y =mu_freq)) +
  geom_boxplot(mapping = aes(fill = type),scale = "width",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values=type)+
  labs(title = "IGHV total_Mut_freq")+#facet_wrap(~c_call,scales = "free_y",ncol=4)+
  xlab("Isotype") + ylab("Mutation frequency") +#coord_cartesian(ylim = c(0, 0.2))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p1

####Isotype
object$c_call %>% as.vector() %>% substr(1,4) ->object$c_call2
object@meta.data[object$c_call2=="",]$c_call2<-"Unknow"
Idents(object)<-"c_call2"
levels(object)
object<-subset(object,idents = c("IGHA","IGHD","IGHM","IGHG")) 

object$type=factor(object$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
objectB$c_call=factor(objectB$c_call,levels = c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHD","IGHM"))

p2=dittoBarPlot(object, "c_call",group.by="type",retain.factor.levels = T,main = "B",color.panel= color.B$BCR_col)+ylim(0,1);p2

###############################################################################
#'                          Manuscipt: figure1E                              '#
###############################################################################
###OR code
if(T){
  do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                            meta.cluster = cellInfo.tb$meta.cluster,
                            colname.patient = "patient",
                            loc = cellInfo.tb$loc,
                            out.prefix,
                            pdf.width=3,
                            pdf.height=5,
                            verbose=0){
    ##input data 
    library(data.table)
    dir.create(dirname(out.prefix),F,T)
    
    cellInfo.tb = data.table(cellInfo.tb)
    cellInfo.tb$meta.cluster = as.character(meta.cluster)
    
    if(is.factor(loc)){
      cellInfo.tb$loc = loc
    }else{cellInfo.tb$loc = as.factor(loc)}
    
    loc.avai.vec <- levels(cellInfo.tb[["loc"]])
    count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
    freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)
    
    {
      count.dist.melt.ext.tb <- test.dist.table(count.dist)
      p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
      OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
      OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
      rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
    }
    
    sscVis::plotMatrix.simple(round(OR.dist.mtx,2),
                              out.prefix=sprintf("%s.OR.dist",out.prefix),
                              show.number=T,
                              waterfall.row=T,
                              par.warterfall = list(score.alpha = 2,do.norm=T),
                              exp.name=expression(italic(OR)),
                              z.hi=3,
                              palatte=viridis::viridis(10),
                              pdf.width = pdf.width, pdf.height = pdf.height)
    if(verbose==1){
      return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                  "p.dist.tb"=p.dist.tb,
                  "OR.dist.tb"=OR.dist.tb,
                  "OR.dist.mtx"=OR.dist.mtx))
    }else{
      return(OR.dist.mtx)
    }
  }
  
  test.dist.table <- function(count.dist,min.rowSum=0)
  {
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
      this.row <- count.dist.melt.tb$rid[i]
      this.col <- count.dist.melt.tb$cid[i]
      this.c <- count.dist.melt.tb$count[i]
      other.col.c <- sum.col[this.col]-this.c
      this.m <- matrix(c(this.c,
                         sum.row[this.row]-this.c,
                         other.col.c,
                         sum(sum.col)-sum.row[this.row]-other.col.c),
                       ncol=2)
      res.test <- fisher.test(this.m)
      data.frame(rid=this.row,
                 cid=this.col,
                 p.value=res.test$p.value,
                 OR=res.test$estimate)
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    return(count.dist.melt.ext.tb)
  }
}

out.prefix <- "./"
OR.B.list <- do.tissueDist(cellInfo.tb=objN@meta.data,meta.cluster=objN@meta.data$celltype_l3,loc=objN@meta.data$type,
                           out.prefix=sprintf("%s/objN_type",out.prefix),
                           pdf.width=6,pdf.height=8,verbose=1)


##OR value
or=round(OR.B.list$OR.dist.mtx,2)
colnames(or)=colnames(OR.B.list$OR.dist.mtx)

or=ifelse(or >2,2,or)

sscVis::plotMatrix.simple(or,
                          out.prefix=sprintf("%s/objN_type",out.prefix),
                          show.number=T,
                          do.clust = F,
                          clust.column = F,
                          clust.row = F,
                          show.dendrogram=T,
                          waterfall.row=F,
                          par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(OR)),
                          z.hi=2,
                          palatte=brewer.pal(9, "YlGnBu"),
                          pdf.width = 6, pdf.height = 8)

                          
###############################################################################
#'                          Manuscipt: figureS2F&G                           '#
###############################################################################
# muscat ------------------------------------------------------------------

Atm <- readRDS("Atm_3type_panC_20cells_new_0922.rds")
DefaultAssay(Atm) <- "RNA"
Idents(Atm) <- Atm$TYPE
table(Atm$TYPE)
patient <- as.data.frame(table(Atm$patient))
#Cancer vs AID
Atm1 <- subset(Atm,idents=c("AID","panC"))
#Atm1 <- Atm
dim(Atm1)
table(Atm1$patient)
AtM_sce <- as.SingleCellExperiment(Atm1) 
AtM_sce$group_id <- AtM_sce$TYPE
AtM_sce$sample_id <- AtM_sce$patient
AtM_sce$cluster_id <- "AtM"
AtM_sce <- prepSCE(AtM_sce,
                   kid = "cluster_id",
                   sid = "sample_id",
                   gid = "group_id",
                   drop = FALSE)

table(AtM_sce$group_id,AtM_sce$sample_id)
pb <- aggregateData(AtM_sce,
                    assay="counts",
                    fun="sum",
                    by=c('cluster_id','sample_id'))
assayNames(pb)
t(head(assay(pb)))
metadata(pb)
#pb_mds <- pbMDS(pb)

res <- pbDS(pb,method = "DESeq2",min_cells = 10)
tmp <- AtM_sce
counts(tmp) <- as.matrix(counts(tmp))
result_table <- resDS(tmp,res,bind="row",frq=FALSE,spm=FALSE)
result_table$compare <- ifelse(result_table$logFC>0,"panC","AID")
write.csv(result_table, file = "result_table_panC_vs_AID_min10.csv", row.names = F)

#Cancer vs HBV
Atm1 <- subset(Atm,idents=c("HBV","panC"))
AtM_sce <- as.SingleCellExperiment(Atm1) 
AtM_sce$group_id <- AtM_sce$TYPE
AtM_sce$sample_id <- AtM_sce$patient
AtM_sce$cluster_id <- "AtM"
AtM_sce <- prepSCE(AtM_sce,
                   kid = "cluster_id",
                   sid = "sample_id",
                   gid = "group_id",
                   drop = FALSE)

table(AtM_sce$group_id,AtM_sce$sample_id)
pb <- aggregateData(AtM_sce,
                    assay="counts",
                    fun="sum",
                    by=c('cluster_id','sample_id'))
metadata(pb)
#pb_mds <- pbMDS(pb)
pb_mds 
res <- pbDS(pb,method = "DESeq2",min_cells = 10)
tmp <- AtM_sce
counts(tmp) <- as.matrix(counts(tmp))
result_table <- resDS(tmp,res,bind="row",frq=FALSE,spm=FALSE)

write.csv(result_table, file = "result_table_panC_vs_HBV_min10.csv", row.names = F)

#volcano

file="result_table_panC_vs_HBV_min10.csv"){

  #Atm
  Atm <- NULL
  for(i in 1:length(file)){
    temp <- read.csv(file[i],row.names = 1)
    #temp$compare <- colsplit(file[i],"marker",names=c("n1","n2"))$n1
    #temp$compare <- colsplit(temp$compare,"_VS_",names = c("n1","n2"))$n2
    temp$gene <- rownames(temp)
    temp <- temp[-grep("^MT-|^IGKV",temp$gene),]
    Atm <- rbind(Atm,temp)
  }
  Atm$cluster <- 'Atm'
  Atm$sig=""
  Atm$sig[abs(Atm$logFC) > 0.25 & Atm$p_adj.loc < 0.05] = "sig"
  Atm$sig2=paste(Atm$cluster,Atm$sig,sep = "_")
  Atm$sig2[str_detect(Atm$sig2,"_$")]="not_sig"
  Atm$sig2=str_replace(Atm$sig2,"_sig","")
  
  Atm$sig2=factor(Atm$sig2,levels = c("not",sort(unique(Atm$cluster))))
  Atm$cluster=factor(Atm$cluster,levels = sort(unique(Atm$cluster)))
  Atm=Atm%>%arrange(cluster,sig2)
  
  max(Atm$logFC);min(Atm$logFC)
  Atm$logFC[Atm$logFC > 7]=7
  Atm$logFC[Atm$logFC < c(-7)]= -7

  library(RColorBrewer)
  library(scales)
  color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],
             brewer.pal(5, "Set1")[2],
             brewer.pal(3, "Dark2")[1])
  names(color_ct)=sort(unique(as.character(Atm$cluster)))
  
  Atm %>% ggplot(aes(x=cluster,y=logFC,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
    scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
    scale_y_continuous("average log2FC",expand = c(0.02,0))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
      axis.text.y.left = element_text(size = 14,color = "black"),
      axis.title.x.bottom = element_blank(),
      axis.title.y.left = element_text(size = 16)
    )
  
  Atm2=Atm
  Atm2$padj_log10_neg= -log10(Atm2$p_adj.loc)
  Atm2$padj_log10_neg=ifelse(Atm2$logFC > 0,
                             Atm2$padj_log10_neg,
                             -Atm2$padj_log10_neg)
  
  Atm2$label <- ifelse(abs(Atm2$logFC)>0.8,Atm2$gene,"")
  
  Atm2 %>% ggplot(aes(x=cluster,y=logFC,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
    scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
    scale_y_continuous("average log2FC",expand = c(0.02,0))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
      axis.text.y.left = element_text(size = 14,color = "black"),
      axis.title.x.bottom = element_blank(),
      axis.title.y.left = element_text(size = 16)
    )+geom_label_repel( aes(label = label),
                        size = 3,box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE,
                        max.overlaps =50)
  
  ###################
  plot.list=list()
  for (ci in sort(unique(as.character(Atm2$cluster)))) {
    tmpdf=Atm2 %>% filter(cluster == ci)
    minabs=abs(min(tmpdf$padj_log10_neg))
    maxabs=max(tmpdf$padj_log10_neg)
    thre=0
    if(minabs < maxabs) {
      tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > minabs] = minabs
      thre=minabs
    }
    if(minabs > maxabs) {
      tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-maxabs)] = -maxabs
      thre=maxabs
    }
    if(minabs == maxabs & maxabs == Inf) {
      thre = min(
        abs(
          range(
            tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < Inf & tmpdf$padj_log10_neg > -Inf]
          )
        )
      )
      tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-thre)] = -thre
      tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > thre] = thre
    }
    
    plotdata = tmpdf
    tmpdf=tmpdf%>%filter(sig2 != "not") 
    tmpdf=tmpdf%>%arrange(desc(logFC))
    tmpdf.a=head(tmpdf%>%filter(logFC > 0),15)
    tmpdf.a$d=thre*2*0.05+(-thre)-tmpdf.a$padj_log10_neg
    tmpdf.b=tail(tmpdf%>%filter(logFC < 0),15)
    tmpdf.b$d=thre*2*0.95-thre  - tmpdf.b$padj_log10_neg
    textdata.down = tmpdf.b
    textdata.up   = tmpdf.a
    
    tmpplot=plotdata%>%ggplot(aes(x=padj_log10_neg,y=logFC))+
      geom_point(aes(color=sig2),size=1)+
      geom_hline(yintercept = c(-0.25,0.25),linetype="dashed")+
      geom_text_repel(data = textdata.down,
                      mapping = aes(label=gene),
                      nudge_x=textdata.down$d,
                      direction = "y", hjust = 1,segment.size = 0.2,max.overlaps = 100)+
      geom_text_repel(data = textdata.up,
                      mapping = aes(label=gene),
                      nudge_x=textdata.up$d,
                      direction = "y", hjust = 0,segment.size = 0.2,max.overlaps = 100)+
      labs(title = ci)+
      scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
      scale_y_continuous("average log2FC",expand = c(0.02,0),limits = c(-8,8))+
      theme_bw()+
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks.x.bottom = element_blank(),
        #axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.y.left = element_text(size = 16),
        plot.title = element_text(size = 12,hjust = 0.5)
      )
    
    index=which(ci == sort(unique(as.character(Atm2$cluster))))
    if (index!=1) {
      tmpplot=tmpplot+theme(
        axis.title.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.left = element_blank()
      )
    }
    if (index == length(sort(unique(as.character(Atm2$cluster))))) {
      segment.df=data.frame(x=c(0 - thre / 5,0 + thre / 5),
                            xend=c(-thre,thre),
                            y=c(-3,-3),
                            yend=c(-3,-3))
      tmpplot=tmpplot+geom_segment(data = segment.df,
                                   mapping = aes(x=x,xend=xend,y=y,yend=yend),
                                   arrow = arrow(length=unit(0.3, "cm")))
      
    }
    plot.list[[get("index")]]=tmpplot
  }
  wrap_plots(plot.list,ncol = 1)&theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
  ggsave("Atm_diffgene_HBV_mincell10_20231011.pdf",width = 4,height = 8)


#GSVA

Atm <- readRDS("Atm_3type_panC_20cells_new_0922.rds")
table(Atm$TYPE)

Idents(Atm) <- Atm$TYPE

##########AID###############
AID <- subset(Atm,idents = c("AID","panC"))
expr <- AID@assays$RNA@counts
dir.create("./review/DEG/GSVA")


genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
head(genesets)


gsva.res<- gsva(expr,genesets, method="gsva",parallel.sz=16,kcdf="Poisson") 
saveRDS(gsva.res, "gsva_hallmarker_AID.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_hallmarker_AID.csv", row.names = F)
gsva.df[1:3,1:3]

gsva.res=gsva.res
gsva=gsva.df


ac <- data.frame(group=AID$TYPE) 
design <- model.matrix(~ 0 + factor(ac$group))
colnames(design) <- levels(factor(ac$group))
# rownames(design) <- colnames(sub_regulonAUC)
head(design)

library(limma)
contrast.matrix <- makeContrasts(panC - AID, levels = design)


fit <- lmFit(gsva.res, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)



pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
head(df)
#write.csv(df, "./GSVA/enrich_B_Hall.csv", quote = F, row.names = F)


cutoff <- 0
df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf),labels = c(1,2))
head(df)

sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
top2=rbind(head(sortdf,n=10),tail(sortdf,n=10))

##  
p1=ggplot(top2, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('#bf79a8', '#f2cdc1'), guide = FALSE) + 
  
  geom_hline(yintercept = c(-1,1), 
             color="white",
             linetype = 2, 
             size = 0.3) + 
  
  geom_text(data = subset(top2, score < 0),
            aes(x=ID, y= 0.1, label=ID, color = group),
            size = 4, 
            hjust = "inward" ) +  
  geom_text(data = subset(top2, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = group),
            size = 4, hjust = "outward") +  
  scale_colour_manual(values = c("black","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank());p1 

p1


ggsave("gsva_hallmarker_AID_TOP20.pdf", width = 10, height = 8)


############HBV#############
HBV <- subset(Atm,idents = c("HBV","panC"))
expr <- HBV@assays$RNA@counts
#dir.create("./diffgene/GSVA")

###
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
head(genesets)

#
gsva.res<- gsva(expr,genesets, method="gsva",parallel.sz=20) 
saveRDS(gsva.res, "gsva_hallmarker_HBV.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_hallmarker_HBV.csv", row.names = F)
gsva.df[1:3,1:3]


gsva.res=gsva.res
gsva=gsva.df
olnames(objectTPN)


ac <- data.frame(group=HBV$TYPE) 
design <- model.matrix(~ 0 + factor(ac$group))
colnames(design) <- levels(factor(ac$group))
# rownames(design) <- colnames(sub_regulonAUC)
head(design)


library(limma)

contrast.matrix <- makeContrasts(panC-HBV, levels = design)

fit <- lmFit(gsva.res, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)


pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
head(df)
#write.csv(df, "./GSVA/enrich_B_Hall.csv", quote = F, row.names = F)

cutoff <- 0
df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf),labels = c(1,2))
head(df)

sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
top2=rbind(head(sortdf,n=10),tail(sortdf,n=10))


p1=ggplot(top2, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('#bf79a8', '#f2cdc1'), guide = FALSE) + 

  geom_hline(yintercept = c(-1,1), 
             color="white",
             linetype = 2, 
             size = 0.3) + 
  geom_text(data = subset(top2, score < 0),
            aes(x=ID, y= 0.1, label=ID, color = group),
            size = 4, 
            hjust = "inward" ) +  
  geom_text(data = subset(top2, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = group),
            size = 4, hjust = "outward") +  
  scale_colour_manual(values = c("black","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank());p1 

p1

ggsave("gsva_hallmarker_HBV_TOP20.pdf", width = 10, height = 8)

###############################################################################
#'                          Manuscipt: figureS3A                             '#
###############################################################################
diversity.norm <- function(x)
{
  -sum(ifelse(x>0,x*log2(x),0))/log2(length(x))
}
dat.plot=freqB
dat.plot$group="B"
dat.plot$group=ifelse(dat.plot$group.var %in% c( "B14.PB","B15.PC"),"PB","B")#,   

dat.plot$type=factor(dat.plot$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
dat.diversity.norm.perSample.tb <- dat.plot[,.(NMCls=.N,
                                               NTotal=NTotal[1],
                                               diversity=diversity.norm(freq)),
                                            by=c("group","cmp.var","donor.var","type")]




p <- ggboxplot(dat.diversity.norm.perSample.tb,x="type",y="diversity",xlab="",
               add = "jitter",outlier.shape=NA,legend="none",
               color="type") +
  scale_color_manual(values=type) +
  facet_wrap(~group,scales = "free_y",ncol=3)+
  stat_compare_means(comparisons=my_comparisons, method="wilcox.test",label = "p.signif") +
  #facet_wrap(~group,nrow=1) +
  theme(strip.background=element_blank(),strip.text=element_text(size=12),axis.text.x = element_text(angle=45, hjust=1, vjust=1));p

###############################################################################
#'                          Manuscipt: figureS3B                             '#
###############################################################################
####Diversity analysis###########
Idents(obj_bcr)="celltype_level3_0912";table(Idents(obj_bcr))
objectnPC=subset(obj_bcr,idents=c(c("B14.PB","B15.PC")));table(objectnPC$celltype_level3_0912)

Idents(objectnPC)="groupn";table(Idents(objectnPC))
objectnPC=subset(objectnPC,idents= c("EF",     "GC"))

table(objectnPC$celltype_level3_0912)
new<-objectnPC@meta.data
colnames(new)
select(new,c("cell_id","c_call","clone_id")) %>% head()
new$clone_id<-as.vector(new$clone_id)
q<-as.data.frame(table(new$clone_id))
q<-q[order(q$Freq,decreasing = T),]
# get top 10 genes
clone_id_150=head(q,n=150)$Var1

top150=new[new$clone_id %in% clone_id_150,]


# Partitions the data on the sample column
# Calculates a 95% confidence interval via 200 bootstrap realizations
curve <- estimateAbundance(top150, group="group", ci=0.95, nboot=1000, clone="clone_id")
# Plots a rank abundance curve of the relative clonal abundances
type=c( "Blood"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
class=c("GC"="#377eb8" ,"EF"="#e41a1c")
plot(curve, colors = class, legend_title="Type diversity")

ggsave("Top150 type objectASC(PB_PC) EF vs GC abundance231007.pdf",width = 4,height = 5)

sample_curve <- alphaDiversity(top150, group="groupn", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,min_n=10,
                               ci=0.95, nboot=1000)

sample_main <- paste0("Type diversity")
type=c( "Blood"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
plot(sample_curve, colors=class, main_title=sample_main, 
     legend_title="Type diversity")


ggsave("Top150 type objectASC(PB_PC) EF vs GC diversity231007.pdf",width = 4,height = 5)


###############################################################################
#'                          Manuscipt: figureS3C                             '#
###############################################################################
object <- readRDS("Pancancer_all_remove_use_final_dim18_20230921.rds")
anno <- object@meta.data
rm(object)
gc()

ASC <- readRDS("ASC_rmtype_20230922.rds")
cellname.ASC <- colnames(ASC)
rm(ASC)
gc()
# anno.PB <- anno[anno$celltype_level3_0912%in%c("c14_PB"),]
anno.ASC <- anno[anno$celltype_level3_0912%in%c("c15_MZB1+ASC"),]
anno.ASC <- anno.ASC[cellname.ASC,]
anno.other <- anno[anno$celltype_level3_0912!=c("c15_MZB1+ASC"),]
anno.use <- rbind(anno.ASC,anno.other)

#pb
use <- anno.use[,c("type","celltype_level3_0912","patient")]
use$type <- gsub("^T$","Cancer",use$type)
use <- as.data.frame(table(use$patient,use$type,use$celltype_level3_0912))
use2 <- matrix(nrow = length(unique(use$Var1)),ncol=5)
colnames(use2) <- c("patient","Blood","Adjacent","LN_Met","Cancer")
use2 <- as.data.frame(use2)
for(i in 1:length(unique(use$Var1))){
  temp <- use[use$Var1==unique(use$Var1)[i],]
  use2[i,1] <- as.character(unique(temp$Var1))
  for(j in 1:length(unique(temp$Var2))){
    temp2 <- temp[temp$Var2==unique(temp$Var2)[j],]
    if(sum(temp2$Freq)==0){
      next()
    }
    else{
      freq <- temp2$Freq[temp2$Var3=="c14_PB"]/(sum(temp2$Freq))
      use2[i,as.character(unique(temp$Var2)[j])] <- freq
    }
  }

}

use2.melt <- melt(use2)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
col_flg <- col_flg[c(3,2,4,1)]
use2.melt <- na.omit(use2.melt)
p1 <- ggplot(use2.melt,aes(x=variable,y=value*100,fill=variable))+geom_boxplot()+
  theme_classic()+xlab("")+ylab("Frequency")+
  scale_fill_manual(values=col_flg)+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("c14_PB")+
  #geom_line(aes(group=patient),color='gray',lwd=0.1)+
  stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Blood","Adjacent"),
                                                                c("Blood","LN_Met"),
                                                                c("Blood","Cancer"),
                                                                c("Adjacent","LN_Met"),
                                                                c("Adjacent","Cancer"),
                                                                c("LN_Met","Cancer")))
p1


#ASC
use <- anno.use[,c("type","celltype_level3_0912","patient")]
use$type <- gsub("^T$","Cancer",use$type)
use <- as.data.frame(table(use$patient,use$type,use$celltype_level3_0912))
use2 <- matrix(nrow = length(unique(use$Var1)),ncol=5)
colnames(use2) <- c("patient","Blood","Adjacent","LN_Met","Cancer")
use2 <- as.data.frame(use2)
for(i in 1:length(unique(use$Var1))){
  temp <- use[use$Var1==unique(use$Var1)[i],]
  use2[i,1] <- as.character(unique(temp$Var1))
  for(j in 1:length(unique(temp$Var2))){
    temp2 <- temp[temp$Var2==unique(temp$Var2)[j],]
    if(sum(temp2$Freq)==0){
      next()
    }
    else{
      freq <- temp2$Freq[temp2$Var3=="c15_MZB1+ASC"]/(sum(temp2$Freq))
      use2[i,as.character(unique(temp$Var2)[j])] <- freq
    }
  }
  
}

use2.melt <- melt(use2)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
col_flg <- col_flg[c(3,2,4,1)]
use2.melt <- na.omit(use2.melt)
p2 <- ggplot(use2.melt,aes(x=variable,y=value*100,fill=variable))+geom_boxplot()+
  theme_classic()+xlab("")+ylab("Frequency")+
  scale_fill_manual(values=col_flg)+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("c15_MZB1+ASC")+
  #geom_line(aes(group=patient),color='gray',lwd=0.1)+
  stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Blood","Adjacent"),
                                                                c("Blood","LN_Met"),
                                                                c("Blood","Cancer"),
                                                                c("Adjacent","LN_Met"),
                                                                c("Adjacent","Cancer"),
                                                                c("LN_Met","Cancer")))
p2


p1|p2
ggsave("FigS3C_20230923.pdf",width = 10,height = 5)

###############################################################################
#'                          Manuscipt: figureS3E                             '#
###############################################################################
set.seed(123)
objN <- readRDS("Pancancer_all_remove_use_final_dim18_20230921.rds")
input.num = 100000
cellid<-sample(1:ncol(objN), input.num, replace=F); length(cellid)
obj_random2w<-objN[,cellid]
dim(obj_random2w)
obj.2_expr=obj_random2w@assays$RNA@counts
metadata=obj_random2w@meta.data

ent.res <- SE_fun(obj.2_expr)
head(ent.res)
SEplot(ent.res)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

dim(obj.2_expr)
dim(metadata)
rogue.res <- rogue(obj.2_expr, labels = metadata$celltype_level3_0912,
                   samples = metadata$patient, platform = "UMI", span = 1.2)
rogue.res=rownames_to_column(rogue.res)
rogue.res1=rogue.res[!is.na(rogue.res),]
rogue.res2=rogue.res1[,-1]

rogue.boxplot(rogue.res2)


rogue.res1=melt(rogue.res1,
                id.vars="rowname",
                variable.names="cluster",
                value.name = "ROGUE")
rogue.res2=na.omit(rogue.res1)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

#my_comparisons=list(c("EF","GC"))
unique(rogue.res2$variable)
rogue.res2$variable <- factor(rogue.res2$variable,levels = c("c01_TCL1A+naiveB",
                                                             "c02_IFIT3+B","c03_HSP+B",
                                                             "c04_MT1X+B","c05_EGR+ACB",
                                                             "c06_NR4A2+ACB","c07_CCR7+ACB",
                                                             "c08_ITGB1+SwBm","c09_AtM",
                                                             "c10_PreGC","c11_DZ GCB",
                                                             "c12_LZ GCB","c13_Cycling GCB",
                                                             "c14_PB","c15_MZB1+ASC"))

mycompare <- list(c("c01_TCL1A+naiveB","c14_PB"),
                  c("c02_IFIT3+B","c14_PB"),
                  c("c03_HSP+B","c14_PB"),
                  c("c04_MT1X+B","c14_PB"),
                  c("c05_EGR+ACB","c14_PB"),
                  c("c06_NR4A2+ACB","c14_PB"),
                  c("c07_CCR7+ACB","c14_PB"),
                  c("c08_ITGB1+SwBm","c14_PB"),
                  c("c09_AtM","c14_PB"),
                  c("c10_PreGC","c14_PB"),
                  c("c11_DZ GCB","c14_PB"),
                  c("c12_LZ GCB","c14_PB"),
                  c("c13_Cycling GCB","c14_PB"))
mycompare2 <- list(c("c01_TCL1A+naiveB","c15_MZB1+ASC"),
                  c("c02_IFIT3+B","c15_MZB1+ASC"),
                  c("c03_HSP+B","c15_MZB1+ASC"),
                  c("c04_MT1X+B","c15_MZB1+ASC"),
                  c("c05_EGR+ACB","c15_MZB1+ASC"),
                  c("c06_NR4A2+ACB","c15_MZB1+ASC"),
                  c("c07_CCR7+ACB","c15_MZB1+ASC"),
                  c("c08_ITGB1+SwBm","c15_MZB1+ASC"),
                  c("c09_AtM","c15_MZB1+ASC"),
                  c("c10_PreGC","c15_MZB1+ASC"),
                  c("c11_DZ GCB","c15_MZB1+ASC"),
                  c("c12_LZ GCB","c15_MZB1+ASC"),
                  c("c13_Cycling GCB","c15_MZB1+ASC"))
ggboxplot(rogue.res2,x="variable",y="ROGUE",
          color = "variable",palette = my36colors,
          add = "jitter")+ 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+xlab("")+
  stat_compare_means(comparisons = mycompare,aes(label=..p.signif..))

ggboxplot(rogue.res2,x="variable",y="ROGUE",
          color = "variable",palette = my36colors,
          add = "jitter")+ 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+xlab("")+
  stat_compare_means(comparisons = mycompare2,aes(label=..p.signif..))

p=ggboxplot(rogue.res2,x="variable",y="ROGUE",
            color = "variable",palette = my36colors,
            add = "jitter")+ 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  coord_cartesian(ylim = c(0, 1))+xlab("");p

ggsave("FigS3E.pdf",width = 10,height = 6)

###############################################################################
#'                          Manuscipt: figureS3F                             '#
###############################################################################
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "none",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(0, 4), "lines")
    )
}
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

umapColor <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
my12colors<-c("#9ED9EA","#B1D896","#DAC0DD","#FDC37E","#F7A7A8","#008FC5","#ABA9AA","#06AC4B","#F68C1F","#D897C2","#F9AC8A","#42BC99")
maumapColor=my12colors
objN_ASC <- readRDS("ASC_rmtype_20230922.rds")
object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
anno <- object@meta.data
anno.EF <- anno[anno$groupn=="EF",]
anno.GC <- anno[anno$groupn=="GC",]
objN_ASC$groupn <- ifelse(colnames(objN_ASC)%in%rownames(anno.EF),"EF",
                          ifelse(colnames(objN_ASC)%in%rownames(anno.GC),"GC","others"))
table(objN_ASC$groupn)
DimPlot(objN_ASC,split.by = 'groupn')

Idents(objN_ASC) <- objN_ASC$groupn
objN_ASC <- subset(objN_ASC,idents = c("EF","GC"))
allObj <- list(ASC = objN_ASC)
#allObj <- list(ASC = objN)

Idents(objN_ASC)="celltype_l4";table(Idents(objN_ASC))
figurePath <- c("./ASC")
for(oi in 1:length(allObj)){
  tempName <- names(allObj[oi])
  tempObj <- allObj[[oi]]
  Idents(tempObj) <- tempObj$seurat_clusters
  coord = Embeddings(object = tempObj, reduction = "umap")
  coord = coord[,c(1,2)]
  colnames(coord) = c("UMAP_1", "UMAP_2")
  coord = data.frame(ID = rownames(coord), coord)
  meta = tempObj@meta.data
  meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
  meta = left_join(meta, coord, by = 'ID')
  # randomly sample cells of large tissuetype ###############################
  print(tempName)
  print(table(meta$groupn))
  minNum <- min(table(meta$groupn))
  print("minNum")
  print(minNum)
  meta <- meta %>%
    group_by(groupn) %>%
    slice_sample(n = minNum)
  for(ti in 1: length(unique(meta$groupn))){
    tin <- sort(unique(meta$groupn))[ti]
    meta_tin <- meta[meta$groupn == tin,]
    meta_ntin <- meta[meta$groupn != tin,]
    g <- ggplot(data = meta_tin, mapping = aes(x = UMAP_1, y = UMAP_2)) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) +
      geom_point(color = 'white', size = .05) +
      scale_fill_viridis(option="magma") +
      theme_black()
    pdf(file.path(figurePath, paste0(tempName, "_type-", tin, "_umap_density.pdf")),width = 10,height = 8)
    print(g)
    dev.off()
    g <- ggplot(data = meta_tin, mapping = aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = celltype_l4), size = 0.5) +
      scale_color_manual(values = umapColor) +
      theme_void() +
      theme(text = element_text(size = 6),
            legend.position = "none")
    pdf(file.path(figurePath, paste0(tempName, "_type-",tin, "_umap.pdf")),width = 10,height = 8)
    print(g)
    dev.off()
  }
}


###############################################################################
#'                          Manuscipt: figureS3G                             '#
###############################################################################
data_raw <- readRDS("Pancancer_all_remove_use_final_dim18_20230921.rds")
cellname <- sample(colnames(data_raw),50000)
data.sample <- data_raw[,cellname]
DimPlot(data.sample,label=T,group.by = "celltype_level3_0912")
Idents(data.sample) <- data.sample$celltype_level3_0912

anno_sample <- data.sample@meta.data
###monocle3
data <- GetAssayData(data.sample,assay = "RNA",slot = "counts")
data <- data[rowSums(data>0)>=3,]
dim(data)
cell_metadata <- data.sample@meta.data
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 50)

cds <- reduce_dimension(cds,preprocess_method = "PCA")
p1 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "type")
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data.sample,reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "type")
p2

cds <- cluster_cells(cds,resolution = 0.01,k=20,random_seed=18,verbose=T)
cds@clusters$UMAP$clusters <- data.sample$celltype_level3_0912
plot_cells(cds,color_cells_by = "partition")

p1 <- plot_cells(cds,group_cells_by = 'cluster')
p1

cds <- learn_graph(cds, verbose =T,
                   use_partition=T,close_loop=F)
p <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=T,
                label_leaves=T, label_branch_points=T,cell_size = 0.5,group_label_size=4)

p
get_earliest_principal_node <- function(cds, time_bin="2"){
  cell_ids <- which(cds@clusters@listData[["UMAP"]][["clusters"]] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

cds <- order_cells(cds)

p1 <- plot_cells(cds,color_cells_by = "pseudotime",label_branch_points = FALSE,label_leaves = F)
p1
p2 <- plot_cells(cds,color_cells_by = "celltype_level3_0912",
                 label_branch_points = FALSE,label_leaves = F)+
  scale_color_manual(values=my36colors)
p2|p1

ggsave("monocle3_0922_FigS3E.pdf",width = 20,height = 8)

###############################################################################
#'                         Manuscipt: figureS3J&K                            '#
###############################################################################
# (Jaccard index)

BCR <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_20230921.rds")
BCR@meta.data[1:3,]

jaccard_TCR = function(TCR1, TCR2){  
  intersection <- length(intersect(TCR1, TCR2))
  union <- length(union(TCR1, TCR2))
  jaccard_index <- intersection / union
  return(jaccard_index)
}


metadata = BCR@meta.data
all.ident = unique(metadata$celltype_level3_0912)
all.type = unique(metadata$type)
jaccard_result = c()

for (j in 1:length(all.ident)){
  for (i in 1:length(all.ident)){
    for (k in 1:length(all.type)){
      
      
      metadata_tmp = subset(metadata, metadata$type == all.type[k])
      
      TCR1 <- subset(metadata_tmp, metadata_tmp$celltype_level3_0912 == all.ident[j])$cdr3
      TCR2 <- subset(metadata_tmp, metadata_tmp$celltype_level3_0912 == all.ident[i])$cdr3
      
      tmp = c(all.ident[j], all.ident[i], jaccard_TCR(TCR1, TCR2), all.type[k])
      jaccard_result = rbind(jaccard_result, tmp)
    }
  }
}
jaccard_result = data.frame(jaccard_result)
colnames(jaccard_result) = c("cell1", "cell2", "jaccard","type")
jaccard_result$jaccard = as.numeric(as.character(jaccard_result$jaccard))
jaccard_result_bak = jaccard_result

jaccard_result$jaccard = round(jaccard_result$jaccard, 3)

jaccard_result = subset(jaccard_result, jaccard_result$jaccard != 1)
jaccard_result = subset(jaccard_result, jaccard_result$cell1 !=  jaccard_result$cell2)

cc <- colorRampPalette(c("#352a86", "#095cd8", "#46b896", "#e7ba4a", "#f8fa0d"))
cc <- colorRampPalette(c("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
# cc = colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
cc <- colorRampPalette(c("grey90", "#fc58a6"))
cc <- colorRampPalette(c("white", "white", "#F47E5D", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"
cc <- colorRampPalette(c("white", "#F47E5D", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"
# cc <- colorRampPalette(c("white", "#F47E5D", "#CA3D74", "#7F2880", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"

ggplot(jaccard_result, aes(cell1, cell2, fill = jaccard)) + 
  geom_tile(aes(fill = jaccard, fill = jaccard),size=1)+
  # geom_point(aes(size = pct.exp, fill =  avg.exp.scaled, color = avg.exp.scaled), shape = 21, colour = "black")+
  scale_fill_gradientn(colours = (cc(100))) +
  # scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  # scale_color_gradient2(low = "#2b8cbe",mid = "gray90",high = my12colors[3])+ ## ("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
  geom_text(aes(label=jaccard),col ="black",size = 3)+
  theme_minimal()+# 
  theme(axis.title.x=element_blank(),#
        axis.ticks.x=element_blank(),#
        axis.title.y=element_blank(),#
        axis.text.x = element_text(angle = 45, hjust = 1),# 
        axis.text.y = element_text(size = 8),
        aspect.ratio=1)+#
  #
  facet_wrap(~ type, ncol = 2, scales = "free") +
  labs(fill =paste0("Jaccard Index",""))



############
# overlap BCR

BCR <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_20230921.rds")
BCR@meta.data[1:3,]


metadata = BCR@meta.data
all.ident = unique(metadata$celltype_level3_0912)
all.type = unique(metadata$type)
# metadata = subset(metadata, type %in% c("Cancer",))



stat_bcr = c()
for (j in 1:length(all.ident)){
  for (i in 1:length(all.ident)){
    metadata_1 = subset(metadata, celltype_level3_0912 == all.ident[j])
    metadata_2 = subset(metadata, celltype_level3_0912 == all.ident[i])
    tmp = intersect(metadata_1$cdr3, metadata_2$cdr3)
    
    number_overlap = length(tmp)
    number_celltype1 = nrow(metadata_1)
    number_celltype2 = nrow(metadata_2)
    stat_bcr = rbind(stat_bcr, c(number_overlap, number_celltype1, number_celltype2, all.ident[j], all.ident[i]))
  }
}
stat_bcr = data.frame(stat_bcr)
colnames(stat_bcr) = c("number_overlap", "number_celltype1", "number_celltype2", "celltype1", "celltype2")
stat_bcr[,1] = as.numeric(as.character(stat_bcr[,1]));stat_bcr[,2] = as.numeric(as.character(stat_bcr[,2]));stat_bcr[,3] = as.numeric(as.character(stat_bcr[,3]));
stat_bcr$prop_overlap = stat_bcr$number_overlap/(stat_bcr$number_celltype1 + stat_bcr$number_celltype2)

stat_bcr = subset(stat_bcr, stat_bcr$celltype1 != stat_bcr$celltype2)


stat_bcr$celltype_merge = paste(stat_bcr$celltype1, stat_bcr$celltype2, sep="-")

stat_bcr$celltype_merge = as.character(stat_bcr$celltype_merge)
stat_bcr = stat_bcr[order(stat_bcr$prop_overlap, decreasing = T),]
stat_bcr$celltype_merge <- factor(stat_bcr$celltype_merge, levels=c(stat_bcr$celltype_merge))
row.names(stat_bcr) = 1:nrow(stat_bcr)

even_columns <- seq(2, nrow(stat_bcr), by = 2)
stat_bcr = stat_bcr[even_columns,]
row.names(stat_bcr) = 1:nrow(stat_bcr)


# stat_bcr = subset(stat_bcr, stat_bcr$celltype_merge %like% "c15_MZB1+ASC" | stat_bcr$celltype_merge %like% "c14_PB" | stat_bcr$celltype_merge %like% "c09_AtM" )
ggplot(stat_bcr, aes(y=prop_overlap, x=celltype_merge)) +
  geom_bar(position="stack", stat="identity")+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


stat_bcr = subset(stat_bcr, stat_bcr$celltype_merge %like% "c09_AtM")



unique(BCR$celltype_level3_0912)


color_select = c("c01_TCL1A+naiveB"=my36colors[1],
                 "c02_IFIT3+B"=my36colors[2],
                 "c03_HSP+B"=my36colors[3],
                 "c04_MT1X+B"=my36colors[4],
                 "c05_EGR+ACB"=my36colors[5],
                 "c06_NR4A2+ACB"=my36colors[6],
                 "c07_CCR7+ACB"=my36colors[7],
                 "c08_ITGB1+SwBm"=my36colors[8],
                 "c09_AtM"=my36colors[9],
                 "c10_PreGC"=my36colors[10],
                 "c11_DZ GCB"=my36colors[11],
                 "c12_LZ GCB"=my36colors[12],
                 "c13_Cycling GCB"=my36colors[13],
                 "c14_PB"=my36colors[14],
                 "c15_MZB1+ASC"=my36colors[15]
)

stat_bcr[1:3,]
ggplot(stat_bcr, aes(y=prop_overlap, x=celltype_merge, fill = celltype1)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = color_select)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(stat_bcr, aes(y=prop_overlap, x=celltype_merge, fill = celltype2)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = color_select)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


#show number
stat_bcr$total_number = stat_bcr$number_celltype1 + stat_bcr$number_celltype2
ggplot(stat_bcr, aes(y=prop_overlap, x=celltype_merge, fill = celltype2)) +
  geom_text(aes(label=number_overlap),col ="black",size = 3)+
  geom_text(aes(label=total_number),col ="black",size = 3)+
  # geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = color_select)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


############
# overlap BCR ids


BCR <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_20230921.rds")
BCR@meta.data[1:3,]

metadata = BCR@meta.data
all.ident = unique(metadata$celltype_level3_0912)
all.type = unique(metadata$type)
all.patient = unique(metadata$patient)

metadata$patient_clonetype = paste(metadata$patient, metadata$clone_id, sep="_")
all.patient_clonetype = unique(metadata$patient_clonetype)
stat.patient_clonetype = data.frame(table(metadata$patient_clonetype))

#

all.pairs = list(c("AtM","ASC"), c("Bm","ASC"), c("GCB","ASC"))

all.pairs = list(
  c("c14_PB","c15_MZB1+ASC"),
  c("c09_AtM","c15_MZB1+ASC"),
  c("c13_Cycling GCB","c14_PB"),
  c("c12_LZ GCB","c14_PB"),
  c("c11_DZ GCB","c14_PB"),
  c("c14_PB","c09_AtM"),
  
  c("c07_CCR7+ACB","c09_AtM"),
  c("c09_AtM","c06_NR4A2+ACB"),
  c("c05_EGR+ACB","c09_AtM"),
  c("c09_AtM","c03_HSP+B"),
  c("c02_IFIT3+B","c09_AtM"),
  c("c09_AtM","c10_PreGC")
)



stat_merge_large = c()
for (j in 1:length(all.pairs)){
  print(j)
  
  metadata_tmp = subset(metadata, metadata$celltype_level3_0912 %in% all.pairs[[j]])
  stat.patient_clonetype = data.frame(table(metadata_tmp$patient_clonetype))
  stat.patient_clonetype = subset(stat.patient_clonetype, stat.patient_clonetype$Freq > 1)
  all.patient_clonetype = as.character(unique(stat.patient_clonetype$Var1))
  
  stat_merge = data.frame(Var1 = "Cancer")
  for (k in 1:length(all.patient_clonetype)){
    if (k %% 200 ==0) print(k)
    metadata_tmp2 = subset(metadata_tmp, metadata_tmp$patient_clonetype == all.patient_clonetype[k])
    stat_tmp = data.frame(table(metadata_tmp2$type))
    colnames(stat_tmp)[2] = all.patient_clonetype[k]
    stat_merge = merge(stat_merge, stat_tmp, by = "Var1", all = TRUE)
  }
  stat_merge = data.frame(stat_merge, row.names = 1)
  stat_merge_t = data.frame(t(stat_merge))
  stat_merge_t$celltype = capture.output(cat(all.pairs[[j]]))
  stat_merge_large = data.frame(rbind(stat_merge_large, stat_merge_t))
}

stat_merge_large_bak = stat_merge_large

stat_merge_large_sub = subset(stat_merge_large, rowSums(is.na(stat_merge_large)) < 3)


stat_merge_large_sub = stat_merge_large_sub[order(stat_merge_large_sub$LN_Met, decreasing = T),]
stat_merge_large_sub = stat_merge_large_sub[order(stat_merge_large_sub$Adjacent, decreasing = T),]
stat_merge_large_sub = stat_merge_large_sub[order(stat_merge_large_sub$Blood, decreasing = T),]
stat_merge_large_sub = stat_merge_large_sub[order(stat_merge_large_sub$Cancer, decreasing = T),]
stat_merge_large_sub = stat_merge_large_sub[order(stat_merge_large_sub$celltype, decreasing = T),]

stat_merge_large_sub_heat = stat_merge_large_sub[,1:4]
stat_merge_large_sub_heat[is.na(stat_merge_large_sub_heat)] = 0
stat_merge_large_sub_heat[stat_merge_large_sub_heat > 8] = 8


library(pheatmap)
cc <- colorRampPalette(c("white", "#F47E5D", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"

save_pheatmap_pdf <- function(x, filename, width=5, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

annotation_row2 = data.frame(cell = row.names(stat_merge_large_sub), celltype = stat_merge_large_sub$celltype, row.names =1)

all.celltype = unique(annotation_row2$celltype)
for (j in 1:length(all.celltype)){
  annotation_row3 = subset(annotation_row2, annotation_row2$celltype == all.celltype[j])
  
  stat_merge_large_sub_heat2 = stat_merge_large_sub_heat[row.names(annotation_row3),c("Cancer","Blood","Adjacent","LN_Met")]
  
  ptmp = pheatmap(stat_merge_large_sub_heat2, cluster_rows = F, cluster_cols = F,
                  color = cc(100),
                  show_rownames = F,
                  annotation_row = annotation_row2,
                  border_color = NA
  )
  save_pheatmap_pdf(ptmp, paste("./", all.celltype[j],".pdf", sep=""))
}


