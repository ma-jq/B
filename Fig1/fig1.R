if(T){
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
  library(alakazam)
  library(UpSetR)
}


###############################################################################
#'                          Manuscipt: figure1B                              '#
###############################################################################

###We modify SelectIntegrationFeatures function from Seurat package
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
# dim(objTN2)
# rm(data)
# rm(data2)
# gc()

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


objTN2 <- objTN2_bak
objTN2@assays$RNA@var.features <- get(feature_type[1])
#remove noncoding
genes<-data.frame(data.table::fread("./Additional_data/biomart_mart_export.txt",header=T,sep="\t"))
  
table(genes$GeneType)
genes_sub<-subset(genes,GeneType!="lncRNA") #processed_pseudogene lncRNA
genes_sub<-subset(genes_sub,GeneType!="processed_pseudogene") #processed_pseudogene lncRNA
genes_sub_all<-c(as.character(genes_sub$gene),as.character(genes_sub$GeneSynonym))
sum(objTN2@assays$RNA@var.features %in% genes_sub_all); length(objTN2@assays$RNA@var.features)
objTN2@assays$RNA@var.features = objTN2@assays$RNA@var.features[(objTN2@assays$RNA@var.features %in% genes_sub_all)]
dim(objTN2)
  
####WYC blacklist
load("./Additional_data/genes_black_WYC.rda")
sum(objTN2@assays$RNA@var.features %in% genes_black); length(objTN2@assays$RNA@var.features)
objTN2@assays$RNA@var.features = objTN2@assays$RNA@var.features[!(objTN2@assays$RNA@var.features %in% genes_black)]
mito.genes <- rownames(objTN2@assays$RNA)[grep("^MT-",rownames(objTN2@assays$RNA))]
objTN2@assays$RNA@var.features =objTN2@assays$RNA@var.features[!(objTN2@assays$RNA@var.features %in% mito.genes)]
  
dim(objTN2);length(objTN2@assays$RNA@var.features)
  
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
objTN2 <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
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
            "CXCR4","GNB2L1","ATP5L","SUGCT",###DZGC and GC
            "LMO2","GMDS","PRPSAP2","MARCKSL1",###LZGC
            "STMN1","TUBB","HMGB2","TUBA1B","MKI67",
            "JCHAIN","PRDM1","XBP1","MZB1"
            
)

DotPlot(object = objTN2, features =  B_marker,scale = T,group.by = "celltype") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks= hypoxia,labels=  hypoxia)


###############################################################################
#'                          Manuscipt: figure1D&F                            '#
###############################################################################
load("./Additional_data/color.B.Rdata")

object_BCR <- readRDS("../../BCR_data/panB_BCR_processed_data.rds")
dim(object_BCR)
object=object_BCR

###SHM
my_comparisons <- list(c("Blood","Adjacent"),
                       c("Blood","LN_Met"),
                       c("Blood","Cancer"),
                       c("Adjacent","LN_Met"),
                       c("Adjacent","Cancer"),
                       c("LN_Met","Cancer"))
p1=ggplot(data = object@meta.data,mapping = aes(x =type,y =mu_freq)) +
  geom_boxplot(mapping = aes(fill = type),scale = "width",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons)+
  #scale_fill_manual(values=type)+
  labs(title = "IGHV total_Mut_freq")+#facet_wrap(~c_call,scales = "free_y",ncol=4)+
  xlab("Isotype") + ylab("Mutation frequency") +#coord_cartesian(ylim = c(0, 0.2))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p1

####Isotype
Idents(object)<-"c_call2"
levels(object)

object$type=factor(object$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
object$c_call=factor(object$c_call,levels = c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHD","IGHM"))

p2=dittoBarPlot(object, "c_call",group.by="type",retain.factor.levels = T,main = "B",color.panel= color.B$BCR_col)+ylim(0,1);p2
p2
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
objN <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
out.prefix <- "./"
OR.B.list <- do.tissueDist(cellInfo.tb=objN@meta.data,meta.cluster=objN@meta.data$celltype,loc=objN@meta.data$type,
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
#'                          Manuscipt: figureS2H                             '#
###############################################################################
object <- readRDS("../../BCR_data/panB_BCR_processed_data.rds")
new<-object@meta.data
colnames(new)
dplyr::select(new,c("cell_id","c_call","clone_id")) %>% head()
new$clone_id<-as.vector(new$clone_id)
q<-as.data.frame(table(new$clone_id))
new$CloneType_Freq<-NA

q=arrange(q,desc(Freq))

library(future)
plan("multisession", workers =12)
options(future.globals.maxSize = 10000 * 1024^2)

for(i in as.vector(q$Var1)){
  new[new$clone_id==i,]$CloneType_Freq<-q[q$Var1==i,]$Freq
}

p
object@meta.data<-new
FeaturePlot(object,features = "CloneType_Freq",cols = c("grey","red"),min.cutoff = 1,max.cutoff = 5)
object$CloneType_Freq_levels<-NA
object@meta.data[object$CloneType_Freq>0,]$CloneType_Freq_levels<-"1"
object@meta.data[object$CloneType_Freq>1,]$CloneType_Freq_levels<-"2"
object@meta.data[object$CloneType_Freq>2,]$CloneType_Freq_levels<-"3"
object@meta.data[object$CloneType_Freq>3,]$CloneType_Freq_levels<-"4"
object@meta.data[object$CloneType_Freq>4,]$CloneType_Freq_levels<-"5"
object@meta.data[object$CloneType_Freq>5,]$CloneType_Freq_levels<-"6"
object@meta.data[object$CloneType_Freq>6,]$CloneType_Freq_levels<-"7"
object@meta.data[object$CloneType_Freq>7,]$CloneType_Freq_levels<-"8"
object@meta.data[object$CloneType_Freq>8,]$CloneType_Freq_levels<-"9"
object@meta.data[object$CloneType_Freq>9,]$CloneType_Freq_levels<-">9"
object$CloneType_Freq_levels=factor(object$CloneType_Freq_levels,levels = c("1","2","3","4","5","6","7","8","9",">9"))
DimPlot(object,group.by="CloneType_Freq_levels",label=F,pt.size=0.7,raster=F,cols = 
          c("#ffffe5","#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")   )#c("grey","#4c88b7","#c085ba","#ec7a1d","#7952a2","#d90dfd","#db2886","#2705f7","#69a69c","#d32623","#f1bebc","#b19c6d","#e6e81e","#bfa7ce","#7ad0d3")

ggsave("CloneType_Freq_levels_umap.pdf",width = 6,height = 6)
anno <- object@meta.data
p1 <- ggplot(anno,aes(x=celltype_1229,fill=CloneType_Freq_levels))+
  scale_fill_manual(values=c("#ffffe5","#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")   )+
  geom_bar(position='fill')+theme(axis.text.x = element_text(angle=90))
p1
p2 <- ggplot(anno,aes(x=celltype_1229,fill=CloneType_Freq_levels))+
  scale_fill_manual(values=c("#ffffe5","#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")   )+
  geom_bar()+theme(axis.text.x = element_text(angle=90))
p2/p1
ggsave("./FigS2H.pdf",width = 10,height = 15)

###############################################################################
#'                          Manuscipt: figureS3A                             '#
###############################################################################
data <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
dir.metaInfo <- "./dataFreq/"
dir.create((dir.metaInfo),F,T)

dir.metaInfo <- "./dataFreq/"
dir.create((dir.metaInfo),F,T)

meta.tb1=objN@meta.data
meta.tb1$stype=ifelse(meta.tb1$celltype %in% c("B.14.Plasmablast","B.15.Plasma cell"),"ASCs","B")
meta.tb1$usedForFreq="Y"
meta.tb1=data.table(meta.tb1)
#meta.tb2=meta.tb1[meta.tb1$loc %in% "Tumor",]

if(T){##### calculate frequency using only baseline samples
  {
    
    getFreqTable <- function(tb,stype,type.return="mtx",group.var="meta.cluster")
    {
      out.tb <- as.data.table(ldply(c("Adjacent","Blood",   "Cancer",   "LN_Met"),function(x){
        x.tb <- plotDistFromCellInfoTable(tb[type==x,], plot.type="none",
                                          cmp.var="cancer",min.NTotal=30,
                                          group.var=group.var,donor.var="patient")
        x.tb[,loc:=x]
        x.tb[,stype:=stype]
      }))
      if(type.return=="tb"){
        return(out.tb)
      }else if(type.return=="mtx"){
        d.tb <- out.tb
        ht.tb <- dcast(d.tb,group.var~loc+donor.var,value.var="freq",fill=0)
        ht.mtx <- as.matrix(ht.tb[,-1])
        rownames(ht.mtx) <- ht.tb[[1]]
        print(ht.mtx[,1:3])
        return(ht.mtx)
      }
    }
    
    
    
    freq.B.ht.tb <- getFreqTable(meta.tb1[usedForFreq=="Y",],
                                 stype="B" ,type.return="tb",group.var = "celltype")
    
    freq.ASC.ht.tb <- getFreqTable(meta.tb1[usedForFreq=="Y",],
                                   stype="ASCs" ,type.return="tb",group.var = "celltype")
    
    
    freq.all.ht.tb <- rbind(freq.B.ht.tb,freq.ASC.ht.tb)
    saveRDS(freq.all.ht.tb,file=sprintf("%s/panC.freq.all.ht.tb.rds",dir.metaInfo))
    
  }
}

diversity.norm <- function(x)
{
  -sum(ifelse(x>0,x*log2(x),0))/log2(length(x))
}

dat.plot=freq.all.ht.tb
dat.plot$group="B"
dat.plot$group=ifelse(dat.plot$group.var %in% c("B.14.Plasmablast","B.15.Plasma cell"),"PB","B")#,   

dat.plot$type=factor(dat.plot$loc,levels = c("Blood","Adjacent","LN_Met","Cancer"))
dat.diversity.norm.perSample.tb <- dat.plot[,.(NMCls=.N,
                                               NTotal=NTotal[1],
                                               diversity=diversity.norm(freq)),
                                            by=c("group","cmp.var","donor.var","type")]


type=c( "Blood"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")

my_comparisons <- list(c("Blood","Adjacent"),
                       c("Blood","LN_Met"),
                       c("Blood","Cancer"),
                       c("Adjacent","LN_Met"),
                       c("Adjacent","Cancer"),
                       c("LN_Met","Cancer"))
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
obj_bcr <- readRDS("../../BCR_data/panB_BCR_processed_data.rds")
####Diversity analysis###########
Idents(obj_bcr)="celltype";table(Idents(obj_bcr))
obj_ASC <- subset(obj_bcr,idents=c("B.15.Plasma cell","B.14.Plasmablast"))
table(obj_ASC$celltype);table(obj_ASC$type)

new<-obj_ASC@meta.data
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
curve <- estimateAbundance(top150, group="type", ci=0.95, nboot=1000, clone="clone_id",min_n = 20)
# Plots a rank abundance curve of the relative clonal abundances
type=c( "Blood"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
class=c("GC"="#377eb8" ,"EF"="#e41a1c")
plot(curve, colors = type, legend_title="Type diversity")

ggsave("Top150 type abundance.pdf",width = 4,height = 5)

sample_curve <- alphaDiversity(top150, group="type", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,min_n=10,
                               ci=0.95, nboot=1000)

sample_main <- paste0("Type diversity")
type=c( "Blood"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
plot(sample_curve, colors=type, main_title=sample_main, 
     legend_title="Type diversity")


ggsave("Top150 type diversity.pdf",width = 4,height = 5)


###############################################################################
#'                          Manuscipt: figureS3C                             '#
###############################################################################
object <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
anno <- object@meta.data
rm(object)
gc()

PC <- readRDS("../../scRNA_data/panB_Plasma_cell_selected_scRNA_processed_data.rds")
cellname.PC <- colnames(PC)
rm(PC)
gc()

anno.PC <- anno[anno$celltype%in%c("B.15.Plasma cell"),]
anno.PC <- anno.PC[cellname.PC,]
anno.other <- anno[anno$celltype!=c("B.15.Plasma cell"),]
anno.use <- rbind(anno.PC,anno.other)

#pb
use <- anno.use[,c("type","celltype","patient")]
use <- as.data.frame(table(use$patient,use$type,use$celltype))
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
      freq <- temp2$Freq[temp2$Var3=="B.14.Plasmablast"]/(sum(temp2$Freq))
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
  ggtitle("B.14.Plasmablast")+
  #geom_line(aes(group=patient),color='gray',lwd=0.1)+
  stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Blood","Adjacent"),
                                                                c("Blood","LN_Met"),
                                                                c("Blood","Cancer"),
                                                                c("Adjacent","LN_Met"),
                                                                c("Adjacent","Cancer"),
                                                                c("LN_Met","Cancer")))
p1


#ASC
use <- anno.use[,c("type","celltype","patient")]
use <- as.data.frame(table(use$patient,use$type,use$celltype))
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
  ggtitle("B.15.Plasma cell")+
  #geom_line(aes(group=patient),color='gray',lwd=0.1)+
  stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Blood","Adjacent"),
                                                                c("Blood","LN_Met"),
                                                                c("Blood","Cancer"),
                                                                c("Adjacent","LN_Met"),
                                                                c("Adjacent","Cancer"),
                                                                c("LN_Met","Cancer")))
p2


p1|p2
ggsave("FigS3C.pdf",width = 10,height = 5)

###############################################################################
#'                          Manuscipt: figureS3E                             '#
###############################################################################
set.seed(123)
objN <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
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
rogue.res <- rogue(obj.2_expr, labels = metadata$celltype,
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
rogue.res2$variable <- factor(rogue.res2$variable,levels = c("B.01.TCL1A+naiveB",
                                                             "B.02.IFIT3+B","B.03.HSP+B",
                                                             "B.04.MT1X+B","B.05.EGR1+ACB",
                                                             "B.06.NR4A2+ACB2","B.07.CCR7+ACB3",
                                                             "B.08.ITGB1+SwBm","B.09.DUSP4+AtM",
                                                             "B.10.ENO1+Pre_GCB","B.11.SUGCT+DZ_GCB",
                                                             "B.12.LMO2+LZ_GCB","B.13.Cycling_GCB",
                                                             "B.14.Plasmablast","B.15.Plasma cell"))


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
objN_PC <- readRDS("../../scRNA_data/panB_Plasma_cell_selected_scRNA_processed_data.rds")

Idents(objN_PC) <- objN_PC$celltype_l4
allObj <- list(PC = objN_PC)

Idents(objN_PC)="celltype_l4";table(Idents(objN_PC))
figurePath <- c("./")
for(oi in 1:length(allObj)){
  tempName <- names(allObj[oi])
  tempObj <- allObj[[oi]]
  Idents(tempObj) <- tempObj$celltype_l4
  coord = Embeddings(object = tempObj, reduction = "umap")
  coord = coord[,c(1,2)]
  colnames(coord) = c("UMAP_1", "UMAP_2")
  coord = data.frame(ID = rownames(coord), coord)
  meta = tempObj@meta.data
  meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
  meta = left_join(meta, coord, by = 'ID')
  # randomly sample cells of large tissuetype ###############################
  print(tempName)
  print(table(meta$type))
  minNum <- min(table(meta$type))
  print("minNum")
  print(minNum)
  meta <- meta %>%
    group_by(type) %>%
    slice_sample(n = minNum)
  for(ti in 1: length(unique(meta$type))){
    tin <- sort(unique(meta$type))[ti]
    meta_tin <- meta[meta$type == tin,]
    meta_ntin <- meta[meta$type != tin,]
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
data_raw <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
cellname <- sample(colnames(data_raw),50000)
data.sample <- data_raw[,cellname]
DimPlot(data.sample,label=T,group.by = "celltype")
Idents(data.sample) <- data.sample$celltype

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

cds <- order_cells(cds)

p1 <- plot_cells(cds,color_cells_by = "pseudotime",label_branch_points = FALSE,label_leaves = F)
p1
p2 <- plot_cells(cds,color_cells_by = "celltype",
                 label_branch_points = FALSE,label_leaves = F)+
  scale_color_manual(values=my36colors)
p2|p1

ggsave("monocle3_FigS3G.pdf",width = 20,height = 8)


###############################################################################
#'                         Manuscipt: figureS3H                              '#
###############################################################################

DimPlot_theme<-theme_bw()+theme(aspect.ratio=0.4, panel.grid.minor = element_blank(), panel.grid.major = element_blank())

set.seed(123)
object <- readRDS("../../BCR_data/panB_BCR_processed_data.rds")
object$celltype <- factor(object$celltype,
                               levels = unique(object$celltype)[order(unique(object$celltype))])
type <- unique(object$type)
result <- list()
for(k in 1:length(type)){
  objtest=object[,object$type %in% type[k]]
  if(T){
    objtest$celltype_clone_id<-paste0(objtest$celltype,"__",objtest$clone_id)
    objtest$count<-1
    data2<-objtest@meta.data %>% select(c("celltype_clone_id","celltype","clone_id","count"))%>%
      group_by(celltype,clone_id)%>%summarise_at(c("count"),funs(sum))  %>% as.data.frame()
    data2$celltype_clone_id<-paste0(data2$celltype,"__",data2$clone_id)
    Idents(objtest)<-"celltype"
    x<-data.frame(row.names = levels(objtest))
    i=1
    for(i in 1:length(rownames(x))){
      label<-rownames(x)[i]
      tmp<-c()
      for(j in rownames(x)){
        table1<-data2[data2$celltype==label,]
        table2<-data2[data2$celltype==j,]
        rownames(table1)<-table1$clone_id
        rownames(table2)<-table2$clone_id
        clone_id_list<-intersect(table1$clone_id,table2$clone_id)
        tmp<-c(tmp,c(sum(table1[clone_id_list,]$count)+sum(table2[clone_id_list,]$count))/c(sum(table1$count)+sum(table2$count)))
      }
      x[,i]<-tmp
    }
    
    
    label<-rownames(x)[1]
    tmp<-c()
    for(j in rownames(x)){
      table1<-data2[data2$celltype==label,]
      table2<-data2[data2$celltype==j,]
      rownames(table1)<-table1$clone_id
      rownames(table2)<-table2$clone_id
      clone_id_list<-intersect(table1$clone_id,table2$clone_id)
      tmp<-c(tmp,c(sum(table1[clone_id_list,]$count)+sum(table2[clone_id_list,]$count))/c(sum(table1$count)+sum(table2$count)))
    }
    x[,i]<-tmp
    
    colnames(x)<-rownames(x)
    get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
    }
    get_lower_tri<-function(cormat){
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
    }
    x <- get_lower_tri(x)
    x<-x[dim(x)[1]:1,]
    library(data.table)
    melted_cormat <- melt(as.matrix(x), na.rm = TRUE)
    melted_cormat$value<-melted_cormat$value
    
    melted_cormat<-melted_cormat[melted_cormat$value>=0 & melted_cormat$value!=1,]
  }
  melted_cormat$type <- type[k]
  result[[k]] <- melted_cormat
  names(result)[k] <- type[k]
}


melted_cormat=rbind(result$Cancer,result$Adjacent,result$LN_Met,result$Blood)
table(melted_cormat$type)
melted_cormat$type=factor(melted_cormat$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))

f<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+geom_text(aes(label = round(value,3)),color = "black", size = 3)+RotatedAxis();f
P3<-f+ #scale_fill_gradientn(colours = (cc(100)));P3
  scale_fill_gradient2(low = "white", mid = "#F47E5D", high = "#463873", ##low = "#fee8c8", mid = "#fdbb84", high = "#e34a33", 
                       midpoint = max(melted_cormat$value)/2,
                       limit = c(0,max(melted_cormat$value)),
                       space = "Lab", 
                       name="Jaccard index")+theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), 
                                                              panel.grid.major = element_blank())+RotatedAxis()+
  facet_wrap(~ type, ncol = 4, scales = "free") +
  labs(fill =paste0("Jaccard Index",""));P3

ggsave("jaccard.pdf",width = 35,height = 8)

###############################################################################
#'                         Manuscipt: figureS3K                              '#
###############################################################################


############
# overlap BCR ids
BCR <- readRDS("../../BCR_data/panB_BCR_processed_data.rds")
BCR@meta.data[1:3,]
table(BCR$celltype)

metadata = BCR@meta.data
all.ident = unique(metadata$celltype)
all.type = unique(metadata$type)
all.patient = unique(metadata$patient)

metadata$patient_clonetype = paste(metadata$patient, metadata$clone_id, sep="_")
all.patient_clonetype = unique(metadata$patient_clonetype)
stat.patient_clonetype = data.frame(table(metadata$patient_clonetype))

#

all.pairs = list(c("AtM","PC"), c("SwBm","PC"), c("GCB","PC"))

all.pairs = list(
  c("B.14.Plasmablast","B.15.Plasma cell"),
  c("B.09.DUSP4+AtM","B.15.Plasma cell"),
  c("B.13.Cycling_GCB","B.14.Plasmablast"),
  c("B.12.LMO2+LZ_GCB","B.14.Plasmablast"),
  c("B.11.SUGCT+DZ_GCB","B.14.Plasmablast"),
  c("B.14.Plasmablast","B.09.DUSP4+AtM"),
  
  # c("B.07.CCR7+ACB3","B.09.DUSP4+AtM"),
  # c("B.09.DUSP4+AtM","B.06.NR4A2+ACB2"),
  # c("B.05.EGR1+ACB","B.09.DUSP4+AtM"),
  # c("B.09.DUSP4+AtM","B.03.HSP+B"),
  # c("B.02.IFIT3+B","B.09.DUSP4+AtM"),
  # c("B.09.DUSP4+AtM","B.10.ENO1+Pre_GCB")
)



stat_merge_large = c()
for (j in 1:length(all.pairs)){
  print(j)
  
  metadata_tmp = subset(metadata, metadata$celltype %in% all.pairs[[j]])
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


###############################################################################
#'                         Manuscipt: figureS3J                              '#
###############################################################################
BCR <- readRDS("../../BCR_data/panB_BCR_processed_data.rds")
dim(BCR)
Idents(BCR) <- BCR$celltype
unique(BCR$celltype)
anno <- BCR@meta.data
cloneid <- list()
for(i in 1:length(unique(anno$celltype))){
  temp <- anno[anno$celltype==unique(anno$celltype)[i],]
  temp.cloneid <- unique(temp$clone_id)
  cloneid[[i]] <- temp.cloneid
  names(cloneid)[i] <- unique(anno$celltype)[i]
}


###plot

names(cloneid)
data <- fromList(cloneid)
colnames(data) <- gsub("\\+","_",colnames(data))
colnames(data) <- gsub(" ","_",colnames(data))
colnames(data)
celltype <- colnames(data)
celltype <- celltype[-14]
celltype <- celltype[order(celltype,decreasing = T)]
intersection <- list()
for(i in 1:length(celltype)){
  intersection[[i]] <- c("B.15.Plasma_cell",celltype[i])
}

pdf("./cloneshare_upset.pdf",width = 10,height = 8)
upset(data,nsets = 15,keep.order = TRUE,
      intersections = intersection,
      sets = colnames(data)[order(colnames(data),decreasing = T)],
      mb.ratio = c(0.4, 0.6)
      #group.by = 'sets'
)


dev.off()