##########################
### fig.S1A. B cell consensus score ###
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Matrix.utils)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(cowplot)
library(harmony)
library(SeuratWrappers)


load(file = "../PC_scRNA_B/TCGA/TCGA/tcga_xcell.rda")
load(file = "../PC_scRNA_B//TCGA/TCGA/tcga_mcp.rda")
load(file = "../PC_scRNA_B//TCGA/TCGA/tcga_quantiseq.rda")
load(file = "../PC_scRNA_B//TCGA/TCGA/tcga_signature.rda")
load(file = "../PC_scRNA_B/TCGA/TCGA/tcga_obj.rda")

tcga_xcell$tmp = substr(row.names(tcga_xcell), 1, 15); tcga_xcell<-tcga_xcell[!duplicated(tcga_xcell$tmp), ]; row.names(tcga_xcell) = tcga_xcell$tmp; tcga_xcell$tmp = NULL
tcga_mcp$tmp = substr(row.names(tcga_mcp), 1, 15); tcga_mcp<-tcga_mcp[!duplicated(tcga_mcp$tmp), ]; row.names(tcga_mcp) = tcga_mcp$tmp; tcga_mcp$tmp = NULL
tcga_quantiseq$tmp = substr(row.names(tcga_quantiseq), 1, 15); tcga_quantiseq<-tcga_quantiseq[!duplicated(tcga_quantiseq$tmp), ]; row.names(tcga_quantiseq) = tcga_quantiseq$tmp; tcga_quantiseq$tmp = NULL
tcga_signature=as.data.frame(tcga_signature);rownames(tcga_signature)=tcga_signature$ID
tcga_signature$tmp = substr(row.names(tcga_signature), 1, 15); tcga_signature<-tcga_signature[!duplicated(tcga_signature$tmp), ]; row.names(tcga_signature) = tcga_signature$tmp; tcga_signature$tmp = NULL

sample_intersect = intersect(colnames(tcga_obj), row.names(tcga_xcell))
sample_intersect = intersect(sample_intersect, row.names(tcga_mcp))
sample_intersect = intersect(sample_intersect, row.names(tcga_quantiseq))
sample_intersect = intersect(sample_intersect, row.names(tcga_signature))
length(sample_intersect)

tcga_obj_sub = subset(tcga_obj, cells = sample_intersect)
sample_intersect = colnames(tcga_obj_sub)

tcga_xcell = tcga_xcell[sample_intersect,]
tcga_mcp = tcga_mcp[sample_intersect,]
tcga_quantiseq = tcga_quantiseq[sample_intersect,]
tcga_signature = tcga_signature[sample_intersect,]

tcga_obj_sub$B_xcell = tcga_xcell$B.cells
tcga_obj_sub$plasma_xcell = tcga_xcell$Plasma.cells
tcga_obj_sub$B_mcp = tcga_mcp$B_lineage_MCPcounter
tcga_obj_sub$B_quantiseq = tcga_quantiseq$B_cells_quantiseq
tcga_obj_sub$B_Rooney = tcga_signature$B_cells_Rooney_et_al
tcga_obj_sub$B_Danaher = tcga_signature$B_cells_Danaher_et_al
tcga_obj_sub$B_Bindea = tcga_signature$B_cells_Bindea_et_al

colnames(tcga_obj_sub@meta.data)
tcga_obj_sub$project_id[tcga_obj_sub$project_id == ""] = NA

tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, !is.na(tcga_obj_sub@meta.data$project_id))))
tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, (tcga_obj_sub@meta.data$project_id != "TCGA-LAML"))))
tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, (tcga_obj_sub@meta.data$project_id != "TCGA-DLBC"))))
tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, !(tcga_obj_sub@meta.data$project_id %like% "TARGET"))))

tcga_obj_sub = SCTransform(tcga_obj_sub)
tcga_obj_sub <- RunUMAP(tcga_obj_sub, return.model = TRUE, dims = 1:50)
tcga_obj_sub <- RunTSNE(tcga_obj_sub, return.model = TRUE, dims = 1:50)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
tcga_obj_sub$B_xcell_norm = range01(tcga_obj_sub$B_xcell)
tcga_obj_sub$plasma_xcell_norm = range01(tcga_obj_sub$plasma_xcell)
tcga_obj_sub$B_mcp_norm = range01(tcga_obj_sub$B_mcp)
tcga_obj_sub$B_quantiseq_norm = range01(tcga_obj_sub$B_quantiseq)
tcga_obj_sub$B_Rooney_norm = range01(tcga_obj_sub$B_Rooney)
tcga_obj_sub$B_Danaher_norm = range01(tcga_obj_sub$B_Danaher)
tcga_obj_sub$B_Bindea_norm = range01(tcga_obj_sub$B_Bindea)

tcga_obj_sub$B_consensus = rowMeans(cbind(tcga_obj_sub$B_xcell_norm, tcga_obj_sub$plasma_xcell_norm,tcga_obj_sub$B_mcp_norm, 
                                          tcga_obj_sub$B_quantiseq_norm,tcga_obj_sub$B_Rooney_norm,tcga_obj_sub$B_Danaher_norm,tcga_obj_sub$B_Bindea_norm))
quantile(tcga_obj_sub$B_consensus)

meta=tcga_obj_sub@meta.data
cc <- colorRampPalette(c("#253494", "#2c7fb8", "#41b6c4", "#bae4bc", "#f0f9e8")) #matched 2
meta$project_id = reorder(meta$project_id, meta$B_consensus, median)
p14 = ggplot(meta, aes(x = project_id,  y = B_consensus, fill = project_id)) + 
  geom_boxplot(outlier.shape = NA) + ylim(0,0.7) +scale_fill_manual(values = cc(length(unique(b_comb2$project_id)))) +
  geom_hline(yintercept = c(0.2303,0.2855,0.3485),colour="#990000",linetype="dashed")+
  # xlab("Isotype") + ylab("Mutation frequency") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none");p14
ggsave("./TCGA/B_consensus.pdf",width = 8,height=6)

###########################
##### Fig.1B  UMAP

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Matrix.utils)
library(data.table)
library(scCancer)
library(harmony)
library(tidyverse)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(sctransform)
library(glmGamPoi)
library(viridis)


objN<- NormalizeData(objN)
objN<- FindVariableFeatures(objN, nfeatures = 3000)


genes<-data.frame(data.table::fread("../PC_scRNA_B/ST/biomart_mart_export.txt",header=T,sep="\t"))

table(genes$GeneType)
genes_sub<-subset(genes,GeneType!="lncRNA") #processed_pseudogene lncRNA
genes_sub<-subset(genes_sub,GeneType!="processed_pseudogene") #processed_pseudogene lncRNA
genes_sub_all<-c(as.character(genes_sub$gene),as.character(genes_sub$GeneSynonym))
sum(objN@assays$RNA@var.features %in% genes_sub_all); length(objN@assays$RNA@var.features)
objN@assays$RNA@var.features = objN@assays$RNA@var.features[(objN@assays$RNA@var.features %in% genes_sub_all)]
dim(objN)

load("../PC_scRNA_B/publisheddata/0725/genes_black_WYC.rda")
sum(objN@assays$RNA@var.features %in% genes_black); length(objN@assays$RNA@var.features)
objN@assays$RNA@var.features = objN@assays$RNA@var.features[!(objN@assays$RNA@var.features %in% genes_black)]
dim(objN);length(objN@assays$RNA@var.features)

orig.ident.stat = data.frame(table(objN$orig.ident))
orig.ident.stat_sub = subset(orig.ident.stat, orig.ident.stat$Freq > 10)

objN = subset(objN, cells = row.names(subset(objN@meta.data, objN@meta.data$orig.ident %in% orig.ident.stat_sub$Var1)))
dim(objN)

objN <- ScaleData(objN)
objN <- RunPCA(objN, verbose = FALSE)

objN <- harmony::RunHarmony(objN,"patient", plot_convergence = TRUE)

objN <- Seurat::RunUMAP(objN,reduction = "harmony", dims = 1:10) 
objN <- Seurat::FindNeighbors(objN,reduction = "harmony", dims = 1:10) 
objN <- Seurat::FindClusters(objN,resolution =1.5,graph.name="RNA_snn")

objN.celltypemarkers <- Seurat::FindAllMarkers(objN, only.pos = TRUE, min.pct = 0.25)

p3=DimPlot(objN, reduction = "umap", label = F, pt.size = .000001, group.by = "celltype_l3", cols = my36colors,
           label.size = 3, raster = T) + theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), 
                                                          panel.grid.major = element_blank());p3


########## Fig1 C Dotplot
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)
library(dittoSeq)
library(viridis)
library(harmony)

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

DotPlot(object = objN, features =  B_marker,scale = T,group.by = "celltype_l3") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks= hypoxia,labels=  hypoxia)

########## Fig1 D and F BCR isotype and SHM
library(dittoSeq)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(Seurat)
library(pheatmap)
library(cowplot)
library(harmony)
library(dplyr)
library(alakazam)
library(shazam)
library(readr)
library(tidyverse)
library(ggpubr)
library(scRepertoire)
library(viridis)
library(ggsci)
library(viridis)
library(RColorBrewer)
library(ggpie)

load("../PC_scRNA_B/color.B.Rdata")

db<-read.table("/home/data/vip13t45/project/PC/PC_scRNABNew/PC_B_arranged221031/combined_BCR/mjq/docker_final230201/total_final230201_heavy_germ-pass.tsv",header=T,sep="\t")

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




########## Fig1. E OR

##OR值绘制
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
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

#############fig1D and E###
f5=dittoBarPlot(objN, "celltype_l3",group.by="type",retain.factor.levels = T,main = "type",color.panel= my36colors)+ylim(0,1);f5
f6=dittoBarPlot(objN, "type",group.by="celltype_l3",retain.factor.levels = T,main = "type",color.panel= type)+ylim(0,1);f6
f7=dittoBarPlot(objN, "cancer",group.by="celltype_l3",retain.factor.levels = T,main = "cancer",color.panel= my36colors)+ylim(0,1);f7

#############fig1f###

FeaturePlot(objN,features = c("TCL1A", ####NaiveB
                              "IFIT3",###IFN
                              "DNAJB1", #"HSPA1A",###Activated
                              "MT1X",
                              "EGR1",####ACB1
                              "NR4A2",####ACB2
                              "CCR7",
                              "ITGB1",
                              "DUSP4",
                              "PSME2",###PreGC
                              "CXCR4",###DZGC and GC
                              "LMO2",###LZGC
                              "STMN1",
                              "MZB1"
                              
),ncol = 7,raster=T,cols =  brewer.pal(9, "YlOrRd"),min.cutoff = 1, max.cutoff = 3)
ggsave("./B_annotation/objTN_BN featureplot for figs1.pdf",height = 8,width = 30)



