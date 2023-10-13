library(sscVis)
library(Startrac)
library(tictoc)
library(muscat)
library(SingleCellExperiment)
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(plyr)
library(ggpubr)
library(pheatmap)
library(Seurat)
library(dittoSeq)
library(patchwork)
library(pheatmap)
library(cowplot)
library(harmony)
library(dplyr)
library(alakazam)
library(shazam)
library(readr)
library(tidyverse)
library(scRepertoire)
library(viridis)
library(ggsci)
library(ggpie)
library(ggraph)
library(ROGUE)
library(DESeq2)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

###filter
object <- readRDS("BCR_add_ASC_select_20230919.rds")
anno <- object@meta.data
unique(anno$celltype_level3_0912)
unique(anno$c_call)
anno.use <- anno[anno$celltype_level3_0912=="c01_TCL1A+naiveB"&anno$c_call%in%c("IGHM","IGHD"),]
anno.use2 <- anno[anno$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC")&anno$c_call%in%c("IGHM","IGHA1","IGHG2","IGHG3","IGHA2","IGHG1", "IGHG4"),]
anno.use3 <- anno[anno$celltype_level3_0912%in%c("c08_ITGB1+SwBm" ,"c10_PreGC","c06_NR4A2+ACB",
                                         "c09_AtM", "c05_EGR+ACB","c03_HSP+B","c07_CCR7+ACB",         
                                         "c02_IFIT3+B","c13_Cycling GCB","c12_LZ GCB","c11_DZ GCB","c04_MT1X+B" )&anno$c_call%in%c("IGHM","IGHA1","IGHG2","IGHG3","IGHA2","IGHG1", "IGHG4","IGHD"),]
anno.final <- rbind(anno.use3,anno.use2) 
anno.final <- rbind(anno.final,anno.use)
object <- object[,rownames(anno.final)]
table(object$c_call)
saveRDS(object,file="BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")


###############################################################################
#'                    Calculate pTrans by Startrac                           '#
###############################################################################


object <- readRDS("./BCR_add_ASC_select_20230919.rds")

dir.startrac <- "./review/BCR"
dir.create(dir.startrac,F,T)

Idents(object)="type";table(Idents(object))

BCR.data=as.data.table(object@meta.data)

BCR.data$celltype_level3_0912=as.character(BCR.data$celltype_level3_0912)
#
BCR.data$celltype_l4=ifelse(BCR.data$celltype_level3_0912 %in% c("c15_MZB1+ASC","c14_PB"),
                           "ASC",BCR.data$celltype_level3_0912)
#split ASC & PB
#BCR.data$celltype_l4=BCR.data$celltype_level3_0912

BCR.data$majorCluster <- as.character(BCR.data$celltype_l4)
BCR.data$clone.id <- BCR.data$clone_id


#### filter out patient.cluster with number of cell < 10
ncell.patient.cluster <- sort(unclass(BCR.data[,table(sprintf("%s.%s",patient,majorCluster))]))
BCR.data <- BCR.data[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
table(BCR.data$patient)
dim(BCR.data)

BCR.data$Cell_Name=BCR.data$BCR_id
BCR.data$clone.id=BCR.data$clone_id
BCR.data$patient=BCR.data$patient
BCR.data$majorCluster=BCR.data$celltype_l4
BCR.data$loc=BCR.data$type

##all type
if(!file.exists(sprintf("%s/B.out.startrac_alltype.nperm100_ASCselect_20231004.rds",dir.startrac))){
  tic("Startrac.run")
  out <- Startrac.run(BCR.data, proj="panC",verbose=T,cores=1,n.perm=100)
  toc()
  out2 <- list()
  out2[['proj']] <- out@proj
  out2[['cluster.data']] <- out@cluster.data
  out2[['cluster.sig.data']] <- out@cluster.sig.data
  out2[['pIndex.migr']] <- out@pIndex.migr
  out2[['pIndex.tran']] <- out@pIndex.tran
  out2[['pIndex.sig.migr']] <- out@pIndex.sig.migr
  out2[['pIndex.sig.tran']] <- out@pIndex.sig.tran
  saveRDS(out2,sprintf("%s/B.out.startrac_alltype.nperm100_ASCselect_20231004.rds",dir.startrac))
})

###Single cancer
colnames(BCR.data)
cancerType.vec <- unique(BCR.data$cancer)
cancerType.vec <- BCR.data[,unique(cancer)];cancerType.vec;length(cancerType.vec)
#BCR.data <- BCR.data[BCR.data$type=="Cancer",]
if(!file.exists(sprintf("%s/B.out.startrac.total.ASCselect_20231005.rds",dir.startrac))){
  res.byCancerType <- llply(cancerType.vec,function(x){
    Startrac.run(BCR.data[cancer==x,],
                 proj=x,verbose=F,cores=1,n.perm=NULL)
  })
  names(res.byCancerType) <- cancerType.vec
  saveRDS(res.byCancerType,sprintf("%s/B.out.startrac.total.ASCselect_20231005.rds",dir.startrac))
}

###############################################################################
#'                      Define EF&GC derived ASCs                            '#
###############################################################################
object <- readRDS("./review/BCR/BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
dim(object)
Idents(object) <- object$celltype_level3_0912

####################ef and gc identificationclone type###########
table(object$celltype_level3_0912);table(object$type)

#object=object_save
object$celltype_level3_0912=as.character(object$celltype_level3_0912)
object$celltype_l4=ifelse(object$celltype_level3_0912 %in% c("c14_PB","c15_MZB1+ASC"),"PCs",object$celltype_level3_0912)

object$celltype_l5=ifelse(object$celltype_level3_0912 %in% c("c15_MZB1+ASC"),"ASC",
                          ifelse(object$celltype_level3_0912 %in% c("c14_PB"),"PB",
                          ifelse(object$celltype_level3_0912 %in% c("c13_Cycling GCB","c12_LZ GCB","c11_DZ GCB","c08_ITGB1+SwBm"),"GC", "EF")))

table(object$celltype_l5)

table(object$celltype_l5)
table(object$celltype_level3_0912)
table(object$celltype_l5,object$celltype_level3_0912)


Idents(object)="type";table(Idents(object))
#object=subset(object,idents="Cancer")
object$type_celltype=paste0("All","-",object$celltype_l5);table(object$type_celltype)
#object$type_celltype=paste0("All","-",object$celltype_level3_0912);table(object$type_celltype)
#object$type_celltype=paste0(object$type,"-",object$celltype_l5);table(object$type_celltype)

patient=object$patient %>% unique();head(patient);length(patient)
table(object$celltype_l5)
table(object$type_celltype)

obj_pro=list()
for(i in 1:length(unique(object$patient))){
  a=object[,object$patient %in% patient[i]];table(a$patient)
  B89pro <- table(a$clone_id,a$type_celltype)%>% as.data.frame()
  B89pro=B89pro[B89pro$Freq !=0 ,]
  B89pro$patient=patient[i]
  obj_pro[[i]]=B89pro
}


obj_proN=do.call(rbind,obj_pro);head(obj_proN)
#obj_proN$Var2 <- str_split(obj_proN$Var2,"-",simplify = T)[,2]
obj_proN$type=str_split(obj_proN$Var2,"-",simplify = T)[,1]
obj_proN$celltype=str_split(obj_proN$Var2,"-",simplify = T)[,2];head(obj_proN)
obj_proN <- obj_proN[,-2]


e=spread(obj_proN,key = "celltype",value = "Freq");head(e)
e[is.na(e)]=0;head(e)

efPB=e[e$EF >0  & e$GC<1  & e$PB >0 ,];dim(efPB);table(efPB$patient);head(efPB)
efASC=e[e$EF >0  & e$GC<1  & e$ASC >0,];dim(efASC);table(efASC$patient);head(efASC)
efPBASC=e[e$EF >0  & e$GC<1 & e$PB >0 & e$ASC >0,];dim(efPBASC);table(efPBASC$patient);head(efPBASC)

ef=rbind(efPB,efASC,efPBASC)
table(duplicated(ef$Var1))


gcPB=e[e$GC >0 &  e$EF<1& e$PB >0 ,];dim(gcPB);table(gcPB$patient);head(gcPB)
gcASC=e[e$GC >0 &  e$EF<1 & e$ASC >0,];dim(gcASC);table(gcASC$patient);head(gcASC)
gcPBASC=e[e$GC >0 &  e$EF<1 & e$PB >0 & e$ASC >0,];dim(gcPBASC);table(gcPBASC$patient);head(gcPBASC)


gc=rbind(gcPB,gcASC,gcPBASC)
table(duplicated(gc$Var1))

object$groupn=ifelse(object$clone_id %in% ef$Var1,"EF",ifelse(object$clone_id %in% gc$Var1,"GC","others"));table(object$groupn)
table(object$groupn,object$celltype_l4)
table(object$groupn)
table(object$groupn,object$cancer)
data.plot <- object@meta.data
data.plot <- data.plot[data.plot$groupn%in%c("EF","GC")&
                         data.plot$celltype_l5%in%c("ASC","PB"),]
saveRDS(data.plot,file = "cell_level_EF_GC_anno_1005.rds")

###############################################################################
#'                         Manuscipt: figure2D                               '#
###############################################################################

out <- readRDS("B.out.startrac_alltype.nperm100_ASCselect_20231004.rds")

pIndex.sig.tran <- as.data.table(out$pIndex.sig.tran)

dat.plot <- pIndex.sig.tran
dat.plot[,index:=as.character(index)]
color.B$cluster.name <- color.B$celltype
names(color.B$cluster.name) <- names(color.B$celltype)
f.col <- !is.na(names(color.B$cluster.name))
color.B$cluster.name <- color.B$cluster.name[f.col]

out.prefix <- c("BCR/")

dat.tmp <- dat.plot[majorCluster=="ASC" & aid=="panC" &
                        !is.na(value),]
dat.med.tmp <- dat.tmp[order(value),]
dat.tmp[,index:=factor(index,levels=dat.med.tmp$index)]
dat.tmp[,index.name:=factor(index,levels=dat.med.tmp$index.name)]
p <- ggplot(dat.tmp, aes(index,value)) +
    geom_col(fill="steelblue",col=NA,width=0.8) +
    geom_hline(yintercept=0.02,linetype="dashed",color="black",alpha=0.8) +
    xlab("") + ylab(sprintf("pTrans Of %s","ASC")) +
    theme_pubr() +
    theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1));p
ggsave(sprintf("%s pTrans.Fig.barplot_nperm1000_B14_B15_merge_alltype.%s.pdf",out.prefix,"PCs"),width=5.5,height=4.5)


###############################################################################
#'                           Manuscipt: figure2G                             '#
###############################################################################
readRDS("BCR_add_ASC_select_20230919.rds")
atm=subset(object,celltype_level3_0912 == "B09.DUSP4+AtM");table(atm$celltype_level3_0912)..B15.PC
b=table(atm$clone_id,atm$type) %>% as.data.frame()
b2=reshape2::dcast(b,Var1~Var2);dim(b2)
rownames(b2)=b2$Var1;head(b2)
b2=b2[,-1]

# Select a clone from the example database
object=obj_bcr
data=object@meta.data

a=table(object$clone_id,object$celltype_level3_0912) %>% as.data.frame()
a1=arrange(a,desc(Freq));head(a1,10)
head(a1,70)

data_1=subset(data, clone_id == "195636_135121");dim(data_1);table(data_1 $celltype_level3_0912);table(data_1 $cancer);table(data_1 $patient);table(data_1 $type);table(data_1 $c_call)

if(T){
  clone <- makeChangeoClone(data_1, text_fields=c("c_call","celltype_level3_0912","type"), 
                            num_fields="consensus_count")
  # 1. Run PHYLIP and process output
  phylip_exec <- "~/phylip-3.697/exe/dnapars"
  PhylipLineage <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
  
  data.frame(clone_id=PhylipLineage$clone,
             junction_length=PhylipLineage$junc_len,
             v_gene=PhylipLineage$v_gene,
             j_gene=PhylipLineage$j_gene)
  
  data.frame(sequence_id=V(PhylipLineage)$name, 
             c_call=V(PhylipLineage)$c_call,
             consensus_count=V(PhylipLineage)$consensus_count)
  
  #plot(PhylipLineage)
  
  # Modify PhylipLineage and plot attributes,majorcelltype
  V(PhylipLineage)$color <-  "steelblue"
  V(PhylipLineage)$color[V(PhylipLineage)$name == "Germline"] <- "black"
  V(PhylipLineage)$color[grepl("Inferred", V(PhylipLineage)$name)] <- "white"
  V(PhylipLineage)$label <- V(PhylipLineage)$celltype_level3_0912
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B01.TCL1A+naiveB"] <- '#E5D2DD'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B02.IFIT3+B"] <- '#53A85F'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B03.HSP+B"] <-'#F1BB72'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B05.EGR1+ACB1"] <-'#D6E7A3'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B06.NR4A2+ACB2"] <-'#57C3F3'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B08.ITGB1+SwBm"] <- '#E95C59'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B09.DUSP4+AtM"] <- '#E59CC4'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B11.SUGCT+DZGC"] <- '#23452F'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B12.LMO2+LZGC"] <- '#BD956A'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B13.Cycling GCB"] <- '#8C549C'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B07.CCR7+ACB3"] <-'#476D87'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B10.ENO1+PreGC"] <-'#AB3282'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B14.PB"] <- '#585658'
  V(PhylipLineage)$color[V(PhylipLineage)$label == "B15.PC"] <- '#9FA3A8'
  V(PhylipLineage)$label <- V(PhylipLineage)$type
  V(PhylipLineage)$color[V(PhylipLineage)$label == "Blood"] <- "#4daf4a"
  V(PhylipLineage)$color[V(PhylipLineage)$label == "LN_Met"] <- "#984ea3"
  V(PhylipLineage)$color[V(PhylipLineage)$label == "Adjacent"] <- "#377eb8"
  V(PhylipLineage)$color[V(PhylipLineage)$label == "Cancer"] <- "#e41a1c"
  
  
  # Remove large default margins
  par(mar=c(0, 0, 0, 0) + 0.1)
  # Plot PhylipLineage
  
  plot(PhylipLineage, layout=layout_as_tree, vertex.frame.color="grey", 
       vertex.label.color="black", edge.label.color="black", 
       edge.arrow.mode=0)
}


###############################################################################
#'                             Manuscipt: figure2H                           '#
###############################################################################

dat.startrac <- readRDS("B.out.startrac.total.ASCselect_20231005.rds")
z.max <- 2

####### pIndex (tran)

cancerType.vec <- names(dat.startrac);length(cancerType.vec)
dat.plot.tb <- as.data.table(ldply(cancerType.vec,function(x){
  as.data.table(dat.startrac[[x]]@pIndex.sig.tran)[aid ==x & majorCluster=="ASC",]
}))

dat.plot.dcast.tb <- dcast(data=dat.plot.tb,aid~index,value.var="value",fill=0)

dat.plot.dcast.tb[is.na(dat.plot.dcast.tb)]=0

dat.plot.mtx <- as.matrix(dat.plot.dcast.tb[,-1])
rownames(dat.plot.mtx) <- dat.plot.dcast.tb[[1]]
colnames(dat.plot.mtx)[colMaxs(dat.plot.mtx) < 0.01]
head(dat.plot.mtx)

dat.plot.mtx=dat.plot.mtx[,-c(1)]

dat.plot.mtx.tmp <- t(scale(t(dat.plot.mtx)))
dat.plot.mtx.tmp[ dat.plot.mtx.tmp > z.max ] <- z.max
dat.plot.mtx.tmp[ dat.plot.mtx.tmp < -z.max ] <- -z.max
f.mcls <- colSds(dat.plot.mtx) > 0
colnames(dat.plot.mtx)[f.mcls]
#dat.plot.mtx.tmp[ abs(dat.plot.mtx) < 0.01] <- 0
f.cancerType <- rowSds(dat.plot.mtx.tmp)==0
dat.plot.mtx <- dat.plot.mtx[!f.cancerType,,drop=F]
dat.plot.mtx.tmp <- dat.plot.mtx.tmp[!f.cancerType,,drop=F]

#pTrans.hclust.row <- run.cutree(dat.plot.mtx.tmp[,f.mcls],
#                                k=3,method.distance="cosine",method.hclust="ward.D2")
coul <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(30))

pdf("./review/BCR/pTrans_cancer_1005.pdf",width = 8,height = 8)
pheatmap::pheatmap(dat.plot.mtx.tmp,z.lo=-z.max,z.hi=z.max,z.len=50,
                   color = coul,clustering_method = "ward.D2")
dev.off()

colnames(dat.plot.mtx.tmp)
rownames(dat.plot.mtx.tmp)
dat.plot.mtx.tmp2 <- dat.plot.mtx.tmp[c(4,14,5,13,8,15,1,11,9,3,7,10,12,2,6),
                                      c(13,10,11,7,8,9,2,4,6,5,1,3)]

pdf("./pTrans_cancer_alltype_0921_Fig2H.pdf",width = 8,height = 8)
pheatmap::pheatmap(dat.plot.mtx.tmp2,z.lo=-z.max,z.hi=z.max,z.len=50,
                   color = coul,clustering_method = "ward.D2",
                   cluster_cols = F,cluster_rows = F,border_color = "NA")
dev.off()

###############################################################################
#'         Manuscipt: figureS4B(define EF&GC domaint patient)                '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_20230919.rds")
in.dat=as.data.table(object@meta.data)
in.dat$celltype_l3=as.character(in.dat$celltype_level3_0912)
unique(in.dat$celltype_level3_0912)
in.dat$celltype_l4=ifelse(in.dat$celltype_level3_0912 %in% c("c15_MZB1+ASC","c14_PB"),"ASC",in.dat$celltype_level3_0912)

in.dat$majorCluster <- as.character(in.dat$celltype_l4)
in.dat$clone.id <- in.dat$clone_id


#### filter out patient.cluster with number of cell < 10
ncell.patient.cluster <- sort(unclass(in.dat[,table(sprintf("%s.%s",patient,majorCluster))]))

in.dat <- in.dat[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
####
dim(in.dat)

out <- readRDS("B.out.startrac_alltype.nperm100_ASCselect_20231004.rds")
cex.point <- 0.5

in.dat.flt <- in.dat
dim(in.dat.flt)

ncell.patient.cluster <- sort(unclass(in.dat.flt[,table(sprintf("%s.%s",patient,majorCluster))]))
in.dat.flt <- in.dat.flt[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
dim(in.dat.flt)


patient2cancerType.tb <- in.dat.flt[,.N,by=c("patient","cancer")][order(cancer,patient),]
TexSize.tb <- in.dat.flt[majorCluster=="ASC" ,.(N.Tex=.N),
                         by=c("patient","cancer")][order(cancer,patient),]
patient2cancerType.tb <- merge(patient2cancerType.tb,TexSize.tb,by=c("patient","cancer"))

i="pair1"
out.prefix <- c("review/BCR/")
out.prefix=sprintf("%s",out.prefix)
obj.startrac.out=out$pIndex.sig.tran
mcls.pairs=list("pair1"=c( "c01_TCL1A+naiveB","c02_IFIT3+B",
                                 "c03_HSP+B","c04_MT1X+B",
                                 "c05_EGR+ACB","c06_NR4A2+ACB",
                                 "c07_CCR7+ACB","c08_ITGB1+SwBm",
                                 "c09_AtM" , "c10_PreGC" ,      
                                 "c11_DZ GCB","c12_LZ GCB",
                                 "c13_Cycling GCB" 
                                 
      ))

  i=1
  
    mcls.1 <- mcls.pairs[[i]][1]
    mcls.2 <- mcls.pairs[[i]][2]
    mcls.3 <- mcls.pairs[[i]][3]
    mcls.4 <- mcls.pairs[[i]][4]
    mcls.5 <- mcls.pairs[[i]][5]
    mcls.6 <- mcls.pairs[[i]][6]
    mcls.7 <- mcls.pairs[[i]][7]
    mcls.8 <- mcls.pairs[[i]][8]
    mcls.9 <- mcls.pairs[[i]][9]
    mcls.10 <- mcls.pairs[[i]][10]
    mcls.11 <- mcls.pairs[[i]][11]
    mcls.12 <- mcls.pairs[[i]][12]
    mcls.13 <- mcls.pairs[[i]][13]

    get_density <- function(x, y, ...) {
      dens <- MASS::kde2d(x, y, ...)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }
    
    dat.Tex.pTran.2path.a.tb <- as.data.table(obj.startrac.out)[aid!="panC",
    ][majorCluster=="ASC" &
        index %in% c(mcls.1,mcls.2,mcls.3,mcls.4,mcls.5,mcls.6,mcls.7,mcls.8,mcls.9,mcls.10,mcls.11,mcls.12,mcls.13),]
    dat.Tex.pTran.2path.b.tb <- dcast(data=dat.Tex.pTran.2path.a.tb,
                                      formula=aid+majorCluster~index,value.var="value")
    colnames(dat.Tex.pTran.2path.b.tb)[c(1,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c("patient","mcls.1","mcls.2","mcls.3","mcls.4","mcls.5","mcls.6",
                                                                               "mcls.7","mcls.8","mcls.9","mcls.10","mcls.11","mcls.12","mcls.13")#、
    dat.Tex.pTran.2path.b.tb=as.data.table(dat.Tex.pTran.2path.b.tb)
    dat.Tex.pTran.2path.b.tb[is.na(mcls.1),mcls.1:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.2),mcls.2:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.3),mcls.3:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.4),mcls.4:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.5),mcls.5:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.6),mcls.6:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.7),mcls.7:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.8),mcls.8:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.9),mcls.9:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.10),mcls.10:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.11),mcls.11:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.12),mcls.12:=0]
    dat.Tex.pTran.2path.b.tb[is.na(mcls.13),mcls.13:=0]
    
    dat.Tex.pTran.2path.b.tb <- merge(dat.Tex.pTran.2path.b.tb,
                                      patient2cancerType.tb,by="patient")
    dat.Tex.pTran.2path.b.tb[,size:=(log10(N.Tex))*cex.point]
    
    dat.plot.mat <- as.matrix(dat.Tex.pTran.2path.b.tb[,c("mcls.1","mcls.2","mcls.3","mcls.4","mcls.5","mcls.6",
                                                          "mcls.7","mcls.8","mcls.9","mcls.10","mcls.11","mcls.12","mcls.13"),with=F])###,"mcls.3","mcls.4","mcls.5","mcls.6"
    rownames(dat.plot.mat) <- dat.Tex.pTran.2path.b.tb$patient
    ## pattern sort
    dat.plot.mat.pattern <- as.data.table(dat.plot.mat > 0.01)
    dat.plot.mat.pattern$rid <- rownames(dat.plot.mat)
    dat.plot.mat.pattern$V.mcls.1 <- dat.plot.mat[,"mcls.1"]
    dat.plot.mat.pattern$V.mcls.2 <- dat.plot.mat[,"mcls.2"]
    dat.plot.mat.pattern$V.mcls.3 <- dat.plot.mat[,"mcls.3"]
    dat.plot.mat.pattern$V.mcls.4 <- dat.plot.mat[,"mcls.4"]
    dat.plot.mat.pattern$V.mcls.5 <- dat.plot.mat[,"mcls.5"]
    dat.plot.mat.pattern$V.mcls.6 <- dat.plot.mat[,"mcls.6"]
    dat.plot.mat.pattern$V.mcls.7 <- dat.plot.mat[,"mcls.7"]
    dat.plot.mat.pattern$V.mcls.8 <- dat.plot.mat[,"mcls.8"]
    dat.plot.mat.pattern$V.mcls.9 <- dat.plot.mat[,"mcls.9"]
    dat.plot.mat.pattern$V.mcls.10 <- dat.plot.mat[,"mcls.10"]
    dat.plot.mat.pattern$V.mcls.11<- dat.plot.mat[,"mcls.11"]
    dat.plot.mat.pattern$V.mcls.12 <- dat.plot.mat[,"mcls.12"]
    dat.plot.mat.pattern$V.mcls.13 <- dat.plot.mat[,"mcls.13"]
    dat.plot.mat.pattern <- dat.plot.mat.pattern[order(-mcls.1,-mcls.2,-mcls.3,-mcls.4,-mcls.5,-mcls.6,
                                                       -mcls.7,-mcls.8,-mcls.9,-mcls.10,-mcls.11,-mcls.12,-mcls.13,
                                                       -V.mcls.1,-V.mcls.2),]
    dat.plot.mat.pattern[,Group:=sprintf("%s%s",as.integer(mcls.3),as.integer(mcls.2))]
    dat.plot.mat.pattern[,Group:=factor(Group,levels=c("11","10","01","00"))]
    dat.plot.mat <- dat.plot.mat[dat.plot.mat.pattern$rid,,drop=F]
    colnames(dat.plot.mat) <- c(mcls.1,mcls.2,mcls.3,mcls.4,mcls.5,mcls.6,mcls.7,mcls.8,mcls.9,mcls.10,mcls.11,mcls.12,mcls.13)
    
    dat.fisher <- matrix(dat.plot.mat.pattern[,.N,by=c("Group")][["N"]],nrow=2)
    res.fisher <- fisher.test(dat.fisher)
    
    ann.cancerType <- rowAnnotation(cate=dat.plot.mat.pattern$Group,
                                    cancer=dat.Tex.pTran.2path.b.tb[dat.plot.mat.pattern$rid,][["cancer"]],
                                    border=F)
    

    dat.plot.mat1=as.data.frame(dat.plot.mat)
    dat.plot.mat2=as.matrix(dat.plot.mat)
    
    dat.plot.mat3=dat.plot.mat2[rowSums(dat.plot.mat2) >0,]
    #dat.plot.mat1[,14]= rowsum(dat.plot.mat1[1,])

    pheatmap::pheatmap(dat.plot.mat3,scale = "row",clustering_method = 'ward.D2')

    dat.plot.mat3 <- dat.plot.mat3[,c(13,11,12,8,9,10,3,1,2,4,5,6,7)]

    
    pdf("pTrans_ASC_per_patients_FigS4B.pdf",width = 6,height = 10)
    pheatmap::pheatmap(dat.plot.mat3,scale = "row",
                       clustering_method = 'ward.D2',
                       cluster_cols = F)
    dev.off()

###############################################################################
#'                           Manuscipt: figure2I&J                           '#
###############################################################################
###Fig 2I
object <- readRDS("BCR_add_ASC_select_20230919.rds")
EF <- c("CESC5","OV5","GIST2","THCA15","CESC8","HCC39",
        "RCC9","HCC40","LC12","OV4","COAD14","OV3","CESC4","HNSC2",
        "HCC19","BRCA1","HCC31T","HCC32T","PAAD2","CESC3","RCC3",
        "BLCA3","LC11","RCC5","OV2","HCC41","HCC21","HCC36","HCC24",
        "LC14")
GC <- c("COAD3","THCA1","COAD16","COAD2","THYM1",
        "STAD7","THCA13","STAD9","COAD20","COAD15",
        "STAD11","GBC01","COAD24","THCA10","LC13","LC18",
        "LC10T","LC15","STAD10","STAD5")

Idents(object) <- object$celltype_level3_0912 


anno_all <- object@meta.data
anno.use <- anno_all
anno.use$type2 <- ifelse(anno.use$patient%in%EF,"EF",
                    ifelse(anno.use$patient%in%GC,"GC","others"))
anno.use <- anno.use[anno.use$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC")&
                       anno.use$type2!="others",]

anno.use <- anno.use[,c("type","type2","mu_freq")]
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
col_flg <- c("#4DAF4A","#E41A1C","#377EB8","#984EA3")
anno.use$type3 <- paste0(anno.use$type,"_",anno.use$type2)
table(anno.use$type3)
table(is.na(anno.use))
anno.use$type3 <- factor(anno.use$type3,levels = c("Blood_EF","Blood_GC",
                                                   "Adjacent_EF","Adjacent_GC",
                                                   "Cancer_EF","Cancer_GC",
                                                   "LN_Met_GC"))
anno.use$type <- factor(anno.use$type,levels = c("Blood","Adjacent",
                                                 "Cancer","LN_Met"))
ggplot(anno.use,aes(x=type3,y=mu_freq*100,fill=type))+geom_boxplot()+
  theme_classic()+
  stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Adjacent_EF","Adjacent_GC"),
                                        c("Blood_EF","Blood_GC"),
                                        c("Cancer_EF","Cancer_GC")),
                     method="wilcox.test")+
  ylab("Mutation frquency")+xlab('')+
  scale_fill_manual(values=col_flg)+theme(axis.text.x = element_text(angle=45,hjust=1))
ggsave("SHM_EF_GC_Fig2I.pdf",width = 5,height = 5)

###Fig2J
anno <- object@meta.data
anno <- anno[anno$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC")&
               anno$groupn%in%c("EF","GC"),]
anno.use <- as.data.frame(table(anno$groupn,anno$c_call,anno$type))
anno.use <- anno.use[anno.use$Var2!="",]
anno.use <- anno.use[-which(anno.use$Var2%in%c("IGHD","IGHE")),]
anno.use$type <- paste0(as.character(anno.use$Var3),"_",as.character(anno.use$Var1))
anno.use <- anno.use[anno.use$type!="LN_Met_EF",]
anno.use$type <- factor(anno.use$type,levels = c("Blood_EF","Blood_GC",
                                                 "Adjacent_EF","Adjacent_GC",
                                                 "Cancer_EF","Cancer_GC",
                                                 "LN_Met_GC"))
my36colors <-c("#3d7bb0" ,"#a2c8db" ,"#5aa554", "#afd195", "#d83f36", "#e99997" ,"#f2bc7c", "#e68740"
)
names(my36colors) <- c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4")
ggplot(anno.use,aes(x=type,fill=Var2,y=Freq))+geom_bar(stat='identity',position='fill')+
  theme_classic()+scale_fill_manual(values=my36colors)+
  xlab("")+ylab("Percent of cells")+theme(axis.text.x = element_text(angle=45,hjust=1))
ggsave("Fig2J.pdf",width = 6,height = 4)

###############################################################################
#'                           Manuscipt: figure2K                             '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_20230919.rds")
EF <- c("CESC5","OV5","GIST2","THCA15","CESC8","HCC39",
        "RCC9","HCC40","LC12","OV4","COAD14","OV3","CESC4","HNSC2",
        "HCC19","BRCA1","HCC31T","HCC32T","PAAD2","CESC3","RCC3",
        "BLCA3","LC11","RCC5","OV2","HCC41","HCC21","HCC36","HCC24",
        "LC14")
GC <- c("COAD3","THCA1","COAD16","COAD2","THYM1",
        "STAD7","THCA13","STAD9","COAD20","COAD15",
        "STAD11","GBC01","COAD24","THCA10","LC13","LC18",
        "LC10T","LC15","STAD10","STAD5")

Idents(object) <- object$celltype_level3_0912 
EF <- object@meta.data
EF <- EF[EF$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC")&EF$groupn=="EF",]

unique(EF$type)
set.seed(123)

j=1
EF_Cancer <- EF[EF$type==unique(EF$type)[j],]
data2 <- as.data.frame(table(EF_Cancer$clone_id,EF_Cancer$c_call))
data2 <- data2[data2$Var2!=""&data2$Var2!="IGHD",]
data <- dcast(data2,Var1~Var2)
rownames(data) <- data$Var1
data <- data[,-1]
#data <- data[,-1]
num <- colSums(data)
num_percent <- as.data.frame(num/sum(num))
colnames(num_percent) <- c('clone')
num_percent$name <- rownames(num_percent)
num_percent$type <- 'EF'
data <- data[rowSums(data>0)>1,]

compare <- as.data.frame(combn(colnames(data),2))
sum <- matrix(nrow = length(colnames(compare)),ncol = 3)
for(i in 1:length(colnames(compare))){
  temp <- data[,compare[,i]]
  temp <- temp[rowSums(temp>0)>1,]
  sum[i,1] <- compare[1,i]
  sum[i,2] <- compare[2,i]
  sum[i,3] <- sum(apply(temp, 1, min))
}
sum <- as.data.frame(sum)
colnames(sum) <- c("iso1","iso2","value")
result <- sum
result$value <- as.numeric(result$value)

edge <- result[,c(1,2,3)]
edge$value=as.numeric(edge$value)
colnames(edge) <- c("from","to","value")
nodes <- data.frame(name=unique(union(edge$from,edge$to)))
nodes <- merge(nodes,num_percent,by='name')
mygraph <- tbl_graph(nodes=nodes,edges=edge)
mygraph

my36colors <-c("#3d7bb0" ,"#a2c8db" ,"#5aa554", "#afd195", "#d83f36", "#e99997" ,"#f2bc7c", "#e68740"
)
names(my36colors) <- c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4")
col2 <- my36colors

ggraph(mygraph,layout = "circle")+
  geom_edge_link(aes(width=value),alpha=0.5)+
  geom_node_point(aes(color=name,size=clone))+
  geom_node_text(aes(label=name))+theme_void()+
  scale_color_manual(values=col2)+
  scale_edge_width(range = c(0,2))+
  ggtitle(paste0("AtM_EF_",unique(EF$type)[j]))+
  scale_size_continuous(range = c(5,10))
ggsave("Fig2K_ASC_EF_cancer.pdf",width=4,
       height=3)

#GC
GC <- object@meta.data
GC <- GC[GC$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC")&GC$groupn=="GC",]

j=1
GC_Cancer <- GC[GC$type==unique(GC$type)[j],]
data2 <- as.data.frame(table(GC_Cancer$clone_id,GC_Cancer$c_call))
data2 <- data2[data2$Var2!="",]
data2 <- data2[-which(data2$Var2%in%c("IGHD","IGHE")),]
data <- dcast(data2,Var1~Var2)
rownames(data) <- data$Var1
data <- data[,-1]

num <- colSums(data)
num_percent <- as.data.frame(num/sum(num))
colnames(num_percent) <- c('clone')
num_percent$name <- rownames(num_percent)
num_percent$type <- "GC"
data <- data[rowSums(data>0)>1,]

compare <- as.data.frame(combn(colnames(data),2))
sum <- matrix(nrow = length(colnames(compare)),ncol = 3)
for(i in 1:length(colnames(compare))){
  temp <- data[,compare[,i]]
  temp <- temp[rowSums(temp>0)>1,]
  sum[i,1] <- compare[1,i]
  sum[i,2] <- compare[2,i]
  sum[i,3] <- sum(apply(temp, 1, min))
}
sum <- as.data.frame(sum)
colnames(sum) <- c("iso1","iso2","value")
result <- sum
result$value <- as.numeric(result$value)

edge2 <- result[,c(1,2,3)]
edge2$value=as.numeric(edge2$value)
colnames(edge2) <- c("from","to","value")
nodes2 <- data.frame(name=unique(union(edge2$from,edge2$to)))
nodes2 <- merge(nodes2,num_percent,by='name')
mygraph2 <- tbl_graph(nodes=nodes2,edges=edge2)
mygraph2
my36colors <-c("#3d7bb0" ,"#a2c8db" ,"#5aa554", "#afd195", "#d83f36", "#e99997" ,"#f2bc7c", "#e68740"
)
names(my36colors) <- c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4")
col2 <- my36colors

ggraph(mygraph2,layout = "circle")+
  geom_edge_link(aes(width=value),alpha=0.5)+
  geom_node_point(aes(color=name,size=clone))+
  geom_node_text(aes(label=name))+theme_void()+
  scale_color_manual(values=col2)+
  scale_edge_width(range = c(0,2))+
  ggtitle(paste0("AtM_GC_",unique(GC$type)[j]))+
  scale_size_continuous(range = c(5,10))
ggsave("Fig2K_ASC_GC_Cancer.pdf"),width=4,height=3)

###############################################################################
#'                           Manuscipt: figure2L                             '#
###############################################################################
ASC <- readRDS("ASC_rmtype_20230922.rds")
#monocle3
ASC_obj_mnc3 = ASC

ASC_obj_mnc3$celltype_l2 = ASC_obj_mnc3$celltype_l3
table(ASC_obj_mnc3$celltype_l2)
table(ASC_obj_mnc3$celltype_l2)
DimPlot(ASC_obj_mnc3, reduction = "umap", label = TRUE, pt.size = .1, group.by = "celltype_l2", cols = my36colors, label.size = 0, raster = F) + 
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())

obj.input = ASC_obj_mnc3
cds <- as.cell_data_set(obj.input)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)


plot_cells(cds, color_cells_by = "partition") & theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())


number_tmp<-floor(ncol(obj.input)*0.05)


#select root cell
unique(obj.input$celltype_l2)
cell_tmp = row.names(subset(obj.input@meta.data, obj.input@meta.data$celltype_l2 == "PC.08.ATP5E" ))

num<-5
cell_tmp<-cell_tmp[sample(1:length(cell_tmp), num, replace=F)]
cds <- order_cells(cds, root_cells = cell_tmp)



cc <- colorRampPalette(c("#352a86", "#095cd8", "#46b896", "#e7ba4a", "#f8fa0d"))
cc <- colorRampPalette(c("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
# cc = colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
cc <- colorRampPalette(c("grey90", "#fc58a6"))
cc <- colorRampPalette(c("white", "white", "#F47E5D", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"
cc <- colorRampPalette(c("white", "#F47E5D", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_roots = F,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 1,
           graph_label_size=1.5)+
  scale_color_gradientn(colours = (cc(100))) +
  DimPlot(ASC_obj_mnc3, reduction = "umap", label = TRUE, pt.size = .1, group.by = "celltype_l2", cols = my36colors, label.size = 0, raster = F) & 
  theme_bw() & theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())

########cytotrace
ASC <- readRDS("ASC_rmtype_20230922.rds")
cellname <- sample(colnames(ASC),30000)
cellname <- colnames(ASC)[ASC$groupn%in%c("EF","GC")]
obj.2_expr=ASC[,cellname]@assays$RNA@counts

obj.2_expr <- obj.2_expr[rowSums(obj.2_expr)>0,]  

obj.2_expr <- as.matrix(obj.2_expr)# %>% as.data.frame()
obj.2_expr[1:4,1:2]
obj.2_pheno=ASC[,cellname]$celltype_l4
obj.2_pheno = as.character(obj.2_pheno) 

names(obj.2_pheno) <- rownames(ASC[,cellname]@meta.data)
object2_emb=Embeddings(object = ASC[["umap"]])



Plasma.CytoTRACE <- CytoTRACE(mat = obj.2_expr,ncores = 15,enableFast = FALSE)

plotCytoTRACE(Plasma.CytoTRACE, phenotype = obj.2_pheno, emb = object2_emb,
              outputDir = "./")

###
cellid=names(Plasma.CytoTRACE$CytoTRACErank)

cytotrace=data.frame(type=obj.2_pheno,
                     CytoTRACE=Plasma.CytoTRACE$CytoTRACE)

p2=ggplot(data = cytotrace,mapping = aes(x =type ,y =CytoTRACE)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type))  + 
  geom_boxplot(mapping = aes(fill = type),outlier.shape = NA)+
  scale_fill_manual(values=my36colors)+ labs(title = "CytoTRACE")+#facet_wrap(~cancerType,scales = "free_y",ncol=3)+
  xlab("type") + ylab("CytoTRACE score") +#coord_cartesian(ylim = c(0, 0.25))+
  theme_classic() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p2

library(ggpubr)
unique(cytotrace$type)
cytotrace$type=factor(cytotrace$type,levels = c("PC.02.RGS13","PC.01.NME2",
                                                "PC.07.CD83","PC.06.SLC3A2",
                                                "PC.03.DUSP5","PC.10.SPINK2",
                                                "PC.08.IFI6",
                                                "PC.04.HSPA1A","PC.09.ATP5E",
                                                "PC.05.NEAT1"
                                                ))
ggboxplot(cytotrace,x="type",y="CytoTRACE",
          color = "type",palette = my36colors,
          add = "jitter",add.params = list(size=0.1))+ 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  coord_cartesian(ylim = c(0, 1))+xlab('')

ggsave("CytoTRACE_Boxplot_modify.pdf",width = 10,height = 6)




###############################################################################
#'                           Manuscipt: figureS5A                            '#
###############################################################################

object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
anno.use <- object@meta.data
anno.use <- anno.use[,c("groupn","celltype_level3_0912","patient","type")]
anno.use <- anno.use[anno.use$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC"),]
anno.use <- anno.use[anno.use$celltype_level3_0912%in%c("c15_MZB1+ASC"),]
anno.use$cellname <- rownames(anno.use)
ASC.use <- readRDS("cell_level_EF_GC_anno_1005.rds")

colnames(ASC.use)

use <- ASC.use[,c("groupn","celltype_level3_0912","patient")]
use <- use[use$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC"),]
use <- use[use$celltype_level3_0912%in%c("c15_MZB1+ASC"),]
use$cellname <- rownames(use)

library(dplyr)
use <- left_join(anno.use,use,by='cellname')
use <- use[,c(2,3,4,6)]
colnames(use) <- c("celltype_level3_0912","patient","type","groupn")
use$groupn <- ifelse(is.na(use$groupn),"others",use$groupn)
use <- as.data.frame(table(use$patient,use$groupn,use$type))
use$Freq <- as.integer(use$Freq)


plot <- list()
for(j in 1:length(unique(use$Var3))){
  temp2 <- use[use$Var3==unique(use$Var3)[j],]
  use2 <- matrix(nrow = length(unique(temp2$Var1)),ncol=3)
  colnames(use2) <- c("patient","EF","GC")
  use2 <- as.data.frame(use2)
  for(i in 1:length(unique(temp2$Var1))){
    temp <- temp2[temp2$Var1==unique(temp2$Var1)[i],]
    use2[i,1] <- as.character(unique(temp$Var1))
    use2[i,2] <- as.numeric(temp$Freq[1]/sum(temp$Freq))
    use2[i,3] <- as.numeric(temp$Freq[2]/sum(temp$Freq))
  }
  use2 <- use2[use2$EF!="NaN"&use2$GC!="NaN",]
  use2 <- use2[use2$EF!=0|use2$GC!=0,]
  use2.melt <- melt(use2)
  col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
  p <- ggplot(use2.melt,aes(x=variable,y=value,fill=variable))+geom_boxplot()+
    theme_classic()+xlab("")+ylab("Frequency")+
    scale_fill_manual(values=col_flg)+
    stat_compare_means(aes(label=..p.signif..),comparisons = list(c("EF","GC")),method = 't.test')+
    #stat_compare_means(comparisons = list(c("EF","GC")),method = 't.test')+
    ggtitle(unique(use$Var3)[j])
  p
  plot[[j]] <- p
}
plot[['nrow']] <- 2
plot[['ncol']] <- 2
library(gridExtra)
pdf('FigS5A_4type_20231008.pdf', width = 10, height = 10)
do.call('grid.arrange', plot)
dev.off()

###############################################################################
#'                           Manuscipt: figureS5B                            '#
###############################################################################
objN <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
anno <- objN@meta.data
anno <- anno[anno$groupn%in%c("EF","GC")&anno$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC"),]
obj_random2w <- objN[,rownames(anno)]
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
rogue.res <- rogue(obj.2_expr, labels = metadata$groupn,
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


class=c("GC"="#377eb8" ,"EF"="#e41a1c")
my_comparisons=list(c("EF","GC"))
unique(rogue.res2$variable)
p=ggboxplot(rogue.res2,x="variable",y="ROGUE",
            color = "variable",palette = class,
            add = "jitter")+ theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ 
  stat_compare_means(aes(label=..p.signif..),comparisons = my_comparisons)+
  coord_cartesian(ylim = c(0, 1))+xlab("");p

ggsave("FigS5B.pdf",width = 4,height = 5)



###############################################################################
#'                       Manuscipt: figureS5C                                '#
###############################################################################
class=c("GC"="#377eb8" ,"EF"="#e41a1c")
####Diversity analysis###########
objectnPC <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_20230921.rds")
anno <- readRDS("cell_level_EF_GC_anno.rds")
Idents(objectnPC) <- objectnPC$celltype_level3_0912
objectnPC <- subset(objectnPC,idents=c("c14_PB","c15_MZB1+ASC"))
anno <- anno[anno$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC"),]
anno <- anno[anno$groupn%in%c("EF","GC"),]
length(intersect(rownames(anno),colnames(objectnPC)))
objectnPC <- objectnPC[,rownames(anno)]
objectnPC$groupn <- anno$groupn  

Idents(objectnPC) <- objectnPC$groupn
objectnPC <- subset(objectnPC,idents=c("EF","GC"))
table(objectnPC$groupn)
new<-objectnPC@meta.data
colnames(new)
select(new,c("cell_id","c_call","clone_id")) %>% head()
new$clone_id<-as.vector(new$clone_id)
q<-as.data.frame(table(new$clone_id))
q<-q[order(q$Freq,decreasing = T),]
# get top 10 genes
clone_id_150=head(q,n=150)$Var1

top200=new[new$clone_id %in% clone_id_150,]

clones <- countClones(objectnPC@meta.data, group="groupn")
head(clones, 5)
EF <- clones[clones$groupn=="EF",]
EF_top5 <- sum(EF$seq_count[1:5])
EF_all <- sum(EF$seq_count)
EF_top5/EF_all

GC <- clones[clones$groupn=="GC",]
GC_top5 <- sum(GC$seq_count[1:5])
GC_all <- sum(GC$seq_count)
GC_top5/GC_all
# Partitions the data on the sample column
# Calculates a 95% confidence interval via 200 bootstrap realizations
curve <- estimateAbundance(top200, group="groupn", ci=0.95, nboot=100, clone="clone_id")
# Plots a rank abundance curve of the relative clonal abundances
type=c( "Cancer_PBMC"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")

plot(curve, colors = class, legend_title="Type diversity")
ggsave("Top150 EF vs GC objectPC clonal Abundance_20230927.pdf",width = 6,height = 6)


sample_curve <- alphaDiversity(top200, group="groupn", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)

sample_main <- paste0("Type diversity")
type=c( "Cancer_PBMC"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
plot(sample_curve, colors=class, main_title=sample_main, 
     legend_title="Type diversity")

ggsave("Top150 EF vs GC objectPC diversity_20230927.pdf",width = 6,height = 6)




###############################################################################
#'                           Manuscipt: figureS5D                            '#
###############################################################################

object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
library(ggpubr)
anno <- object@meta.data
table(anno$patient,anno$groupn)
anno <- anno[anno$celltype_level3_0912%in%c("c14_PB","c15_MZB1+ASC")&
               anno$groupn%in%c("EF","GC"),]
colnames(anno)
anno.use <- anno[,c("groupn","mu_freq","c_call","type")]
#anno.use <- anno.use[-which(anno.use$c_call%in%c("IGHD","IGHE")),]
#anno.use <- anno.use[anno.use$c_call!="",]
anno.use$type2 <- paste0(anno.use$type,"_",anno.use$groupn)
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
col_flg <- c("#4DAF4A","#E41A1C","#377EB8","#984EA3")
plot <- list()
anno.use$c_call <- factor(anno.use$c_call,levels = c("IGHA1","IGHA2",
                                                     "IGHG1","IGHG2",
                                                     "IGHG3","IGHG4",
                                                     "IGHM"))

for(i in 1:length(levels(unique(anno.use$c_call)))){
  temp <- anno.use[anno.use$c_call==levels(unique(anno.use$c_call))[i],]
  temp$type2 <- factor(temp$type2,levels = c("Blood_EF","Blood_GC",
                                              "Adjacent_EF","Adjacent_GC",
                                              "Cancer_EF","Cancer_GC",
                                              "LN_Met_GC"))
  if(i==6){
    p1 <- ggplot(temp,aes(x=type2,y=mu_freq,fill=type))+
      geom_boxplot()+theme_classic()+ggtitle(levels(unique(anno.use$c_call))[i])+
      theme(axis.text.x = element_text(angle=45,hjust=1))+
      xlab("")+ylab("Mutation frequency")+
      scale_fill_manual(values = col_flg)+
      stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Adjacent_EF","Adjacent_GC"),
                                                                    c("Cancer_EF","Cancer_GC")),
                         method="wilcox.test")
  }
  else{
    p1 <- ggplot(temp,aes(x=type2,y=mu_freq,fill=type))+
      geom_boxplot()+theme_classic()+ggtitle(levels(unique(anno.use$c_call))[i])+
      theme(axis.text.x = element_text(angle=45,hjust=1))+
      xlab("")+ylab("Mutation frequency")+
      scale_fill_manual(values = col_flg)+
      stat_compare_means(aes(label=..p.signif..),comparisons = list(c("Adjacent_EF","Adjacent_GC"),
                                                                    c("Blood_EF","Blood_GC"),
                                                                    c("Cancer_EF","Cancer_GC")),
                         method="wilcox.test")
  }

  plot[[i]] <- p1
}

plot[['nrow']] <- 2
plot[['ncol']] <- 4
library(gridExtra)
pdf('FigS5D.pdf', width = 20, height = 8)
do.call('grid.arrange', plot)
dev.off()

###############################################################################
#'                         Manuscipt: figureS5F&G                            '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
Idents(object) <- object$celltype_level3_0912
cellname.PB <- 

object <- subset(object,idents=c("c14_PB","c15_MZB1+ASC"))
ASC <- readRDS("ASC_rmtype_20230922.rds")
DimPlot(ASC,label = T)
dim(ASC)

ASC.use <- object
dim(ASC.use)
table(ASC.use$patient)
ASC.use_sce <- as.SingleCellExperiment(ASC.use) 
ASC.use_sce$group_id <- ASC.use_sce$groupn
ASC.use_sce$sample_id <- paste0(ASC.use_sce$patient,"_",ASC.use_sce$groupn)
ASC.use_sce$cluster_id <- "ASC.use"
ASC.use_sce <- prepSCE(ASC.use_sce,
                       kid = "cluster_id",
                       sid = "sample_id",
                       gid = "group_id",
                       drop = FALSE)
nk <- length(kids <- levels(ASC.use_sce$cluster_id))
ns <- length(sids <- levels(ASC.use_sce$sample_id))
names(kids) <- kids; names(sids) <- sids
table(ASC.use_sce$group_id,ASC.use_sce$sample_id)
pb <- aggregateData(ASC.use_sce,
                    assay="counts",
                    fun="sum",
                    by=c('cluster_id','sample_id'))
assayNames(pb)
t(head(assay(pb)))
metadata(pb)
#pb_mds <- pbMDS(pb)
ei <- metadata(pb)$experiment_info
mm <- model.matrix(~0+ei$group_id)
dimnames(mm) <- list(ei$sample_id,levels(ei$group_id))
contrast <- makeContrasts("EF-GC",levels=mm)
res <- pbDS(pb,method = "DESeq2",min_cells = 10,design = mm,contrast = contrast)

tmp <- ASC.use_sce
counts(tmp) <- as.matrix(counts(tmp))
result_table <- resDS(tmp,res,bind="row",frq=FALSE,spm=FALSE)
result_table$compare <- ifelse(result_table$logFC>0,"EF","GC")
result_table.use <- result_table[result_table$p_adj.loc<0.05,]
write.csv(result_table, file = "muscat_result_table_EF_vs_GC_min10_patient_20230921.csv", row.names = F)

######volcano
aa <- result_table

logFC <-aa$logFC
logFC <- as.numeric(logFC)
adj <- aa$p_adj.loc
gene<- aa$gene

data <- data.frame(logFC=logFC,padj=adj,gene=gene)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC <0.25)& data$logFC > -0.25] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.25] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.25] <- "down"

# 
x_lim <- max(logFC,-logFC)
# 
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
#pdf(file = "./pdf/volcano_MDDvsNormal.pdf",width=8,height=8)
theme_set(theme_classic())
data$gene<-as.character(data$gene)
data$type <- data$sig
data$label=ifelse((data$logFC>1.5|data$logFC<(-1.5))&data$padj<0.05,data$gene,"")

data$type <- as.character(data$type)
data$logFC <- as.numeric(data$logFC)

p <- ggplot(data,aes(round(logFC,4),-1*log10(padj),
                     color = type))+geom_point()+
  xlim(-10,10) +  labs(x="log2(FoldChange)",y="-log10(pvalue_adj)")
p
p <- p + scale_color_manual(values =c("#BC3C28","grey","#BC3C28"))+ #0072B5
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.25,0.25),linetype=4)

p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p  +guides(colour = FALSE)

p<-p+geom_text_repel(data = data, aes(x = data$logFC, 
                                      y = -log10(data$padj), 
                                      label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE,max.overlaps = 20)
p

ggsave("FigS5F_volcano_muscat_EF_GC_strict_patient_20230921.pdf",width = 10,height = 8)


# GSVA --------------------------------------------------------------------

library(GSVA)
library(msigdbr)
library(GSEABase)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- genesets[grep("^HALLMARK_",genesets$gs_name),]

head(genesets)

expr <- ASC.use@assays$RNA@counts
dim(expr)
rownames(expr)[1:3]
gsva.res<- gsva(expr,gset.idx.list=genesets, method="gsva",parallel.sz=16) 
saveRDS(gsva.res, ".gsva_KEGG_ASC_EF_GC_patient_20230921.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_KEGG_ASC_EF_GC_patient_20230921.csv", row.names = F)
gsva.df[1:3,1:3]

gsva=gsva.df


ac <- data.frame(group=ASC.use$groupn) 
design <- model.matrix(~ 0 + factor(ac$group))
colnames(design) <- levels(factor(ac$group))
# rownames(design) <- colnames(sub_regulonAUC)
head(design)

contrast.matrix <- makeContrasts(EF-GC, levels = design)

fit <- lmFit(gsva.res, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
head(df)

cutoff <- 0
df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf),labels = c(1,2))
head(df)

sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
top2=rbind(head(sortdf,n=9),tail(sortdf,n=9))

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
            hjust = "outward" ) +  
  geom_text(data = subset(top2, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = group),
            size = 4, hjust = "inward") +  
  scale_colour_manual(values = c("black","black"), guide = FALSE) +
  xlab("") +ylab("t value of GSVA score")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank());p1 

p1

ggsave("FigS5G_gsva_HALLMARK_20230921.pdf", width = 10, height = 8)


###############################################################################
#'                         Manuscipt: figureS5I                              '#
###############################################################################
ASC.use <- readRDS("ASC_rmtype_20230922.rds")
anno <- ASC.use@meta.data
anno <- anno[anno$groupn!="others",]
table(anno$groupn)
plot.data <- as.data.frame(table(anno$celltype_l4,anno$groupn,anno$type))
plot.data$type2 <- paste0(plot.data$Var3,"_",plot.data$Var2)
plot.data$Var1 <- gsub("\\.",'_',plot.data$Var1)
plot.data <- plot.data[plot.data$Freq!=0,]
#plot.data <- plot.data[-which(plot.data$Var3%in%c("Lymph","LNN","Thymus")),]
plot.data$type2 <- factor(plot.data$type2,levels = c("Blood_EF","Blood_GC",
                                                 "Adjacent_EF","Adjacent_GC",
                                                 "Cancer_EF","Cancer_GC",
                                                 "LN_Met_GC"))
my12colors<-c("#9ED9EA","#B1D896","#DAC0DD","#FDC37E","#F7A7A8","#008FC5","#ABA9AA","#06AC4B","#F68C1F","#D897C2","#F9AC8A","#42BC99")
ggplot(plot.data,aes(x=type2,y=Freq,fill=Var1))+
  geom_bar(stat='identity',position='fill')+
  theme_classic()+scale_fill_manual(values=my12colors)+
  xlab("")+ylab("Percent of cells")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
ggsave("FigS5I.pdf",width = 6,height = 4)

###############################################################################
#'                       Manuscipt: figureS5J&K                              '#
###############################################################################
DimPlot_theme<-theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())

object<-readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
object$sub_cell_type2<-object$celltype_level3_0912

object$v_call_genotyped<-object$v_call
object$sample_id<-object$patient
object$duplicate_count<-1
# colnames(db)
# colnames(object@meta.data)[match(colnames(db),colnames(object@meta.data))]
# head(db$duplicate_count)
Idents(object)<-"sub_cell_type2"
levels(object)<-levels(object)[order(levels(object))]
object$sub_cell_type2<-Idents(object)

Cell_list1<-rownames(object@meta.data)[object$groupn=="EF" & c(object$sub_cell_type2 %in% c("c14_PB","c15_MZB1+ASC")) &  object$c_call=="IGHM"]
Cell_list2<-rownames(object@meta.data)[object$groupn=="EF" & c(object$sub_cell_type2 %in% c("c14_PB","c15_MZB1+ASC")) &  c(object$c_call %in% c("IGHG1","IGHG2","IGHG3","IGHG4"))]
Cell_list3<-rownames(object@meta.data)[object$groupn=="EF" & c(object$sub_cell_type2 %in% c("c14_PB","c15_MZB1+ASC")) &  c(object$c_call %in% c("IGHA1","IGHA2"))]

Cell_list4<-rownames(object@meta.data)[object$groupn=="GC" & c(object$sub_cell_type2 %in% c("c14_PB","c15_MZB1+ASC")) &  object$c_call=="IGHM"]
Cell_list5<-rownames(object@meta.data)[object$groupn=="GC" & c(object$sub_cell_type2 %in% c("c14_PB","c15_MZB1+ASC")) &  c(object$c_call %in% c("IGHG1","IGHG2","IGHG3","IGHG4"))]
Cell_list6<-rownames(object@meta.data)[object$groupn=="GC" & c(object$sub_cell_type2 %in% c("c14_PB","c15_MZB1+ASC")) &  c(object$c_call %in% c("IGHA1","IGHA2"))]
length(Cell_list1)

db<-object@meta.data
type <- paste0("Cell_list",rep(1:6))
filename <- c("EF_IGHM","EF_IGHG","EF_IGHA",
              "GC_IGHM","GC_IGHG","GC_IGHA")
fig <- list()
for(i in 1:length(type)){
  clone_db <- collapseClones(db[get(type[i]),], cloneColumn="clone_id", 
                             sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             nproc=1)
  # Create targeting model in one step using only silent mutations
  # Use consensus sequence input and germline columns
  model <- createTargetingModel(clone_db, sequenceColumn="clonal_sequence", 
                                germlineColumn="clonal_germline", vCallColumn="v_call")
  P1<-plotMutability(model, style="hedgehog")
  Pic1<-cowplot::plot_grid(plotlist = P1)
  Pic1
  fig[[i]] <- P1
  names(fig)[i] <- filename[i]
  ggsave(Pic1,filename = paste0("Hotpoint_",filename[i],"_20230927.pdf"),width = 8,height = 8)
}

###boxplot
data.total2 <- NULL
for(j in 1:length(fig)){
  P<-fig[[j]]
  name <- names(fig)[j]
  tmp<-P$A$data %>% as.data.frame() %>%
    group_by(motif) %>% 
    as.data.frame()
  tmp$type <- paste0(name,"_A")
  data.total2 <- rbind(data.total2,tmp)
  
  tmp<-P$T$data %>% as.data.frame() %>%
    group_by(motif) %>% 
    as.data.frame()
  tmp$type <- paste0(name,"_T")
  data.total2 <- rbind(data.total2,tmp)
  
  tmp<-P$C$data %>% as.data.frame() %>%
    group_by(motif) %>% 
     as.data.frame()
  tmp$type <- paste0(name,"_C")
  data.total2 <- rbind(data.total2,tmp)
  
  tmp<-P$G$data %>% as.data.frame() %>%
    group_by(motif) %>% 
    as.data.frame()
  tmp$type <- paste0(name,"_G")
  data.total2 <- rbind(data.total2,tmp)
}

data.total2$type2 <- data.total2$motif
unique(data.total2$type2)
data.total2$type2 <- gsub("WA\\/TW|WRC\\/GYW","hot",data.total2$type2)
data.total2$type2 <- gsub("SYC\\/GRS","cold",data.total2$type2)
table(data.total2$type2)
data.total2$type3 <- data.total2$type
data.total2$type3 <- gsub("_A|_G|_C|_T","",data.total2$type3)
data.total2$type4 <- paste0(data.total2$type2,"_",data.total2$type3)
data.total2$type5 <- colsplit(data.total2$type3,"_",names = c("n1","n2"))$n2
data.total2$type6 <- colsplit(data.total2$type3,"_",names = c("n1","n2"))$n1
data.total2$type5 <- factor(data.total2$type5,levels = c("IGHM","IGHG","IGHA"))
data.total2$type2 <- factor(data.total2$type2,levels = c("Neutral","cold","hot"))

library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
ggplot(data.total2,aes(x=type5,y=score,fill=type6))+geom_boxplot()+
  theme_classic()+theme(axis.text.x = element_text(angle=45,hjust=1))+
  facet_grid(.~type2)+
  stat_compare_means(method = "wilcox.test",aes(label=..p.signif..))+
  xlab("")+scale_fill_manual(values=col_flg)
ggsave("hotpoint_score_20230927.pdf",width = 10,height = 8)


###############################################################################
#'                       Manuscipt: figureS5L                                '#
###############################################################################
object<-readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")

Idents(object) <- object$groupn
object2 <- subset(object,idents = c("EF","GC"))
Idents(object2) <- object2$celltype_level3_0912
object2 <- subset(object2,idents = c("c14_PB","c15_MZB1+ASC"))
dim(object2)
object2$type_group <- paste0(object2$type,"_",object2$groupn)
object2$group <- object2$groupn
table(object2$c_call);table(object2$type);table(object2$celltype_level3_0912);table(object2$type_group)
anno <- object2@meta.data
#anno <- anno[-which(anno$c_call%in%c("","IGHD","IGHE")),]
object2 <- object2[,rownames(anno)]
Idents(object2) <- object2$type

col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
isotype_color <- c("EF"=col_flg[1],"GC"=col_flg[2])
for(i in 1:length(unique(object2$type))){
  object2.temp <- subset(object2,idents = unique(object2$type)[i])
  # db_sub=ASCGC
  # db_sub <- subset(object2, c_call %in% c("IGHG4")#  celltype_l8 %in% c("B.c12.rASC") 
  #                  )#, "IGHG","IGHD","IGHA"
  object2.temp@meta.data$c_call <- factor(object2.temp@meta.data$c_call,
                                          levels = c("IGHM","IGHG3","IGHG1","IGHG2",
                                                     "IGHG4","IGHA1","IGHA2"))
  # Collapse clonal groups into single sequence
  clones_sub <- collapseClones(object2.temp@meta.data, cloneColumn="clone_id",
                               sequenceColumn="sequence_alignment",
                               germlineColumn="germline_alignment_d_mask",
                               regionDefinition=IMGT_V, 
                               method="thresholdedFreq", minimumFrequency=0.6,
                               includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                               nproc=20)
  
  # Calculate selection scores from scratch
  baseline_sub <- calcBaseline(clones_sub, testStatistic="focused", 
                               regionDefinition=IMGT_V, 
                               #mutationDefinition=HYDROPATHY_MUTATIONS,
                               nproc=20)
  
  # Combine selection scores by time-point and isotype
  # grouped_2 <- groupBaseline(baseline_sub, groupBy=c("celltype_level3_0912", "type"))
  # testBaseline(grouped_2, groupBy="type")
  
  grouped_2 <- groupBaseline(baseline_sub, groupBy=c("group", "c_call"))
  testBaseline(grouped_2, groupBy="c_call")
  grouped_2@stats
  grouped_2@stats$c_call <- factor(grouped_2@stats$c_call,
                                levels = c("IGHM","IGHG3","IGHG1","IGHG2",
                                           "IGHG4","IGHA1","IGHA2"))
  p1=plotBaselineDensity(grouped_2,"group", groupColumn="c_call", colorElement="group", 
                         #groupColors=c("Cancer_PBMC"="#4daf4a","LN_Met"="#984ea3", "Adjacent"="#377eb8", "Cancer"="#e41a1c"),
                         #colorValues=c("Blood"="#4daf4a","LN_Met"="#984ea3", "Adjacent"="#377eb8", "Cancer"="#e41a1c"), 
                         sigmaLimits=c(-1, 1))+geom_vline(xintercept=c(0), linetype="dotted");p1
  p2=plotBaselineSummary(grouped_2, "c_call", "group", #facetBy="group"
                         groupColors=isotype_color#color.B$BCR_col
  );p2
  ggsave(paste0("c_call_PB_ASC_summaryplot_0928_",unique(object2$type)[i],".pdf"),width = 4,height = 5)
}

###############################################################################
#'                       Manuscipt: figureS5M                                '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")

Idents(object) <- object$groupn
object <- subset(object,idents=c("EF","GC"))
Idents(object) <- object$celltype_level3_0912
object <- subset(object,idents=c("c14_PB","c15_MZB1+ASC"))
dim(object)
object$group <- object$groupn
object$sample.label=paste0(object$type,"-",object$cancer,"-",object$patient,"-",object$group)


sample.size=table(object@meta.data[,"sample.label"] %>% as.vector()) %>% as.data.frame()
rownames(sample.size)<-sample.size$Var1
sample.size$type=str_split(sample.size$Var1,"-",simplify = T)[,1]
sample.size$cancer=str_split(sample.size$Var1,"-",simplify = T)[,2]
sample.size$stype=str_split(sample.size$Var1,"-",simplify = T)[,4]
#sample.size<-sample.size[sample.size$type %in% c("Adjacent","Cancer"),]
sample.size<-sample.size[sample.size$type %in% c("Adjacent","Cancer"),]
#sample.size<-sample.size[sample.size$type %in% c("Cancer_PBMC","Cancer"),]

#sample.size<-sample.size[c(GroupA,GroupB,reference),]
keep.samples<-sample.size[sample.size$Freq>=5,]#min.cells
#keep.samples$stype<-ifelse(keep.samples$cancer %in% c("COAD","STAD","GIST","BLCA"),"GC","EF")



new<-object@meta.data
new$IGHV<-NA
i=1
for(i in 1:dim(new)[1]){
  new[,"v_call"][i] %>%as.vector() %>%strsplit(split='\\*') ->tmp
  new[i,"IGHV"]<-tmp[[1]][1]
}

Usage<-prop.table(table(new$IGHV,new[,"sample.label"]));Usage[1:3,1:3]
for(j in 1:dim(Usage)[2]){
  Usage[,j]<-Usage[,j]/sum(Usage[,j])
}
Usage <- as.data.frame.array(t(Usage[,rownames(keep.samples)]))##[,rownames(keep.samples)]

############stype EF vs GC
Usage$group<-keep.samples$stype
Usage[,-dim(Usage)[2]]<-Usage[,-dim(Usage)[2]]*1000

condition <- factor(Usage$group)
colData <- data.frame(row.names=rownames(Usage), condition)
Usage<-floor(Usage[,-dim(Usage)[2]])
Usage<-t(Usage) %>% as.data.frame()
Usage=Usage+1

dds <- DESeqDataSetFromMatrix(countData = Usage,
                              colData = colData,
                              design= ~condition)
dds2 <- DESeq(dds)
res1 <- results(dds2, contrast=c("condition","EF","GC")) %>% as.data.frame()####EF在前，则EF为log2FoldChange大于0
data1<-na.omit(res1)

#######type T vs P
Usage<-prop.table(table(new$IGHV,new[,"sample.label"]));Usage[1:3,1:3]
for(j in 1:dim(Usage)[2]){
  Usage[,j]<-Usage[,j]/sum(Usage[,j])
}
Usage <- as.data.frame.array(t(Usage[,rownames(keep.samples)]))

Usage$group<-keep.samples$type
Usage[,-dim(Usage)[2]]<-Usage[,-dim(Usage)[2]]*1000

condition <- factor(Usage$group)
colData <- data.frame(row.names=rownames(Usage), condition)
Usage<-floor(Usage[,-dim(Usage)[2]])
Usage<-t(Usage) %>% as.data.frame()
Usage=Usage+1

dds <- DESeqDataSetFromMatrix(countData = Usage,
                              colData = colData,
                              design= ~condition)
dds2 <- DESeq(dds)
res2 <- results(dds2, contrast=c("condition","Cancer","Adjacent"))%>% as.data.frame()

data2<-na.omit(res2)


m<-Usage
z<-m %>% apply(2,sum)
for(i in 1:dim(m)[2]){
  m[,i]<-m[,i] / as.numeric(z[i])
}
n<-m %>% apply(1,median) %>% as.data.frame()
y<-intersect(data1 %>% rownames(),data2%>% rownames())
k<-n[y,]

data<-data.frame(X1=data1[y,]$log2FoldChange,X2=data2[y,]$log2FoldChange,row.names = y,Size=k)

data$ID<-rownames(data)
library(ggrepel)
data$color<-"0"
data[data$X1>0 &data$X2>0,]$color<-"1"
data[data$X1<0 &data$X2>0,]$color<-"2"
data[data$X1<0 &data$X2<0,]$color<-"3"
data[data$X1>0 &data$X2<0,]$color<-"4"
P1<-ggplot(data,mapping = aes(x=X1,y=X2,size=Size,color=color))+geom_point()+
  #DimPlot_theme+
  geom_text_repel(data=data[data$X1 >0 & data$X2 >0,],aes(x=X1,y=X2,label=ID))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+scale_color_manual(values = c("black","red","yellow","blue","green"));P1
return(P1)

P2<-ggplot(data,mapping = aes(x=X1,y=X2,size=Size,color=color))+geom_point()+
  #DimPlot_theme+
  theme_classic()+
  geom_text_repel(data=data[abs(data$X1) >0 & abs(data$X2) >0,],aes(x=X1,y=X2,label=ID),max.overlaps=100)+
  geom_hline(aes(yintercept=0))+xlab("EF vs GC") + ylab("Cancer vs adjacent") +
  geom_vline(aes(xintercept=0))+scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"));P2
ggsave("./review/BCR/figures_use/EF vs GC and Cancer vs Adjacent PCs IGHV usage 231003.pdf",width = 8,height = 6)


Usage <- as.data.frame(t(Usage))
Usage$type=keep.samples$type
Usage$stype=keep.samples$stype
Usage$sam=rownames(Usage)
x=melt(Usage,id.vars=c("sam","type","stype"))
x$type2 <- paste0(x$variable,"_",x$sam)
ef=c("GC"="#377eb8" ,"EF"="#e41a1c")
my_comparisons=list(c("Adjacent_EF", "Adjacent_GC"),c("Adjacent_EF" , "Cancer_EF"),
                    c("Adjacent_GC" , "Cancer_GC" ),c("Cancer_EF", "Cancer_GC"))
my_comparisons=list(c("EF","GC"))

x2<-prop.table(table(new$IGHV,new[,"sample.label"]))
x2 <- as.data.frame(x2)
x2$type2 <- paste0(x2$Var1,"_",x2$Var2)
x <- merge(x,x2,by="type2")
data1=arrange(data1,log2FoldChange)

top10 <- data1[abs(data1$log2FoldChange)>1,];top10
#x1=x[x$variable %in% c("IGHV2-5","IGHV1-3","IGHV4-59","IGHV2-70","IGHV7-4-1"),]
x1=x[x$variable %in% rownames(data),]
x1$typeN=paste0(x1$type,"_",x1$stype)

#cancer
library(ggpubr)
data2 <- x1[x1$type=="Cancer",]
p5=ggplot(data = data2,mapping = aes(x =stype ,y =Freq*100)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type)) +
  geom_boxplot(mapping = aes(fill = stype),scale = "width",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons)+#ylim(1,100)+
  scale_fill_manual(values=ef)+ labs(title = "IgHV usage")+
  facet_wrap(~variable,scales = "free_y",ncol=4)+
  xlab("type") + ylab("Freq") +#coord_cartesian(ylim = c(0, 0.25))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p5
ggsave("./review/BCR/figures_use/EF vs GC IgHV cancer_20231003.pdf",width = 20,height =40)
#Adj
data3 <- x1[x1$type=="Adjacent",]
p5=ggplot(data = data3,mapping = aes(x =stype ,y =Freq*100)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type)) +
  geom_boxplot(mapping = aes(fill = stype),scale = "width",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons)+#ylim(1,100)+
  scale_fill_manual(values=ef)+ labs(title = "IgHV usage")+
  facet_wrap(~variable,scales = "free_y",ncol=4)+
  xlab("type") + ylab("Freq") +#coord_cartesian(ylim = c(0, 0.25))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p5
ggsave("FigS5M.pdf",width = 20,height =40)  
