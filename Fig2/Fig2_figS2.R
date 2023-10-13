
###########Fig2E and H ####
library("sscVis")
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("data.table")
library("plyr")
library("ggpubr")


dir.startrac <- "./BCR data/"
dir.create(dir.startrac,F,T)

Idents(object)="type";table(Idents(object))

in.dat1=as.data.table(object@meta.data)

in.dat1$celltype_l3=as.character(in.dat1$celltype_l3)
in.dat1$celltype_l4=ifelse(in.dat1$celltype_l3 %in% c("B_14_MZB1_rASC","B_13_STMN1_PB"),"ASC",in.dat1$celltype_l3)

in.dat1$majorCluster <- as.character(in.dat1$celltype_l4)
in.dat1$clone.id <- in.dat1$clone_id

{
  
  #### filter out patient.cluster with number of cell < 10
  ncell.patient.cluster <- sort(unclass(in.dat1[,table(sprintf("%s.%s",patient,majorCluster))]))
  in.dat1 <- in.dat1[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=1,]
  
  dim(in.dat1)

  in.dat1$Cell_Name=in.dat1$BCR_id
  in.dat1$clone.id=in.dat1$clone_id
  in.dat1$patient=in.dat1$patient
  in.dat1$majorCluster=in.dat1$celltype_l4
  in.dat1$loc=in.dat1$type
  
  
  
  if(!file.exists(sprintf("%s/B.out.startrac_Tumor.nperm10020221228.rds",dir.startrac))){
    tic("Startrac.run")
    out <- Startrac.run(in.dat1, proj="panC",verbose=F,cores=15,n.perm=100)##n.perm=1000 原始是1000，一跑内存就超
    toc()
    dir.startrac <- "./BCR data/"
    saveRDS(out,sprintf("%s/object_total.out.startrac.nperm100_20230501_B13 and B14 merge ASC.rds",dir.startrac))
  }
  out <- readRDS(sprintf("%s/B.out.startrac.nperm1000.rds",dir.startrac))
  
  if(!file.exists(sprintf("%s/B.out.startrac.nperm1000.rds",dir.startrac))){
    tic("Startrac.run")
    out <- Startrac.run(in.dat1[stype=="CD8",], proj="panC",verbose=F,cores=6,n.perm=1000)
    toc()
    saveRDS(out,sprintf("%s/B.out.startrac.nperm1000.rds",dir.startrac))
  }
  out <- readRDS(sprintf("%s/B.out.startrac.nperm1000.rds",dir.startrac))
  
}


############Single cancer##############
cancerType.vec <- in.dat1[,unique(cancer)];cancerType.vec;length(cancerType.vec)

if(!file.exists(sprintf("%s/B.out.startrac.total.rds",dir.startrac)))
{
  res.byCancerType <- llply(cancerType.vec,function(x){
    Startrac.run(in.dat1[cancer==x,],
                 proj=x,verbose=F,cores=15,n.perm=NULL)
  })
  names(res.byCancerType) <- cancerType.vec
  saveRDS(res.byCancerType,sprintf("%s/B.out.startrac total.rds",dir.startrac))
}

# load data
{
    in.startrac.file.list <- list(
                                "MZB1less1rm" ="./BCR data/B.out.startrac total rm MZB1less1 230407.rds"
  )
}

#######  startrac pTrans across cancerTypes ########
## version 1 (used in paper) (Fig. 2D, Fig. 3C ??)
{
  z.max <- 2
  
  dat.startrac <- readRDS(in.startrac.file.list[[ "MZB1less1rm"]])
    
    
    ####### pIndex (tran)
      
      #mcls.x <- "CD8.c12.Tex.CXCL13"
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
      dat.plot.mtx.tmp[ abs(dat.plot.mtx) < 0.01] <- 0
      f.cancerType <- rowSds(dat.plot.mtx.tmp)==0
      dat.plot.mtx <- dat.plot.mtx[!f.cancerType,,drop=F]
      dat.plot.mtx.tmp <- dat.plot.mtx.tmp[!f.cancerType,,drop=F]
      
      pTrans.hclust.row <- run.cutree(dat.plot.mtx.tmp[,f.mcls],
                                      k=3,method.distance="cosine",method.hclust="ward.D2")
      coul <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(30))
      
      pdf("./BCR_plot/B_c09_DUSP4_AtM pTrans onlytumor0326.pdf",width = 8,height = 8)
      pheatmap::pheatmap(dat.plot.mtx.tmp,,z.lo=-z.max,z.hi=z.max,z.len=50,color = coul,clustering_method = "ward.D2")#centroid
      dev.off()
      
     

############pTrans####
      out <- readRDS("./BCR data/objectB.out.startrac.nperm100.rds")
      
      pIndex.sig.tran <- as.data.table(out@pIndex.sig.tran)
      
      dat.plot <- pIndex.sig.tran
      dat.plot[,index:=as.character(index)]
      color.B$cluster.name <- color.B$celltype
      names(color.B$cluster.name) <- names(color.B$celltype)
      f.col <- !is.na(names(color.B$cluster.name))
      color.B$cluster.name <- color.B$cluster.name[f.col]
      
      mcls.moi <- c("B_14_MZB1_rASC")
      
      l_ply(mcls.moi,function(x){
        
        dat.tmp <- dat.plot[majorCluster=="ASC" & aid=="panC" &
                              !is.na(value),]
        dat.med.tmp <- dat.tmp[order(value),]
        dat.tmp[,index:=factor(index,levels=dat.med.tmp$index)]
        dat.tmp[,index.name:=factor(index,levels=dat.med.tmp$index.name)]
        p <- ggplot(dat.tmp, aes(index,value)) +
          geom_col(fill="steelblue",col=NA,width=0.8) +
          geom_hline(yintercept=0.03,linetype="dashed",color="black",alpha=0.8) +
          xlab("") + ylab(sprintf("pTrans Of %s","ASC")) +
          theme_pubr() +
          theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1));p
        ggsave(sprintf("%s pTrans.Fig.barplot nperm100 B13 B14 merge.%s.pdf",out.prefix,"PCs"),width=5.5,height=4.5)
        if(mcls.moi=="B.c12.rASC"){
          ggsave(sprintf("%s.pTrans.Fig.barplot.byPatientF.%s.Fig2C.v2.pdf",out.prefix,x),width=3.8,height=3.8)
        }
      },.parallel=T)


###########Fig2K.BCR Abundance and Diversity analysis ###

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
curve <- estimateAbundance(top150, group="type", ci=0.95, nboot=1000, clone="clone_id")
# Plots a rank abundance curve of the relative clonal abundances
type=c( "Cancer_PBMC"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")

plot(curve, colors = class, legend_title="type")

####Diversity analysis###########
clones <- countClones(objectnPC@meta.data, group="type")
head(clones, 5)

# Partitions the data on the sample column
# Calculates a 95% confidence interval via 200 bootstrap realizations
curve <- estimateAbundance(top150, group="groupn", ci=0.95, nboot=100, clone="clone_id")
# Plots a rank abundance curve of the relative clonal abundances
type=c( "Cancer_PBMC"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")

plot(curve, colors = class, legend_title="Type diversity")
ggsave("./BCR_plot/Top150 EF vs GC objectPC clonal Abundance.pdf",width = 6,height = 6)


sample_curve <- alphaDiversity(top150, group="groupn", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)

sample_main <- paste0("Type diversity")
type=c( "Cancer_PBMC"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
plot(sample_curve, colors=class, main_title=sample_main, 
     legend_title="Type diversity")

ggsave("./BCR_plot/Top150 EF vs GC objectPC diversity.pdf",width = 6,height = 6)



###########Fig2K. CSR###

library(reshape2)
library(ggplot2)
library(ggraph)
library(tidygraph)
EF <- readRDS("./objectASCEF_meta.data230218.rds")
unique(EF$type)
set.seed(123)
for(j in 1:length(unique(EF$type))){
  EF_Cancer <- EF[EF$type==unique(EF$type)[j],]
  # EF_Cancer <- EF_Cancer[sample(rownames(EF_Cancer),3000),]
  data2 <- as.data.frame(table(EF_Cancer$clone_id,EF_Cancer$c_call))
  data <- dcast(data2,Var1~Var2)
  rownames(data) <- data$Var1
  data <- data[,-1]
  num <- colSums(data)
  num_percent <- as.data.frame(num/sum(num))
  colnames(num_percent) <- c('clone')
  num_percent$name <- rownames(num_percent)
  num_percent$type <- 'EF'
  data <- data[rowSums(data>0)>1,]
  
  compare <- as.data.frame(combn(colnames(data),2))
  sum <- matrix(nrow = 28,ncol = 3)
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
  edge[28,3]="5"
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
    scale_edge_width(range = c(0,3))+
    ggtitle(paste0("AtM_EF_",unique(EF$type)[j]))+
    scale_size_continuous(range = c(10,20))
  ggsave(paste0("ASC_EF_",unique(EF$type)[1],".pdf"),width=5,height=4)
}

gridplot_grid(p1, p2, p3, p4)
######

objectASCGC=subset(objectASC,idents="GC");table(Idents(objectASCGC))
GC=objectASCGC@meta.data

GC <- readRDS("./objectASCGC_meta.data230218.rds")
unique(GC$type)

for(j in 1:length(unique(GC$type))){
  GC_Cancer <- GC[GC$type==unique(GC$type)[j],]
  data2 <- as.data.frame(table(GC_Cancer$clone_id,GC_Cancer$c_call))
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
  sum <- matrix(nrow = 28,ncol = 3)
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
  edge2[28,3]="5"
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
    scale_edge_width(range = c(0,3))+
    ggtitle(paste0("AtM_GC_",unique(GC$type)[j]))+
    scale_size_continuous(range = c(10,20))
  ggsave(paste0("ASC_GC_",unique(GC$type)[j],".pdf"),width=5,height=4)
}


###########fig.S2G. Clonetype frequency###
new<-objectT@meta.data
colnames(new)
select(new,c("cell_id","c_call","clone_id")) %>% head()
new$clone_id<-as.vector(new$clone_id)
q<-as.data.frame(table(new$clone_id))
new$CloneType_Freq<-NA

q=arrange(q,desc(Freq))

library(future)
plan("multiprocess", workers =12)
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

f4=dittoBarPlot(object, "CloneType_Freq_levels",group.by="celltype_l3",retain.factor.levels = T,main = "CloneType_Freq_levels",color.panel=
                  c("#ffffe5","#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"))+ylim(0,1);f4

ggsave("./BCR_plot/objectB Tumor CloneType_Freq_levels_Type freq.pdf",width = 6,height = 6)

