
################fig.S4H SHM model################
table(objectASC$c_call);table(objectASC$celltype_l3);table(objectASC$type);table(objectASC$group)

Idents(objectASC)="group";table(Idents(objectASC))
ASCEF=subset(objectASC,group=="EF")
ASCGC=subset(objectASC,group=="GC")

ASCEFT=subset(ASCEF,type %in% c("Cancer"))
ASCEFTIGG=subset(ASCEFT,c_call %in% c("IGHG1","IGHG2","IGHG3","IGHG4"))
ASCEFTIGA=subset(ASCEFT,c_call %in% c("IGHA1","IGHA2"))
ASCEFTIGM=subset(ASCEFT,c_call %in% c("IGHM"))

ASCGCT=subset(ASCGC,type %in% c("Cancer"))
ASCGCTIGG=subset(ASCGCT,c_call %in% c("IGHG1","IGHG2","IGHG3","IGHG4"))
ASCGCTIGA=subset(ASCGCT,c_call %in% c("IGHA1","IGHA2"))
ASCGCTIGM=subset(ASCGCT,c_call %in% c("IGHM"))
# Collapse sequences into clonal consensus
clone_db <- collapseClones(ASCGCT@meta.data, cloneColumn="clone_id", 
                           sequenceColumn="sequence_alignment",
                           germlineColumn="germline_alignment_d_mask",
                           nproc=10)
# Create targeting model in one step using only silent mutations
# Use consensus sequence input and germline columns
model <- createTargetingModel(clone_db, model="rs", sequenceColumn="clonal_sequence", 
                              germlineColumn="clonal_germline", vCallColumn="v_call")

# Generate hedgehog plot of mutability model
p1=plotMutability(model, nucleotides="A", style="hedgehog");p1
ggsave("./BCR_plot/ASCGCT SHM A model.pdf",width = 6,height = 6)

p2=plotMutability(model, nucleotides="C", style="hedgehog");p2
ggsave("./BCR_plot/ASCGCT SHM C model.pdf",width = 6,height = 6)

p3=plotMutability(model, nucleotides="T", style="hedgehog");p3
ggsave("./BCR_plot/ASCGCT SHM T model.pdf",width = 6,height = 6)

p4=plotMutability(model, nucleotides="G", style="hedgehog");p4
ggsave("./BCR_plot/ASCGCT SHM G model.pdf",width = 6,height = 6)

dist <- calcTargetingDistance(model)

################fig.S4I selection pressure################

# Collapse clonal groups into single sequence
clones_sub <- collapseClones(objectASCT@meta.data, cloneColumn="clone_id",
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
grouped_2 <- groupBaseline(baseline_sub, groupBy=c("celltype_l3", "type"))
testBaseline(grouped_2, groupBy="type")

grouped_2 <- groupBaseline(baseline_sub, groupBy=c("group", "c_call"))
testBaseline(grouped_2, groupBy="c_call")
grouped_2@stats
p1=plotBaselineDensity(grouped_2,"group", groupColumn="c_call", colorElement="group", 
                       sigmaLimits=c(-1, 1))+geom_vline(xintercept=c(0), linetype="dotted");p1
p2=plotBaselineSummary(grouped_2, "c_call", "group", #facetBy="group"
                       groupColors=ef#color.B$BCR_col
);p2
ggsave("./BCR_plot/objectASC c_call summaryplot.pdf",width = 4,height = 5)

################fig.S4J IgH gene usage################
object$celltype=ifelse(object$celltype_l7 %in% "B.c06.ITGB1+SwBm","Bm","AtM")
object$sample.label=paste0(object$type,"-",object$celltype,"-",object$patient)
head(str_split(object$sample.label,"-",simplify = T))

library(DESeq2)

runVDJ.diffUsage<-function(object,v.call="germline_v_call",
                           sample.label="samples",group.label="group",
                           GroupA,GroupB,reference,
                           min.cells=10){
  sample.size=table(object@meta.data[,"sample.label"] %>% as.vector()) %>% as.data.frame()
  rownames(sample.size)<-sample.size$Var1
  sample.size$type=str_split(sample.size$Var1,"-",simplify = T)[,1]
  sample.size$celltype=str_split(sample.size$Var1,"-",simplify = T)[,2]
  sample.size<-sample.size[sample.size$type %in% c("Adjacent","Cancer"),]
  
  keep.samples<-sample.size[sample.size$Freq>=5,]#min.cells
  
  
  
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
  
  ############stype AtM vs Bm
  Usage$group<-keep.samples$celltype;table( Usage$group)
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
  res1 <- results(dds2, contrast=c("condition","AtM","Bm")) %>% as.data.frame()
  data1<-na.omit(res1)
  data1=arrange(data1,log2FoldChange)
  
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
  data2=arrange(data2,log2FoldChange)
  
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
  
  P2<-ggplot(data,mapping = aes(x=X1,y=X2,size=Size,color=color))+geom_point()+
    #DimPlot_theme+
    theme_classic()+
    geom_text_repel(data=data,aes(x=X1,y=X2,label=ID))+
    geom_hline(aes(yintercept=0))+xlab("AtM vs Bm") + ylab("Cancer vs Adjacent") +
    geom_vline(aes(xintercept=0))+scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"));P2
  ggsave("./BCR_plot/AtM vs Bm and Cancer vs Adjacent IGHV usage.pdf",width = 8,height = 6)
}


