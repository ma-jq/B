if(T){
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(dittoSeq)
  library(reshape2)
  library(dplyr)
  library(stringr)
  library(ggpubr)
  library(RColorBrewer)
  library(alakazam)
  library(ggprism)
}


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

###############################################################################
#'                   Manuscipt: figure4D&E;S8J&K                             '#
###############################################################################

data <- readRDS("../scRNA_data/panB_scRNA_processed_data.rds")
Idents(data) <- data$celltype
data <- subset(data,idents=c("B.14.Plasmablast","B.15.Plasma cell"),invert=TRUE)
##
set.seed(123)
cellname <- sample(colnames(data),50000)
scdata_sample <- data[,cellname]
DimPlot(scdata_sample)
scdata_sample <- FindVariableFeatures(scdata_sample)
#remove noncoding
genes<-data.frame(data.table::fread("../Additional_data/biomart_mart_export.txt",header=T,sep="\t"))

table(genes$GeneType)
genes_sub<-subset(genes,GeneType!="lncRNA") #processed_pseudogene lncRNA
genes_sub<-subset(genes_sub,GeneType!="processed_pseudogene") #processed_pseudogene lncRNA
genes_sub_all<-c(as.character(genes_sub$gene),as.character(genes_sub$GeneSynonym))
sum(scdata_sample@assays$RNA@var.features %in% genes_sub_all); length(scdata_sample@assays$RNA@var.features)
scdata_sample@assays$RNA@var.features = scdata_sample@assays$RNA@var.features[(scdata_sample@assays$RNA@var.features %in% genes_sub_all)]
dim(scdata_sample)

####WYC blacklist
load("../Additional_data/genes_black_WYC.rda")
sum(scdata_sample@assays$RNA@var.features %in% genes_black); length(scdata_sample@assays$RNA@var.features)
scdata_sample@assays$RNA@var.features = scdata_sample@assays$RNA@var.features[!(scdata_sample)]
mito.genes <- rownames(scdata_sample@assays$RNA)[grep("^MT-",rownames(scdata_sample@assays$RNA))]
scdata_sample@assays$RNA@var.features =scdata_sample@assays$RNA@var.features[!(scdata_sample@assays$RNA@var.features %in% mito.genes)]

dim(scdata_sample);length(scdata_sample@assays$RNA@var.features)

scdata_sample<- NormalizeData(scdata_sample)
scdata_sample <- ScaleData(scdata_sample)
scdata_sample <- RunPCA(scdata_sample, verbose = FALSE)

scdata_sample <- harmony::RunHarmony(scdata_sample,"patient", plot_convergence = TRUE)

ElbowPlot(objTN2, ndims = 50,reduction = "harmony")

scdata_sample <- Seurat::RunUMAP(scdata_sample,reduction = "harmony", dims = 1:15) 

DimPlot(scdata_sample,group.by = "celltype",label=T,cols = my36colors)
dim(scdata_sample)
anno_sample <- scdata_sample@meta.data
###monocle3
data <- GetAssayData(scdata_sample,assay = "RNA",slot = "counts")
data <- data[rowSums(data>0)>=3,]
dim(data)
cell_metadata <- scdata_sample@meta.data
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
int.embed <- Embeddings(scdata_sample,reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "type")
p2

cds <- cluster_cells(cds,resolution = 0.001,k=40,random_seed=18,verbose=T)
cds@clusters$UMAP$clusters <- scdata_sample$celltype
plot_cells(cds,color_cells_by = "partition")

p1 <- plot_cells(cds,group_cells_by = 'cluster')
p1

cds <- learn_graph(cds, verbose =T,
                   use_partition=T,close_loop=F,learn_graph_control=
                     list(minimal_branch_len=10,rann.k=40)  )
p <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=T,
                label_leaves=T, label_branch_points=T,cell_size = 0.5,group_label_size=4)+
  scale_color_manual(values=my36colors)
p
cds <- order_cells(cds)

unique(cds@colData$celltype)

p1 <- plot_cells(cds,color_cells_by = "pseudotime",label_branch_points = FALSE,
                 label_leaves = F)
p1
p2 <- DimPlot(scdata_sample,group.by = "celltype",cols = my36colors,label = T)
p2
p2|p1

p3 <- plot_cells(cds = cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           trajectory_graph_color = "white",
           trajectory_graph_segment_size = 0.5,
           graph_label_size = 2,
           cell_size = 1,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F) 
p2|p3
ggsave("monocle3.pdf",width = 18,height = 8)


#####boxplot#####
df2 <- data.frame(pseudotime=pseudotime(cds),
                  celltype=cds@colData$celltype)
unique(df2$celltype)

df2 <- df2[df2$celltype%in%c("B.08.ITGB1+SwBm" ,"B.09.DUSP4+AtM" , "B.07.CCR7+ACB3", 
                             "B.06.NR4A2+ACB2" ,"B.10.ENO1+Pre_GCB","B.11.SUGCT+DZ_GCB","B.12.LMO2+LZ_GCB",
                             "B.01.TCL1A+naiveB", "B.05.EGR1+ACB","B.13.Cycling_GCB" ),]
ggplot(df2,aes(y=celltype,x=pseudotime,color=celltype))+geom_boxplot()+
  theme_bw()+scale_color_manual(values = my36colors[c(1,5:13)])
ggplot(df2,aes(y=celltype,x=pseudotime,color=celltype))+geom_boxplot()+
  theme_bw()+scale_color_manual(values = my36colors[c(1:13)])
ggsave("boxplot_pseudotime.pdf",width = 9,height = 6)

#######trajectory plot######
Track_genes_sig <- c("PDCD1","ENTPD1","HAVCR2",
                     "CD24","CD38","SELL")
gene2 <- c("XRCC6","XRCC5","APEX1","POLD2","AICDA","APEX2")
Track_genes_sig <- c(Track_genes_sig,gene2)
Track_genes_sig <- c("GLS","GLA","PTGES3")
####

colData(cds)$pseudotime <- pseudotime(cds)
colnames(scdata_sample)[1:5]
colData(cds)[1:5]
data <- scdata_sample[Track_genes_sig,]
data <- GetAssayData(data@assays$RNA,slot='data')
data <- as.data.frame(data)

df <- t(data)
df <- as.data.frame(df)
df$pseudotime <- cds@colData$pseudotime 
df$celltype <- cds@colData$celltype
unique(df$celltype)

##
library(dplyr)
md <- as_tibble(cds@colData, rownames = NA) %>% select(celltype, pseudotime)
path1 <- c("B.01.TCL1A+naiveB","B.05.EGR1+ACB","B.06.NR4A2+ACB2","B.07.CCR7+ACB3","B.09.DUSP4+AtM")
path2 <- c("B.01.TCL1A+naiveB","B.05.EGR1+ACB","B.06.NR4A2+ACB2","B.07.CCR7+ACB3","B.09.DUSP4+AtM",
           "B.11.SUGCT+DZ_GCB","B.12.LMO2+LZ_GCB","B.13.Cycling_GCB")
#path2 <- c("c01_TCL1A+naiveB","c05_EGR+ACB","c06_NR4A2+ACB","c07_CCR7+ACB","c08_ITGB1+SwBm")
pathL <- list(path1 = path1, path2 = path2)
library(dplyr)
library(tidyr)
totalMD <-NULL
for(li in 1:length(pathL)){
  tempName <- names(pathL[li])
  tempPath <- pathL[[li]]
  tempdf <- df[df$celltype%in%tempPath,]
  tempdf$path <- names(pathL)[li]
  totalMD <- rbind(totalMD,tempdf)
}

###remove outlier
detect_outlier <- function(x) {
  Quantile1 <- quantile(x, probs=.25)
  Quantile3 <- quantile(x, probs=.75)
  IQR = Quantile3-Quantile1
  #x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
  x > Quantile3  | x < Quantile1 
}

# create remove outlier function
remove_outlier <- function(dataframe,
                            columns=names(dataframe)) {
  for (col in columns) {
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  return(dataframe)
  print("Remove outliers")
}
colnames(totalMD)
unique(totalMD$celltype)

# totalMD.clean <- NULL
# for(i in 1:length(unique(totalMD$celltype))){
#   temp <- totalMD[totalMD$celltype==unique(totalMD$celltype)[i],]
#   temp.clean <- remove_outlier(temp,c("pseudotime"))
#   totalMD.clean <- rbind(totalMD.clean,temp.clean)
# }

totalMD.clean2 <- NULL
for(i in 1:length(unique(totalMD.clean$path))){
  temp.clean <- totalMD.clean[totalMD.clean$path==unique(totalMD.clean$path)[i],]
  temp.clean$pseudotime <- (temp.clean$pseudotime-min(temp.clean$pseudotime))/(max(temp.clean$pseudotime)-min(temp.clean$pseudotime))
  totalMD.clean2 <- rbind(totalMD.clean2,temp.clean)
}

library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
colnames(totalMD.clean2)
totalMD.clean.use <- melt(totalMD.clean2,id.vars = c("celltype", "path",
                                                    "pseudotime")) 
fig <- list()
for(i in 1:length(unique(totalMD.clean.use$variable))){
  temp <- totalMD.clean.use[totalMD.clean.use$variable==unique(totalMD.clean.use$variable)[i],]
  p <- ggplot(temp) +
    stat_smooth(aes(x = pseudotime, y = value,color=path),method = "lm", formula = y ~ poly(x, 3),se=TRUE)+
    theme_classic()+scale_color_manual(values=col_flg)+ggtitle(unique(totalMD.clean.use$variable)[i])
  fig[[i]] <- p
  }
fig[['nrow']] <- 3
fig[['ncol']] <- 4
library(gridExtra)
pdf('features_trajectory.pdf', width = 20, height = 10)
do.call('grid.arrange', fig)
dev.off()

high_affinity=c("BATF","GARS","GART","LER3","MIF","MYC","SPP1","UCK2",
                "CD320","TIMD2","TNFRSF8",
                "AURKA","BUB1","CCNA2","CCNB1","CCNB2","CCND2","CDC20","CDC25C","KIF22","PLK1",
                "NFIL3","PML")
Low_affinity=c("PPP1R15A","RGS1","CCR6","CD22","CD38","CD72","FCER2A","ICOSL","PACAM1","SIGLECG","TLRL","TNFRSF18",
               "BACH2","EGR3","ELK4","FOXP1","JUN","NR4A1","REL")
exhaustion_genes = c('PDCD1','CD160','FASLG','CD244','LAG3','TNFRSF1B','CCR5','CCL3',
                     'CCL4','CXCL10','CTLA4','LGALS1','LGALS3','PTPN13','RGS16','ISG20',
                     'MX1','IRF4','EOMES','PBX3','NFATC1','NR4A2','CKS2','GAS2',
                     'ENTPD1','CA2',"CD52", "APOE", "PTLP", "PTGDS", "PIM2", "DERL3")


Bactivated_genes=c("CD69","CD83","IER2","DUSP2","IL6","NR4A2","JUN","CCR7","GPR183")
BCSR_genes=c("APEX1","APEX2","XRCC5","XRCC6","POLD2","AICDA")
CSR_m=c("APEX1","XRCC5","XRCC6","POLD2","POLE3", #CSR machinery
        "NCL","NME2","DDX21",#IgH locus
        "NPM1","SERBP1",#CSR interactors
        "MIR155HG","HSP90AB1",#AICDA/AICDA stability
        "BATF","HIVEP3","BHLHE40","IRF4")#TF

feature <- list(high_affinity=high_affinity,
                Low_affinity=Low_affinity,
                exhaustion_genes=exhaustion_genes,
                Bactivated_genes=Bactivated_genes,
                BCSR_genes=BCSR_genes,
                CSR_m=CSR_m)
# names(feature) <- c("high_affinity","Low_affinity","exhaustion_genes",
#                     "Bactivated_genes","BCSR_genes","CSR_m")
scdata_sample <- AddModuleScore(scdata_sample,features = feature,name=c("high_affinity","Low_affinity","exhaustion_genes",
                                                                        "Bactivated_genes","BCSR_genes","CSR_m"))
ggplotdata <- data.frame(UMAP1=scdata_sample@reductions$umap@cell.embeddings[,1],
                         UMAP2=scdata_sample@reductions$umap@cell.embeddings[,2])
scdata_sample$pseudotime <- cds@colData$pseudotime

ggplotdata2 <- scdata_sample@meta.data[,c(40:45)]
ggplotdata <- cbind(ggplotdata,ggplotdata2)
colnames(ggplotdata)
max(ggplotdata$high_affinity1)
#ggplotdata$high_affinity1[ggplotdata$high_affinity1<0.4] <- 0
p1 <- ggplot() +
  geom_point(data = ggplotdata, 
             aes(x = UMAP1, y = UMAP2, color = high_affinity1), size = 0.01)+
  viridis::scale_color_viridis(option = "viridis")+
  #scale_color_gradient(low="grey",high="red")+
  #scale_color_gradientn(colors = c( "#0000FF", "#8888FF", "#AAAAFF","#FFFFFF", "#FF8888", "#FF5555", "#FF0000")) +
  theme_void()

p2 <- ggplot() +
  geom_point(data = ggplotdata, 
             aes(x = UMAP1, y = UMAP2, color = Low_affinity2), size = 0.01)+
  viridis::scale_color_viridis(option = "viridis")+
  #scale_color_gradient(low="grey",high="red")+
  #scale_color_gradientn(colors = c("#0000FF", "#8888FF", "#AAAAFF", "#FFFFFF", "#FF8888", "#FF5555", "#FF0000")) +
  theme_void()
p2
p3 <- ggplot() +
  geom_point(data = ggplotdata, 
             aes(x = UMAP1, y = UMAP2, color = exhaustion_genes3), size = 0.01)+
  viridis::scale_color_viridis(option = "viridis")+
  #scale_color_gradient(low="grey",high="red")+
  #scale_color_gradientn(colors = c( "#0000FF", "#8888FF", "#AAAAFF", "#FFFFFF", "#FF8888", "#FF5555", "#FF0000")) +
  theme_void()
p3

p4 <- ggplot() +
  geom_point(data = ggplotdata, 
             aes(x = UMAP1, y = UMAP2, color = Bactivated_genes4), size = 0.01)+
  viridis::scale_color_viridis(option = "viridis")+
  #scale_color_gradient(low="grey",high="red")+
  #scale_color_gradientn(colors = c( "#0000FF", "#8888FF", "#AAAAFF", "#FFFFFF", "#FF8888", "#FF5555", "#FF0000")) +
  theme_void()
p4

p5 <- ggplot() +
  geom_point(data = ggplotdata, 
             aes(x = UMAP1, y = UMAP2, color = BCSR_genes5), size = 0.1)+
  viridis::scale_color_viridis(option = "viridis")+
  #scale_color_gradient(low="grey",high="red")+
  #scale_color_gradientn(colors = c( "#0000FF", "#8888FF", "#AAAAFF", "#FFFFFF", "#FF8888", "#FF5555", "#FF0000")) +
  theme_void()
p5

ggplotdata$CSR_m6 <- (ggplotdata$CSR_m6-min(ggplotdata$CSR_m6))/(max(ggplotdata$CSR_m6)-min(ggplotdata$CSR_m6))
max(ggplotdata$CSR_m6);min(ggplotdata$CSR_m6)
p6 <- ggplot() +
  geom_point(data = ggplotdata, 
             aes(x = UMAP1, y = UMAP2, color = CSR_m6), size = 0.1)+
  viridis::scale_color_viridis(option = "viridis")+
  #scale_color_gradientn(colors = c( "#0000FF", "#8888FF", "#AAAAFF", "#FFFFFF", "#FF8888", "#FF5555", "#FF0000"))+
  #values = c(0, -a/3/(b-a), -a*2/3/(b-a), -a/(b-a), (1-a/(b-a))/3, (1-a/(b-a))*2/3, 1)) +
  theme_void()
p6

p1|p2|p3|p4|p5|p6
ggsave("./module_featureplot.pdf",width = 40,height = 6)

###############################################################################
#'                       Manuscipt: figureS8B                                '#
###############################################################################
object <- readRDS("../BCR_data/panB_BCR_processed_data.rds")
Idents(object) <- object$celltype
unique(object$celltype)
object <- subset(object,idents=c("B.01.TCL1A+naiveB","B.08.ITGB1+SwBm",
                                 "B.09.DUSP4+AtM","B.12.LMO2+LZ_GCB","B.11.SUGCT+DZ_GCB",
                                 "B.13.Cycling_GCB"))
object$celltype <- as.character(object$celltype)
my36colors <-c("#3d7bb0" ,"#a2c8db" ,"#5aa554", "#afd195", "#d83f36", "#e99997" ,"#f2bc7c", "#e68740")
names(my36colors) <- c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4")

object$c_call <- factor(object$c_call,levels=c("IGHA1","IGHA2",
                                               "IGHG1","IGHG2","IGHG3","IGHG4",
                                               "IGHD","IGHM"
))

fig <- list()
for(i in 1:length(unique(object$celltype))){
  Idents(object) <- object$celltype
  temp <- subset(object,idents=unique(object$celltype)[i])
  table(temp$c_call)
  temp$type <- factor(temp$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
  temp$celltype <- factor(temp$celltype
                                        ,levels=c("B.01.TCL1A+naiveB","B.08.ITGB1+SwBm",
                                                  "B.09.DUSP4+AtM","B.12.LMO2+LZ_GCB","B.11.SUGCT+DZ_GCB",
                                                  "B.13.Cycling_GCB"))
  #temp <- subset(temp,idents=c("IGHD","IGHE",""),invert=TRUE)
  f5=dittoBarPlot(temp, "c_call",group.by="type",
                  retain.factor.levels = T,main = "type",color.panel= my36colors)+
    ylim(0,1)+ggtitle(unique(object$celltype)[i]);f5
  fig[[i]] <- f5
  
}
fig[['nrow']] <- 1
fig[['ncol']] <- 6
pdf('FigS8B.pdf', width = 20, height = 8)
do.call('grid.arrange', fig)
dev.off()


###############################################################################
#'                       Manuscipt: figureS8C                                '#
###############################################################################
object <- readRDS("../BCR_data/panB_BCR_processed_data.rds")
Idents(object)="celltype";table(Idents(object))
B89=subset(object,idents=c("B.09.DUSP4+AtM","B.08.ITGB1+SwBm"))
table(B89$c_call)
B89$type_celltype <- paste0(B89$type,"-",B89$celltype)
Idents(B89) <- B89$c_call
unique(B89$c_call)
#B89 <- subset(B89,idents=c("IGHA1","IGHG3","IGHA2" ,"IGHG1", "IGHG2", "IGHM","IGHG4"))
B89$c_call <- factor(B89$c_call,levels = c("IGHA1","IGHA2",
                                           "IGHG1","IGHG2",
                                           "IGHG3","IGHM",
                                           "IGHD","IGHG4"
                                           ))
table(B89$c_call)
data <- as.data.frame(table(B89$type,as.character(B89$celltype),B89$c_call))
ggplot(data,aes(x=Var2,fill=Var3,y=Freq))+geom_bar(stat='identity',position = 'fill')+
  theme_classic()+facet_grid(.~Var1)
patient=B89$patient %>% unique();head(patient)
B89$type <- factor(B89$type,levels = c("Blood","Adjacent","Cancer","LN_Met"))
B89proN=list()
#i=1
for(i in 1:61){
  a=B89[,B89$patient %in% patient[i]];table(B89$patient)
  B89pro <- prop.table(table(a$c_call,a$type_celltype), margin = 2) %>% as.data.frame()
  B89pro=B89pro[!(is.na(B89pro$Freq)) ,]
  B89pro$patient=patient[i]
  B89proN[[i]]=B89pro
}
B89_pro=do.call(rbind,B89proN);head(B89_pro)
B89_pro$type=str_split(B89_pro$Var2,"-",simplify = TRUE)[,1]
B89_pro$celltype=str_split(B89_pro$Var2,"-",simplify = TRUE)[,2]
unique(B89_pro$Var2)
B89_pro$Var2=factor(B89_pro$Var2,levels = c("Blood-B.08.ITGB1+SwBm","Blood-B.09.DUSP4+AtM",
                                            "Adjacent-B.08.ITGB1+SwBm" ,"Adjacent-B.09.DUSP4+AtM",
                                            "LN_Met-B.08.ITGB1+SwBm","LN_Met-B.09.DUSP4+AtM" ,
                                            "Cancer-B.08.ITGB1+SwBm","Cancer-B.09.DUSP4+AtM"
))
#my_comparisons=list(c("c08_ITGB1+SwBm","c09_AtM"))
my_comparisons=list(c("Blood-B.09.DUSP4+AtM","Blood-B.08.ITGB1+SwBm"),
                      c("Adjacent-B.09.DUSP4+AtM","Adjacent-B.08.ITGB1+SwBm") ,
                      c("Cancer-B.09.DUSP4+AtM" , "Cancer-B.08.ITGB1+SwBm"),
                      c("LN_Met-B.09.DUSP4+AtM" , "LN_Met-B.08.ITGB1+SwBm")
)
#col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(8)
col_flg <- c("#4DAF4A","#377EB8","#984EA3","#E41A1C")
B89_pro$type <- factor(B89_pro$type,levels =  c("Blood","Adjacent","LN_Met","Cancer"))
p1=ggplot(data = B89_pro,mapping = aes(x =Var2,y =Freq)) +
  geom_boxplot(mapping = aes(fill =type))+
  stat_compare_means(comparisons = my_comparisons,aes(label=..p.signif..))+
  #scale_fill_npg()+
  scale_fill_manual(values=col_flg)+#
  #labs(title = "IGHV total_Mut_freq")+
  facet_wrap(~Var1,scales = "free_y",ncol=4)+
  xlab("Isotype") + ylab("Mutation frequency") +#coord_cartesian(ylim = c(0, 0.35))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p1
ggsave("FigS8C.pdf",width = 12,height = 8)

###############################################################################
#'                       Manuscipt: figureS8F                                '#
###############################################################################
object <- readRDS("../BCR_data/panB_BCR_processed_data.rds")
Idents(object) <- object$celltype
unique(object$celltype)
object <- subset(object,idents=c("B.01.TCL1A+naiveB","B.08.ITGB1+SwBm",
                                 "B.09.DUSP4+AtM","B.12.LMO2+LZ_GCB","B.11.SUGCT+DZ_GCB",
                                 "B.13.Cycling_GCB"))
object$celltype <- as.character(object$celltype)
my_comparisons=list(c("B.08.ITGB1+SwBm_Blood","B.09.DUSP4+AtM_Blood"),
                    c("B.08.ITGB1+SwBm_Adjacent","B.09.DUSP4+AtM_Adjacent"),
                    c("B.08.ITGB1+SwBm_Cancer","B.09.DUSP4+AtM_Cancer"),
                    c("B.08.ITGB1+SwBm_LN_Met","B.09.DUSP4+AtM_LN_Met"),
                    c("B.09.DUSP4+AtM_Cancer","B.11.SUGCT+DZ_GCB_Cancer"),
                    c("B.09.DUSP4+AtM_Cancer","B.12.LMO2+LZ_GCB_Cancer"),
                   c("B.09.DUSP4+AtM_Cancer","B.13.Cycling_GCB_Cancer"),
                   c("B.09.DUSP4+AtM_Adjacent","B.11.SUGCT+DZ_GCB_Adjacent"),
                   c("B.09.DUSP4+AtM_Adjacent","B.12.LMO2+LZ_GCB_Adjacent"),
                   c("B.09.DUSP4+AtM_Adjacent","B.13.Cycling_GCB_Adjacent"))
object$type_group <- paste0(object$celltype,"_",object$type)
object$type <- factor(object$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
unique(object$type_group)
type <- paste0("B.01.TCL1A+naiveB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type2 <-  paste0("B.08.ITGB1+SwBm","_",c("Blood","Adjacent","LN_Met","Cancer"))
type3 <-  paste0("B.09.DUSP4+AtM","_",c("Blood","Adjacent","LN_Met","Cancer"))
type4 <-  paste0("B.11.SUGCT+DZ_GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type5 <-  paste0("B.12.LMO2+LZ_GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type6 <-  paste0("B.13.Cycling_GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type <- c(type,type2,type3,type4,type5,type6)
object$type_group <- factor(object$type_group,
                            levels = type)
p1=ggplot(data = object@meta.data,mapping = aes(x =type_group,y =mu_freq)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type)) +
  geom_boxplot(mapping = aes(fill = type),scale = "width")+
  stat_compare_means(comparisons = my_comparisons,aes(label=..p.signif..))+
  #scale_fill_npg()+
  scale_fill_manual(values=col_flg)+
  labs(title = "IGHV total_Mut_freq")+#facet_wrap(~c_call,scales = "free_y",ncol=4)+
  xlab("Isotype") + ylab("Mutation frequency") +#coord_cartesian(ylim = c(0, 0.2))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p1
ggsave("FigS8F.pdf",width = 8,height = 8)

###############################################################################
#'                      Manuscipt: figureS8G&H                               '#
###############################################################################

object <- readRDS("../BCR_data/panB_BCR_processed_data.rds")
Idents(object) <- object$celltype
object <- subset(object,idents=c("B.01.TCL1A+naiveB","B.08.ITGB1+SwBm",
                                 "B.09.DUSP4+AtM","B.12.LMO2+LZ_GCB","B.11.SUGCT+DZ_GCB",
                                 "B.13.Cycling_GCB"))
cdr3aa=aminoAcidProperties(object@meta.data, seq="cdr3", trim=TRUE,
                           label="cdr3")
# Define a ggplot theme for all plots
tmp_theme <- theme_bw() + theme(legend.position="bottom")
cdr3aa$type_group <- paste0(cdr3aa$celltype,"_",cdr3aa$type)

my_comparisons=list(c("B.08.ITGB1+SwBm_Blood","B.09.DUSP4+AtM_Blood"),
                    c("B.08.ITGB1+SwBm_Adjacent","B.09.DUSP4+AtM_Adjacent"),
                    c("B.08.ITGB1+SwBm_Cancer","B.09.DUSP4+AtM_Cancer"),
                    c("B.08.ITGB1+SwBm_LN_Met","B.09.DUSP4+AtM_LN_Met"),
                    c("B.09.DUSP4+AtM_Cancer","B.11.SUGCT+DZ_GCB_Cancer"),
                    c("B.09.DUSP4+AtM_Cancer","B.12.LMO2+LZ_GCB_Cancer"),
                    c("B.09.DUSP4+AtM_Cancer","B.13.Cycling_GCB_Cancer"),
                    c("B.09.DUSP4+AtM_Adjacent","B.11.SUGCT+DZ_GCB_Adjacent"),
                    c("B.09.DUSP4+AtM_Adjacent","B.12.LMO2+LZ_GCB_Adjacent"),
                    c("B.09.DUSP4+AtM_Adjacent","B.13.Cycling_GCB_Adjacent"))
type <- paste0("B.01.TCL1A+naiveB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type2 <-  paste0("B.08.ITGB1+SwBm","_",c("Blood","Adjacent","LN_Met","Cancer"))
type3 <-  paste0("B.09.DUSP4+AtM","_",c("Blood","Adjacent","LN_Met","Cancer"))
type4 <-  paste0("B.11.SUGCT+DZ_GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type5 <-  paste0("B.12.LMO2+LZ_GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type6 <-  paste0("B.13.Cycling_GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type <- c(type,type2,type3,type4,type5,type6)
cdr3aa$type_group <- factor(cdr3aa$type_group,
                            levels = type)
cdr3aa$type <- factor(cdr3aa$type,levels =  c("Blood","Adjacent","LN_Met","Cancer"))
# Generate plots for all four of the properties
col_flg <- c("#4DAF4A","#377EB8","#984EA3","#E41A1C")
g1=ggplot(data = cdr3aa,mapping = aes(x =type_group ,y =cdr3_aa_length)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type)) +
  geom_boxplot(mapping = aes(fill = type),scale = "width",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons,aes(label=..p.signif..))+
  scale_fill_manual(values=col_flg)+ labs(title = "total_iso_CDR3_aa")+#facet_wrap(~type,scales = "free_y",ncol=3)+
  xlab("celltype") + ylab("CDR3 aa length") +#coord_cartesian(ylim = c(0, 0.25))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));g1

ggsave("FigS8G.pdf",width = 8,height = 8)

# FigS8H ------------------------------------------------------------------

cdr3aa2 <- cdr3aa[cdr3aa$type=="Cancer",]
a=table(cdr3aa2$cdr3_aa_length,cdr3aa2$celltype) %>% as.data.frame();head(a)

colnames(a)=c("cdr3_aa_length","celltype","Freq")
unique(a$celltype)

a1 <- a
p1<-ggplot(a1,aes(cdr3_aa_length,Freq,color=celltype,fill=celltype))+
  geom_histogram(stat = "identity")+ 
  labs(x="Samples",y="Freq")+
  geom_smooth(aes(cdr3_aa_length,Freq,color=celltype))+theme_classic();p1
for(i in 1:4){
  b=((a1[a1$celltype == unique(a1$celltype)[i],]$Freq)/sum(a1[a1$celltype == unique(a1$celltype)[i],]$Freq)) *100
}

Naive=((a1[a1$celltype == "B.01.TCL1A+naiveB",]$Freq)/sum(a1[a1$celltype == "B.01.TCL1A+naiveB",]$Freq)) *100
swbm=((a1[a1$celltype == "B.08.ITGB1+SwBm",]$Freq)/sum(a1[a1$celltype == "B.08.ITGB1+SwBm",]$Freq)) *100
atm=((a1[a1$celltype == "B.09.DUSP4+AtM",]$Freq)/sum(a1[a1$celltype == "B.09.DUSP4+AtM",]$Freq)) *100
DZ=((a1[a1$celltype == "B.11.SUGCT+DZ_GCB",]$Freq)/sum(a1[a1$celltype == "B.11.SUGCT+DZ_GCB",]$Freq)) *100
LZ=((a1[a1$celltype == "B.12.LMO2+LZ_GCB",]$Freq)/sum(a1[a1$celltype == "B.12.LMO2+LZ_GCB",]$Freq)) *100
CG=((a1[a1$celltype == "B.13.Cycling_GCB",]$Freq)/sum(a1[a1$celltype == "B.13.Cycling_GCB",]$Freq)) *100

a1$Freq1=c(Naive,swbm,atm,DZ,LZ,CG)
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
my36colors <- my36colors[c(1,8,9,11,12,13)]
p1<-ggplot(a1,aes(cdr3_aa_length,Freq1,fill=celltype))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  stat_summary(geom = "errorbar",fun.data = 'mean_sd', width = 0.3)+
  geom_density(alpha=.2) +
  scale_fill_manual(values=my36colors)+
  labs(x="Samples",y="Freq")+theme_classic2() + theme();p1

ggsave("FigS8H_CDR3 anmio length total cancer.pdf",width = 12,height = 6)
