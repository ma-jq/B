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

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

###############################################################################
#'                   Manuscipt: figure4F&G;S8K&L                             '#
###############################################################################

pbmc <- readRDS("Pancancer_all_remove_use_final_dim18_20230921.rds")
pbmc <- subset(pbmc,idents=c("c14_PB","c15_MZB1+ASC"),invert=TRUE)
##
set.seed(123)
cellname <- sample(colnames(pbmc),50000)
scdata_sample <- pbmc[,cellname]
scdata_sample <- Seurat::RunUMAP(scdata_sample,reduction = "harmony", dims = 1:10) 
scdata_sample <- FindNeighbors(scdata_sample,reduction = "harmony", dims = 1:13) 

DimPlot(scdata_sample,group.by = "celltype_level3_0912",label=T,cols = my36colors)
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
ggsave("monocle3_type2_20231005.pdf",width = 18,height = 8)


#####boxplot#####
df2 <- data.frame(pseudotime=pseudotime(cds),
                  celltype=cds@colData$celltype_level3_0912)
unique(df2$celltype)

df2 <- df2[df2$celltype%in%c("c08_ITGB1+SwBm" ,"c09_AtM" , "c07_CCR7+ACB", 
                             "c06_NR4A2+ACB" ,"c10_PreGC","c11_DZ GCB","c12_LZ GCB",
                             "c01_TCL1A+naiveB", "c05_EGR+ACB","c13_Cycling GCB" ),]
ggplot(df2,aes(y=celltype,x=pseudotime,color=celltype))+geom_boxplot()+
  theme_bw()+scale_color_manual(values = my36colors[c(1,5:13)])
ggplot(df2,aes(y=celltype,x=pseudotime,color=celltype))+geom_boxplot()+
  theme_bw()+scale_color_manual(values = my36colors[c(1:13)])
ggsave("./review/monocle3/boxplot_pseudotime_1004.pdf",width = 9,height = 6)

#######trajectory plot######
Track_genes_sig <- c("PDCD1","ENTPD1","HAVCR2",
                     "CD24","CD38","SELL")
gene2 <- c("XRCC6","XRCC5","APEX1","POLD2","AICDA","APEX2")
Track_genes_sig <- c(Track_genes_sig,gene2)
#Track_genes_sig <- c("GLS","GLA","PTGES3")
####
library(reshape2)
library(Seurat)
colData(cds)$pseudotime <- pseudotime(cds)
colnames(scdata_sample)[1:5]
colData(cds)[1:5]
data <- scdata_sample[Track_genes_sig,]
data <- GetAssayData(data@assays$RNA,slot='data')
data <- as.data.frame(data)

df <- t(data)
df <- as.data.frame(df)
df$pseudotime <- cds@colData$pseudotime 
df$celltype <- cds@colData$celltype_level3_0912
unique(df$celltype)

##
library(dplyr)
md <- as_tibble(cds@colData, rownames = NA) %>% select(celltype_level3_0912, pseudotime)
path1 <- c("c01_TCL1A+naiveB","c05_EGR+ACB","c06_NR4A2+ACB","c07_CCR7+ACB","c09_AtM")
path2 <- c("c01_TCL1A+naiveB","c05_EGR+ACB","c06_NR4A2+ACB","c07_CCR7+ACB","c08_ITGB1+SwBm",
           "c11_DZ GCB","c12_LZ GCB","c13_Cycling GCB")
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

totalMD.clean <- NULL
for(i in 1:length(unique(totalMD$celltype))){
  temp <- totalMD[totalMD$celltype==unique(totalMD$celltype)[i],]
  temp.clean <- remove_outlier(temp,c("pseudotime"))
  totalMD.clean <- rbind(totalMD.clean,temp.clean)
}

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


###############################################################################
#'                       Manuscipt: figureS8B                                '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
Idents(object) <- object$celltype_level3_0912
unique(object$celltype_level3_0912)
object <- subset(object,idents=c("c01_TCL1A+naiveB","c08_ITGB1+SwBm",
                                 "c09_AtM","c12_LZ GCB","c11_DZ GCB",
                                 "c13_Cycling GCB"))
object$celltype_level3_0912 <- as.character(object$celltype_level3_0912)
my36colors <-c("#3d7bb0" ,"#a2c8db" ,"#5aa554", "#afd195", "#d83f36", "#e99997" ,"#f2bc7c", "#e68740")
names(my36colors) <- c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4")

object$c_call <- factor(object$c_call,levels=c("IGHA1","IGHA2",
                                               "IGHG1","IGHG2","IGHG3","IGHG4",
                                               "IGHD","IGHM"
))

fig <- list()
for(i in 1:length(unique(object$celltype_level3_0912))){
  Idents(object) <- object$celltype_level3_0912
  temp <- subset(object,idents=unique(object$celltype_level3_0912)[i])
  table(temp$c_call)
  temp$type <- factor(temp$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
  temp$celltype_level3_0912 <- factor(temp$celltype_level3_0912
                                        ,levels=c("c01_TCL1A+naiveB","c08_ITGB1+SwBm",
                                                  "c09_AtM","c12_LZ GCB","c11_DZ GCB"))
  #temp <- subset(temp,idents=c("IGHD","IGHE",""),invert=TRUE)
  f5=dittoBarPlot(temp, "c_call",group.by="type",
                  retain.factor.levels = T,main = "type",color.panel= my36colors)+
    ylim(0,1)+ggtitle(unique(object$celltype_level3_0912)[i]);f5
  fig[[i]] <- f5
  
}
fig[['nrow']] <- 1
fig[['ncol']] <- 6
library(gridExtra)
pdf('FigS8B_new_20231012.pdf', width = 20, height = 8)
do.call('grid.arrange', fig)
dev.off()


###############################################################################
#'                       Manuscipt: figureS8C                                '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
Idents(object)="celltype_level3_0912";table(Idents(object))
B89=subset(object,idents=c("c09_AtM","c08_ITGB1+SwBm"))
table(B89$c_call)
B89$type_celltype <- paste0(B89$type,"-",B89$celltype_level3_0912)
Idents(B89) <- B89$c_call
unique(B89$c_call)
#B89 <- subset(B89,idents=c("IGHA1","IGHG3","IGHA2" ,"IGHG1", "IGHG2", "IGHM","IGHG4"))
B89$c_call <- factor(B89$c_call,levels = c("IGHA1","IGHA2",
                                           "IGHG1","IGHG2",
                                           "IGHG3","IGHM",
                                           "IGHD","IGHG4"
                                           ))
table(B89$c_call)
data <- as.data.frame(table(B89$type,as.character(B89$celltype_level3_0912),B89$c_call))
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
B89_pro$Var2=factor(B89_pro$Var2,levels = c("Blood-c08_ITGB1+SwBm","Blood-c09_AtM",
                                            "Adjacent-c08_ITGB1+SwBm" ,"Adjacent-c09_AtM",
                                            "LN_Met-c08_ITGB1+SwBm","LN_Met-c09_AtM" ,
                                            "Cancer-c08_ITGB1+SwBm","Cancer-c09_AtM"
))
#my_comparisons=list(c("c08_ITGB1+SwBm","c09_AtM"))
my_comparisons=list(c("Blood-c09_AtM","Blood-c08_ITGB1+SwBm"),
                      c("Adjacent-c09_AtM","Adjacent-c08_ITGB1+SwBm") ,
                      c("Cancer-c09_AtM" , "Cancer-c08_ITGB1+SwBm"),
                      c("LN_Met-c09_AtM" , "LN_Met-c08_ITGB1+SwBm")
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
ggsave("FigS8C_SwBm and AtM isotype comparison_20231010.pdf",width = 12,height = 8)

###############################################################################
#'                       Manuscipt: figureS8F                                '#
###############################################################################
object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
Idents(object) <- object$celltype_level3_0912
unique(object$celltype_level3_0912)
object <- subset(object,idents=c("c01_TCL1A+naiveB","c08_ITGB1+SwBm",
                                 "c09_AtM","c12_LZ GCB","c11_DZ GCB","c13_Cycling GCB"))
object$celltype_level3_0912 <- as.character(object$celltype_level3_0912)
my_comparisons=list(c("c08_ITGB1+SwBm_Blood","c09_AtM_Blood"),
                    c("c08_ITGB1+SwBm_Adjacent","c09_AtM_Adjacent"),
                    c("c08_ITGB1+SwBm_Cancer","c09_AtM_Cancer"),
                    c("c08_ITGB1+SwBm_LN_Met","c09_AtM_LN_Met"),
                    c("c09_AtM_Cancer","c11_DZ GCB_Cancer"),
                    c("c09_AtM_Cancer","c12_LZ GCB_Cancer"),
                   c("c09_AtM_Cancer","c13_Cycling GCB_Cancer"),
                   c("c09_AtM_Adjacent","c11_DZ GCB_Adjacent"),
                   c("c09_AtM_Adjacent","c12_LZ GCB_Adjacent"),
                   c("c09_AtM_Adjacent","c13_Cycling GCB_Adjacent"))
object$type_group <- paste0(object$celltype_level3_0912,"_",object$type)
object$type <- factor(object$type,levels = c("Blood","Adjacent","LN_Met","Cancer"))
unique(object$type_group)
type <- paste0("c01_TCL1A+naiveB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type2 <-  paste0("c08_ITGB1+SwBm","_",c("Blood","Adjacent","LN_Met","Cancer"))
type3 <-  paste0("c09_AtM","_",c("Blood","Adjacent","LN_Met","Cancer"))
type4 <-  paste0("c11_DZ GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type5 <-  paste0("c12_LZ GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type6 <-  paste0("c13_Cycling GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
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
ggsave("FigS8F_20231010.pdf",width = 8,height = 8)

###############################################################################
#'                      Manuscipt: figureS8G&H                               '#
###############################################################################

object <- readRDS("BCR_add_ASC_select_EF_GC_patient_level_rm_isotype_20231002.rds")
cdr3aa=aminoAcidProperties(object@meta.data, seq="cdr3", trim=TRUE,
                           label="cdr3")
# Define a ggplot theme for all plots
tmp_theme <- theme_bw() + theme(legend.position="bottom")
cdr3aa$type_group <- paste0(cdr3aa$celltype_level3_0912,"_",cdr3aa$type)
my_comparisons=list(c("c08_ITGB1+SwBm_Blood","c09_AtM_Blood"),
                    c("c08_ITGB1+SwBm_Adjacent","c09_AtM_Adjacent"),
                    c("c08_ITGB1+SwBm_Cancer","c09_AtM_Cancer"),
                    c("c08_ITGB1+SwBm_LN_Met","c09_AtM_LN_Met"),
                    c("c09_AtM_Cancer","c11_DZ GCB_Cancer"),
                    c("c09_AtM_Cancer","c12_LZ GCB_Cancer"),
                    c("c09_AtM_Cancer","c13_Cycling GCB_Cancer"),
                    c("c09_AtM_Adjacent","c11_DZ GCB_Adjacent"),
                    c("c09_AtM_Adjacent","c12_LZ GCB_Adjacent"),
                    c("c09_AtM_Adjacent","c13_Cycling GCB_Adjacent"))
type <- paste0("c01_TCL1A+naiveB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type2 <-  paste0("c08_ITGB1+SwBm","_",c("Blood","Adjacent","LN_Met","Cancer"))
type3 <-  paste0("c09_AtM","_",c("Blood","Adjacent","LN_Met","Cancer"))
type4 <-  paste0("c11_DZ GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type5 <-  paste0("c12_LZ GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type6 <-  paste0("c13_Cycling GCB","_",c("Blood","Adjacent","LN_Met","Cancer"))
type <- c(type,type2,type3,type4,type5,type6)
cdr3aa$type_group <- factor(cdr3aa$type_group,
                            levels = type)
cdr3aa$type <- factor(cdr3aa$type,levels =  c("Blood","Adjacent","LN_Met","Cancer"))
# Generate plots for all four of the properties
g1=ggplot(data = cdr3aa,mapping = aes(x =type_group ,y =cdr3_aa_length)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type)) +
  geom_boxplot(mapping = aes(fill = type),scale = "width",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons,aes(label=..p.signif..))+
  scale_fill_manual(values=col_flg)+ labs(title = "total_iso_CDR3_aa")+#facet_wrap(~type,scales = "free_y",ncol=3)+
  xlab("celltype") + ylab("CDR3 aa length") +#coord_cartesian(ylim = c(0, 0.25))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));g1

ggsave("FigS8G_20231010.pdf",width = 8,height = 8)

# FigS8H ------------------------------------------------------------------

cdr3aa2 <- cdr3aa[cdr3aa$type=="Cancer",]
a=table(cdr3aa2$cdr3_aa_length,cdr3aa2$celltype_level3_0912) %>% as.data.frame();head(a)

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

Naive=((a1[a1$celltype == "c01_TCL1A+naiveB",]$Freq)/sum(a1[a1$celltype == "c01_TCL1A+naiveB",]$Freq)) *100
swbm=((a1[a1$celltype == "c08_ITGB1+SwBm",]$Freq)/sum(a1[a1$celltype == "c08_ITGB1+SwBm",]$Freq)) *100
atm=((a1[a1$celltype == "c09_AtM",]$Freq)/sum(a1[a1$celltype == "c09_AtM",]$Freq)) *100
DZ=((a1[a1$celltype == "c11_DZ GCB",]$Freq)/sum(a1[a1$celltype == "c11_DZ GCB",]$Freq)) *100
LZ=((a1[a1$celltype == "c12_LZ GCB",]$Freq)/sum(a1[a1$celltype == "c12_LZ GCB",]$Freq)) *100
CG=((a1[a1$celltype == "c13_Cycling GCB",]$Freq)/sum(a1[a1$celltype == "c13_Cycling GCB",]$Freq)) *100

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

ggsave("FigS8H_CDR3 anmio length total cancer_20231010.pdf",width = 12,height = 6)
