setwd(".")
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

library(Seurat)
library(ggplot2)
library(monocle3)
set.seed(12345)

scdata <- readRDS("./objTN_B_GCB230320.rds")
DimPlot(scdata,label=T,group.by = "celltype")

set.seed(123)
scdata_sample <- scdata[,sample(colnames(scdata),50000)]
dim(scdata_sample)
DimPlot(scdata_sample,label=T,group.by = "celltype")
Idents(scdata_sample) <- scdata_sample$celltype
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

cds <- cluster_cells(cds,resolution = 0.001,k=80,random_seed=18,verbose=T)
cds@clusters$UMAP$clusters <- scdata_sample$celltype
plot_cells(cds,color_cells_by = "partition")

p1 <- plot_cells(cds,group_cells_by = 'cluster')
p1

cds <- learn_graph(cds, verbose =T,
                   use_partition=T,close_loop=F)
p <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=T,
                label_leaves=T, label_branch_points=T,cell_size = 0.5,group_label_size=4)

p


cds <- order_cells(cds)

unique(cds@colData$celltype)
scdata_sample$celltype <- gsub("B_12_LMO2_LZGC","B_c12_LMO2_LZGC",scdata_sample$celltype)
scdata_sample$celltype <- gsub("B_11_CXCR4_DZGC","B_c11_CXCR4_DZGC",scdata_sample$celltype)
scdata_sample$celltype <- gsub("B_10_PSME2_PreGC","B_c10_PSME2_PreGC",scdata_sample$celltype)
scdata_sample$celltype <- gsub("B_02_IFIT3_B","B_c02_IFIT3_B",scdata_sample$celltype)
scdata_sample$celltype <- gsub("B_01_TCL1A_naïveB","B_c01_TCL1A_naïveB",scdata_sample$celltype)
p1 <- plot_cells(cds,color_cells_by = "pseudotime",label_branch_points = FALSE,label_leaves = F)
p2 <- DimPlot(scdata_sample,group.by = "celltype",cols = my36colors )
p2|p1
dir.create("monocle3")
ggsave("./monocle3/monocle3_0404.pdf",width = 25,height = 8)

save(cds,scdata_sample,file="./monocle3/monocle3_new.RData")



load("./monocle3/monocle3_new.RData")

Track_genes_sig <- c("APEX1","APEX2","XRCC5","XRCC6","POLD2","AICDA")
Track_genes_sig <- c("TCL1A","FGR","ITGAX","DUSP4")
Track_genes_sig <- c("TCL1A","ITGB1","CRIP1", "ITGAX",
                     "TBX21","GLS","GLA","PTGES3",
                     "BCL6", "S1PR2","SUGCT","AICDA","PRDM1",
                     "IRF4","IRF8","EBI2","CXCR3")
Track_genes_sig <- c("PDCD1","CD274","CTLA4","ENTPD1","LAG3",
                     "HAVCR2","TOX","TOX2","ZBED2","BATF","RBPJ","VDR",
                     "CD24","CD38","SELL")

genes <- c("TCL1A","FGR","ITGB1","SUGCT","BCL6","NR4A2")
genes <- c("TCL1A","IFI44L",
           "HSPA1A", "EGR1",
           "CD69","CCR7","ITGB1",
           "CRIP1","ITGAX","TBX21")
genes <- c("GLA","GLS","PTGES")

#基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)
#FeaturePlot图
pdf("./monocle3/featureplot_0316.pdf",width = 16,height = 10)
plot_cells(cds, genes=genes, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
dev.off()


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
df$celltype <- cds@colData$celltype_l3
#df <- df[order(df$pseudoTime),]
df2 <- df[sample(rownames(df),10000),]
df2$cellid <- rownames(df2)
df2 <- melt(df2,id.vars = c('pseudotime','cellid','celltype'))
df2$value <- as.numeric(df2$value)
colnames(df2)[4] <- "gene"

pdf("./monocle3/pseudotime_CSR_0404.pdf",width = 10,height = 8)
ggplot(df2, aes(x=pseudotime,y=value,color=gene))+
  geom_smooth(fullrange=T,method="loess",span=1)+ theme_classic()+
  ylab('Expression')
dev.off()

ggplot(df2, aes(x=pseudotime,y=value))+
geom_smooth(fullrange=T,method="loess",span=1)+ theme_classic()+
  ylab('Expression')+facet_wrap(~gene,scales="free")
ggsave("./monocle3/monocle3_gene_0406.pdf",width=20,height=15)

ggplot(df2, aes(x=pseudotime,y=value,color=celltype))+
  geom_smooth(method="loess",span=1)+ theme_classic()+
  ylab('Expression') +geom_jitter()+facet_wrap(~gene,scales="free")
#####
df2 <- data.frame(pseudotime=pseudotime(cds),
                  celltype=cds@colData$celltype_l3)
unique(df2$celltype)
df2$celltype <- gsub("B_01_TCL1A_naïveB","B_c01_TCL1A_naïveB",df2$celltype)
df2$celltype <- gsub("B_02_IFIT3_B","B_c02_IFIT3_B",df2$celltype)
df2$celltype <- gsub("B_10_PSME2_PreGC","B_c10_PSME2_PreGC",df2$celltype)
df2$celltype <- gsub("B_11_CXCR4_DZGC","B_c11_CXCR4_DZGC",df2$celltype)
df2$celltype <- gsub("B_12_LMO2_LZGC","B_c12_LMO2_LZGC",df2$celltype)

df2 <- df2[df2$celltype%in%c("B_c08_ITGB1_SwBm" ,"B_c09_DUSP4_AtM" , "B_c07_CCR7_ACB3", 
                             "B_c06_NR4A2_ACB2" ,"B_10_PSME2_PreGC","B_11_CXCR4_DZGC","B_12_LMO2_LZGC",
                             "B_01_TCL1A_naïveB", "B_c05_EGR1_ACB1" ),]
ggplot(df2,aes(y=celltype,x=pseudotime,color=celltype))+geom_boxplot()+
  theme_bw()+scale_color_manual(values = my36colors[c(1,5:12)])
ggplot(df2,aes(y=celltype,x=pseudotime,color=celltype))+geom_boxplot()+
  theme_bw()+scale_color_manual(values = my36colors[c(1:12)])
ggsave("./monocle3/boxplot_pseudotime_0404.pdf",width = 9,height = 6)


