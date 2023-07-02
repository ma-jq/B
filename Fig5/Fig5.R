

#####################
Idents(cd45)="celltype_l4";table(Idents(cd45))
cd45=subset(cd45,idents=c("B","Doublet"),invert=TRUE)

set.seed(100)
input.num = 10000
cellid<-sample(1:ncol(cd45), input.num, replace=F); length(cellid)
cd45N<-cd45[,cellid]
dim(cd45N)
saveRDS(cd45N,file="/home/data/vip13t45/project/PC/PC_scRNA_B/CD45/CD45/cd45N 1w 20230526.rds")


write.table(as.matrix(cd45N@assays$RNA@data), './cellphoneDB/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(cd45N@meta.data), cd45N@meta.data[,'celltype_l4', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA

write.table(meta_data, './cellphoneDB/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

#############linux##########
cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt      --counts-data=gene_name  --database out/cellphonedb_user_2023-05-24-08_46.db 
cellphonedb plot dot_plot --means-path ./out/means.txt --pvalues-path ./out/pvalues.txt --output-path ./Outplot

#######
sel_pval = all_pval[match(selected_pairs, intr_pairs), selected_celltype]#将上述需要呈现的受配体选出来
sel_means = all_means[match(selected_pairs, intr_pairs), selected_celltype]
df_names = expand.grid(selected_pairs, selected_celltype)
pval = unlist(sel_pval)
pval[pval==0] = 0.00001
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means))
plot.data = cbind(plot.data,pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
plot.data$clusters <- gsub('[|]', '_', plot.data$clusters)

#-------------------------------------------------------------------------------------------

#设置颜色+作图
plot.data$clusters=factor(plot.data$clusters,levels = c("B_c08_ITGB1_SwBm_CD8_IFI","B_c09_DUSP4_AtM_CD8_IFI",
                                                        "B_c08_ITGB1_SwBm_CD4_FOXP3", "B_c09_DUSP4_AtM_CD4_FOXP3", 
                                                        "B_c08_ITGB1_SwBm_CD8_CTLA4", "B_c09_DUSP4_AtM_CD8_CTLA4",
                                                        "B_c08_ITGB1_SwBm_CD4_PDCD1","B_c09_DUSP4_AtM_CD4_PDCD1"))

plot.data$meann=plot.data$mean
plot.data$meann=ifelse(plot.data$meann >0.5,0.5,plot.data$meann)
my_palette <- colorRampPalette(c("#fef0d9","#fc8d59","#b30000"))(n=10)
ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue + 0.00001),color=meann)) +
  scale_size_continuous(range=c(1,5))+
  scale_color_gradientn("log2(mean)", colors=my_palette,limits=c(0,0.5)) +
  coord_flip()+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size=10),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
ggsave("./cellphoneDB/Figure/dotplot.pdf",width = 10,height = 5)



##########nichenet

library(nichenetr)
library(Seurat)
library(tidyverse)

ligand_target_matrix <- readRDS("../PC_scRNA_B//Nichenet/ligand_target_matrix.rds")
lr_network = readRDS("../PC_scRNA_B//Nichenet/lr_network.rds")
weighted_networks = readRDS("../PC_scRNA_B//Nichenet/weighted_networks.rds")
scRNA <- UpdateSeuratObject(cd45)
head(scRNA)
head(scRNA@meta.data)
Idents(scRNA) <- "celltype_l4"
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = scRNA, 
                                               top_n_ligands = 50,
                                               receiver = c("B_c08_ITGB1_SwBm"), 
                                               sender = c( "CD4_PDCD1"),
                                               condition_colname = "type", 
                                               condition_oi = "Cancer", 
                                               condition_reference = "Adjacent", 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human")

nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = scRNA, 
                                               top_n_ligands = 50,
                                               receiver = c("B_c09_DUSP4_AtM"), 
                                               sender = c( "CD4_PDCD1"),
                                               condition_colname = "type", 
                                               condition_oi = "Cancer", 
                                               condition_reference = "Adjacent", 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human")


###########ST

library(glue, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

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
library(stCancer)
library(grid)


setwd("/www/data/PC")

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)



###########
#batch check
files = list.files(path = "/www/data/PC/publisheddata/treatedST/")
files = files[substr(files, 1,2) == "ST"]
cat(files)

all.sample = files

file.sum = c()
for (j in 1:length(all.sample)){
  
  print(j)
  
  input.sample = all.sample[j]
  if (!file.exists(paste("/www/data/PC/publisheddata/treatedST/", input.sample, "/data/spatial_anno.RDS", sep=""))) {
    print(paste(input.sample, " not exist!!!!!!!!", sep="")); 
    file.sum = rbind(file.sum, c(input.sample, NA, paste("/www/data/PC/publisheddata/treatedST/", input.sample, "/data/spatial_anno.RDS", sep="")))
    next
  }
  
  file.sum = rbind(file.sum, c(input.sample, "YES", paste("/www/data/PC/publisheddata/treatedST/", input.sample, "/data/spatial_anno.RDS", sep="")))
  
}

markerlist_integrate = list()

# #define PD1hiCD4 signature
load(file =  "/www/data/PC/PC_B/treateddata/markerlisttcell_ssgsea_cosg_10000.rda")
markerlist_integrate[["CD4_PDCD1_tcell"]] = markerlist[["CD4_PDCD1"]]

markerlist_integrate = list()
markerlist_integrate[["CD4_PDCD1"]] = markerlist[["CD4_PDCD1"]]

#define AtM marker
load("/home/data/data/PC/PC_B/treateddata/markerlist_ssgsea_cosg_10000.rda")
markerlist_integrate[["B_c09_DUSP4_AtM"]] = markerlist[["B_c09_DUSP4_AtM"]]

#TLS signature
signature <- qusage::read.gmt("/www/data/signature/tls.gmt"); names(signature)
markerlist_integrate[["TLS.NC"]] = signature[["TLS.NC"]]
markerlist_integrate[["TLS.Immunity"]] = signature[["TLS.Immunity"]]

markerlist = markerlist_integrate

for (k in 1:length(markerlist)){
  markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "RPL"]
  markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "RPS"]
  markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "HSP"]
  markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "MT-"]
  markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "EEF"]
  markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "HIST"]
  # markerlist[[k]] = markerlist[[k]][!markerlist[[k]] %like% "CD79"]
}


input.num = 500
for (k in 1:length(markerlist)){
  markerlist[[k]] = markerlist[[k]][1:input.num]
}


list_Bssgsea = list()
# list_xCell = list()


for (j in 1:length(all.sample)){
  # for (j in 1:10){
  
  print(j)
  t1 = Sys.time()
  
  obj_ST = readRDS(file.sum[j,3])
  print(dim(obj_ST))
  
  input.mat = as.matrix(obj_ST@assays$SCT@data)
  print(dim(input.mat))  
  
  
  library(xCell)
  PC_Bssgsea <- GSVA::gsva(as.matrix(input.mat), markerlist, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz = 20) #
  # PC_xCell = xCell::xCellAnalysis(as.matrix(input.mat), parallel.sz = 30)
  
  
  list_Bssgsea[[j]] = PC_Bssgsea
  # list_xCell[[j]] = PC_xCell
  
  print(Sys.time() - t1)
  print("\n")
}

save(list_Bssgsea, file = "/home/data/data/PC/PC_B/treateddata/list_ST_B.CD4T.TLS_Bssgsea_500.rda") 


list_PC_immunemat = list()
sum_B = c()
for (j in 1:length(all.sample)){

  PC_immunemat = t(data.frame(list_Bssgsea[[j]]))
  
  PC_immunemat = round(PC_immunemat, 3)
  
  cor_mat = cor(PC_immunemat, method = "spearman")

  
  list_PC_immunemat[[j]] = cor_mat
  
  tmp = data.frame(cor_mat[,2])
  colnames(tmp) = "spearman"
  tmp$celltype = row.names(tmp)
  tmp$sampleid = j
  tmp$samplename = all.sample[j]
  
  
  sum_B = data.frame(rbind(sum_B, tmp))
}


median_sum = c()
for (k in 1:length(unique(sum_B$celltype))){
  tmp = subset(sum_B, sum_B$celltype == unique(sum_B$celltype)[k])
  tmp2 = c(unique(sum_B$celltype)[k], median(tmp$spearman))
  
  median_sum = data.frame(rbind(median_sum, tmp2))
}
median_sum

sum_B.cd4 = subset(sum_B, sum_B$celltype == "B_c09_DUSP4_AtM")
quantile(sum_B.cd4$spearman)

input.study = all.sample
input.sample = c()
for (j in 1:length(input.study)){
  input.sample = c(input.sample, all.sample[all.sample %like% input.study[j]])
}

sum_B_select = subset(sum_B, sum_B$samplename %in% input.sample)
sum_B_select = subset(sum_B_select, !sum_B_select$celltype %like% "AtM")


sum_B_select <- sum_B_select[order(sum_B_select$spearman, decreasing = T),]
sum_B_select$celltype <- factor( sum_B_select$celltype, levels= unique(as.character(sum_B_select$celltype)))

ggplot(sum_B_select , aes(x = reorder(celltype, spearman, FUN=median, y = spearman), y=spearman,fill = celltype) )+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=my36colors) +
  geom_hline(yintercept=0.3, linetype="dashed", color = "red") +
  theme_bw()+theme(aspect.ratio=2, axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none", legend.box = "none") +
  # theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # facet_wrap(~label, ncol = 4, scales = "free") +
  # stat_compare_means(label = "p.format",method = "wilcox.test",
  #                    comparisons = list(c("Adjacent/Healthy", "Cancer"),
  #                                       c("Infection-COVID", "Cancer"),
  #                                       # c("Adjacent/Healthy", "Cancer"),
  #                                       c("PBMC", "Cancer"))) +
  NULL



#number of spot

nspot = 0

for (j in 1:length(all.sample)){
  # for (j in 1:10){
  
  if (!file.sum[j,1] %in% input.sample) next
  
  print(j)
  obj_ST = readRDS(file.sum[j,3])
  
  nspot = nspot + ncol(obj_ST)
  print(nspot)
}


#spatial plot

nspot = 0

for (j in 1:length(all.sample)){
  # for (j in 1:10){
  
  if (!file.sum[j,1] %in% input.sample) next
  
  print(j)
  t1 = Sys.time()
  
  obj_ST = readRDS(file.sum[j,3])
  print(dim(obj_ST))
  
  
  PC_mat = list_Bssgsea[[j]]
  
  obj_ST@meta.data = cbind(obj_ST@meta.data, t(PC_mat))
  
  input.features = row.names(PC_mat)
  p = SpatialPlot(obj_ST, features = input.features,  pt.size.factor = 2, image.alpha = 0, ncol = 3, 
                  crop = TRUE, stroke = NA) +#& scale_fill_gradientn(colours = (cc(100))) #+
    SpatialPlot(obj_ST,  pt.size.factor = 0, image.alpha = 1, #group.by = "SCT_snn_res.0.4",
                crop = TRUE, stroke = NA)
  
  ggsave(paste("/home/data/data/PC/PC_B/treateddata/treatedST/",j,"_", all.sample[j], ".png", sep=""),p,  width = 3000, height = 3000, units = "px")
  
  
  
}


selectj = 48
j = selectj
obj_ST = readRDS(file.sum[j,3])

PC_mat = list_Bssgsea[[j]]
obj_ST@meta.data = cbind(obj_ST@meta.data, t(PC_mat))
input.features = row.names(PC_mat)


p = SpatialPlot(obj_ST, features = input.features,  pt.size.factor = 2, image.alpha = 0, ncol = 3, 
                crop = TRUE, stroke = NA) +#& scale_fill_gradientn(colours = (cc(100))) #+
  SpatialPlot(obj_ST,  pt.size.factor = 0, image.alpha = 1, #group.by = "SCT_snn_res.0.4",
              crop = TRUE, stroke = NA)
p




