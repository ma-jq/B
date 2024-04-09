

library(Seurat)
library(SeuratData)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(magrittr)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(clusterProfiler)

###############################################################################
#'                          Manuscipt: figure6A                             '#
###############################################################################
objN <- readRDS("../../scRNA_data/panB_scRNA_processed_data.rds")
Idents(objN)="celltype";table(Idents(objN))
objectTP=subset(objN,idents=c("B.08.ITGB1+SwBm",   "B.09.DUSP4+AtM"))
####
set.seed(100)
input.num = 10000
cellid<-sample(1:ncol(objectTP), input.num, replace=F); length(cellid)
objectTPN<-objectTP[,cellid]
dim(objectTPN)

#devtools::install_github("YosefLab/VISION@v2.1.0")

countexp.Seurat <- sc.metabolism.Seurat(obj = objectTPN, 
                                        method = "AUCell", # supports VISION, AUCell, ssgsea, and gsva, which VISION is the default method.
                                        imputation = F, ncores = 10, 
                                        metabolism.type = "KEGG") # supports KEGG and REACTOME, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.
saveRDS(countexp.Seurat,file = "./scMetabolism/countexp.Seurat for AtM and Bm.rds")



countexp.Seurat$type_celltype=paste0(countexp.Seurat$type,"_",countexp.Seurat$celltype_l3);table(countexp.Seurat$type_celltype)
countexp.Seurat$type_celltype=factor(countexp.Seurat$type_celltype,levels = c("Blood_B.08.ITGB1+SwBm" ,    "Blood_B.09.DUSP4+AtM",
                                                                              "Adjacent_B.08.ITGB1+SwBm",  "Adjacent_B.09.DUSP4+AtM",
                                                                              "LN_Met_B.08.ITGB1+SwBm" ,   "LN_Met_B.09.DUSP4+AtM",
                                                                              "Cancer_B.08.ITGB1+SwBm","Cancer_B.09.DUSP4+AtM"))
input.pathway<-c("D-Glutamine and D-glutamate metabolism")
input.pathway<-rownames(significance_pathway)
input.pathway<-significance_pathway
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "cancer", norm = "y")+
  scale_x_discrete(limits=as.character(c("Blood_B.08.ITGB1+SwBm" ,    "Blood_B.09.DUSP4+AtM",
                                         "Adjacent_B.08.ITGB1+SwBm",  "Adjacent_B.09.DUSP4+AtM",
                                         "LN_Met_B.08.ITGB1+SwBm" ,   "LN_Met_B.09.DUSP4+AtM",
                                         "Cancer_B.08.ITGB1+SwBm","Cancer_B.09.DUSP4+AtM")))+
  scale_y_discrete(limits=as.character(c("beta-Alanine metabolism","Pyrimidine metabolism" ,"Purine metabolism","Thiamine metabolism","Tyrosine metabolism","Phenylalanine metabolism",
                                         "Other types of O-glycan biosynthesis","Glycosaminoglycan biosynthesis - keratan sulfate" ,"N-Glycan biosynthesis",
                                         "Mannose type O-glycan biosynthesis" ,"Glycosphingolipid biosynthesis - lacto and neolacto series",
                                         "Ether lipid metabolism" ,"Glycerophospholipid metabolism","Terpenoid backbone biosynthesis","Other glycan degradation","Glycosphingolipid biosynthesis - ganglio series",
                                         "Folate biosynthesis" ,"Amino sugar and nucleotide sugar metabolism","Glycosaminoglycan degradation",
                                         "Glutathione metabolism",
                                         "Glycerolipid metabolism" ,
                                         "Riboflavin metabolism",
                                         "Oxidative phosphorylation","Pentose phosphate pathway","Inositol phosphate metabolism","Glycolysis / Gluconeogenesis" ,"Fructose and mannose metabolism",
                                         "Steroid hormone biosynthesis","Starch and sucrose metabolism" ,"Metabolism of xenobiotics by cytochrome P450" ,
                                         "Arachidonic acid metabolism",
                                         "Arginine biosynthesis","Alanine, aspartate and glutamate metabolism","Galactose metabolism","Sphingolipid metabolism","D-Glutamine and D-glutamate metabolism"
  )))


###############

input.pathway<-c("beta-Alanine metabolism","Pyrimidine metabolism" ,"Purine metabolism","Thiamine metabolism","Tyrosine metabolism","Phenylalanine metabolism",
                 "Other types of O-glycan biosynthesis","Glycosaminoglycan biosynthesis - keratan sulfate" ,"N-Glycan biosynthesis","Mucin type O-glycan biosynthesis",
                 "Galactose metabolism",  "Mannose type O-glycan biosynthesis" ,"Glycosphingolipid biosynthesis - lacto and neolacto series",
                 "Ether lipid metabolism" ,"Steroid hormone biosynthesis","Starch and sucrose metabolism" ,"Metabolism of xenobiotics by cytochrome P450" ,
                 "Arachidonic acid metabolism","Terpenoid backbone biosynthesis","Other glycan degradation","Glycosphingolipid biosynthesis - ganglio series","Glycosphingolipid biosynthesis - globo and isoglobo series",
                 "Amino sugar and nucleotide sugar metabolism","Arginine biosynthesis","Alanine, aspartate and glutamate metabolism","Glycerophospholipid metabolism","Glycerolipid metabolism" ,
                 "Glycosaminoglycan degradation","Folate biosynthesis" ,"Glutathione metabolism","Riboflavin metabolism",
                 "Oxidative phosphorylation","Pentose phosphate pathway","Inositol phosphate metabolism","Glycolysis / Gluconeogenesis" ,"Fructose and mannose metabolism",
                 "Sphingolipid metabolism","D-Glutamine and D-glutamate metabolism"
)
top10$gene
# 
kegg_geneset <- clusterProfiler::read.gmt("~/R/x86_64-pc-linux-gnu-library/4.2/scMetabolism/data/KEGG_metabolism_nc.gmt")
Reac_geneset <- clusterProfiler::read.gmt("~/R/x86_64-pc-linux-gnu-library/4.2/scMetabolism/data/REACTOME_metabolism.gmt")

# 
kegg_gene <- subset(kegg_geneset, term %in% input.pathway) %>% .$gene %>% as.vector(.)

reac_gene <- subset(Reac_geneset, term %in% top10$gene) %>% .$gene %>% as.vector(.)

meta_gene=c(kegg_gene,reac_gene)
meta_gene=meta_gene[!duplicated(meta_gene)]


####
obj_random2w$type_stype=paste0(obj_random2w$type,"_",obj_random2w$stype)
Idents(objectTPN) = "type_celltype";table(Idents(objectTPN))

###
library(future)
plan("multicore", workers =1)
options(future.globals.maxSize = 10000 * 1024^2)

diffpath_gene= Seurat::FindAllMarkers(objectTPN, only.pos = TRUE, min.pct = 0,logfc.threshold = 0)
write.csv(diffpath_gene,file = "./scMetabolism/diffpath_gene.csv")

diffpath_gene_meta = subset(diffpath_gene, diffpath_gene$gene %in% kegg_gene)


top10 <- diffpath_gene_meta%>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10=top10[!duplicated(top10$gene),]

vag_exp=AverageExpression(objectTPN,assays = "RNA",features = top10$gene,#top10$gene,
                          group.by = "type_celltype",slot="data")
vag_exp1=vag_exp$RNA

vag_exp2 <- as.data.frame(t(apply(vag_exp1,1,scale)))
colnames(vag_exp2) <- colnames(vag_exp1)

vag_exp2[vag_exp2 >=2]=2

vag_exp3=vag_exp2[,c("Blood_B.08.ITGB1+SwBm" ,    "Blood_B.09.DUSP4+AtM",
                     "Adjacent_B.08.ITGB1+SwBm",  "Adjacent_B.09.DUSP4+AtM",
                     "LN_Met_B.08.ITGB1+SwBm" ,   "LN_Met_B.09.DUSP4+AtM",
                     "Cancer_B.08.ITGB1+SwBm","Cancer_B.09.DUSP4+AtM")]

vag_exp4=vag_exp3[c("GLS" ,"GLA" ,    "UGCG",    "B4GALT1",    "PLD4" ,   "PIKFYVE", "IDI1",    "INPP5F",  "ASAH1",   "ACP5",
                    "PTGES3",  "IDS" ,    "ODC1",    "PDE4D","LPIN1",   "CERS4" ,  "POLD4",   "MGAT4A",  "NT5E",    "PIP4K2A",
                    "GPX4",    "TYMP","GAPDH",   "ATP6V1F", "ENO1",    "TPI1",    "ATP6V0C", "NDUFS5",  "PKM",
                    "ALPL" ,   "MGAT5",   "POLR2J3", "PNP",     "CHKB" ,   "DDOST" ,  "NDUFA6",
                    "GPCPD1",  "ENTPD1" , "HEXA",    "SYNJ2" ,  "NEU1",
                    "GPX1" ,      "GDA",        "COX7C" ,  "SGMS2",   "ALDOB",  
                    "PLPP3",   "CHST2","DGKD",    "GSTP1",   "MANBA",   "PFKL",    "CHPT1",
                    "AK8","COX4I1",  "ALOX5", "NME2", "COX6C", "TKT", "ATP5PB",  "NDUFB10", "PLPP5",   "UQCRFS1"
),]

vag_exp3=vag_exp2[,c("Blood_EF","Adjacent_EF", "LN_Met_EF" ,"Cancer_EF", "Blood_GC", "Adjacent_GC"   ,"LN_Met_GC",  "Cancer_GC"   )]
vag_exp3=vag_exp2[,c("Blood_EF", "Blood_GC","Adjacent_EF","Adjacent_GC"   , "LN_Met_EF" ,"LN_Met_GC", "Cancer_EF",  "Cancer_GC"   )]
pdf("./metabolism/EF and GC rASC Metabolism gene.pdf",width = 12,height = 20)

pdf("./Metabolism/AtM vsSwBm Metabolism gene230420.pdf",width = 8,height = 10)
p1= pheatmap::pheatmap(vag_exp4,cutree_cols = 5 ,
                       cluster_cols = F,cluster_rows = F, 
                       color=viridis::magma(20),
                       #breaks=seq(-2,2,0.1),
                       #color= scale_colour_gradientn(colors=brewer.pal(9, "BuPu")),
                       #color = viridis::plasma(40),
                       #color = colorRampPalette(c("navy", "white", "firebrick3"))(40),
                       #color = colorRampPalette(c("#ff59ff", "#080704", "#fffb06"))(40), #seurat color
                       #color = c(colorRampPalette(colors = c("#ff59ff","#080704"))(length(40)/2),colorRampPalette(colors = c("#080704","#fffb06"))(length(40)/2)),
                       #colorRampPalette(brewer.pal(9, "RdYIBu"))(50),
                       border=FALSE,#
                       cellwidth = 20,
                       cellheight = 8,
                       fontsize = 8,
                       gaps_col = c(1,2,3,4,5,6,7,8#,9,10,11#,12,
                                    #13,14,15,16,17,18,19,20,
                                    #21,22,23,24,25,26,27,28,29,30
                       )
);p1
dev.off()



#########cancer type######

input.pathway<-c("D-Glutamine and D-glutamate metabolism","Sphingolipid metabolism","Arachidonic acid metabolism")
input.pathway<-rownames(significance_pathway)
input.pathway<-significance_pathway
countexp.Seurat$cancer_celltype=paste0(countexp.Seurat$cancer,"-",countexp.Seurat$celltype_l3)
countexp.SeuratT=countexp.Seurat[,countexp.Seurat$type %in% "Cancer"]
DotPlot.metabolism(obj = countexp.SeuratT, pathway = input.pathway, phenotype = "cancer_celltype", norm = "y")

input.pathway<-c("Sphingolipid metabolism")
input.pathway<-c("Arachidonic acid metabolism")

DefaultAssay(countexp.SeuratT)="RNA"
VlnPlot(countexp.SeuratT,features = input.pathway,group.by = "cancer_celltype")
a=BoxPlot.metabolism(obj = countexp.SeuratT, pathway = input.pathway, phenotype = "cancer_celltype",ncol=1);a# +scale_fill_manual(values=c("B_c08_ITGB1_SwBm"='#E95C59', "B_c09_DUSP4_AtM"='#E59CC4') )

b=a$data
b$cancer_celltype=b$X1
b$celltype=str_split(b$X1,"-",simplify = T)[,2]

b$X2=factor(b$X2,levels = c("D-Glutamine and D-glutamate metabolism","Sphingolipid metabolism","Arachidonic acid metabolism"))
p1=ggplot(data = b,mapping = aes(x =cancer_celltype,y =X3)) +
  #geom_boxplot(scale = "width",adjust =1,trim = TRUE, mapping = aes(fill = type))  + 
  geom_boxplot(mapping = aes(fill = celltype),scale = "width",outlier.shape = NA)+
  #stat_compare_means(comparisons = my_comparisons)+
  #scale_fill_npg()+
  scale_fill_manual(values=c("B.08.ITGB1+SwBm"='#E95C59', "B.09.DUSP4+AtM"='#E59CC4') )+
  labs(title = "D-Glutamine and D-glutamate metabolism")+facet_wrap(~X2,scales = "free_y",ncol=1)+
  xlab("celltype") + ylab("score") +#coord_cartesian(ylim = c(0, 0.3))+
  theme_classic2() + theme(#legend.position = "none",
    axis.text.x =element_text(angle=45, hjust=1, vjust=1));p1

ggsave("./Metabolism/AtM vs Bm by cancer.pdf",width = 10,height = 10)






