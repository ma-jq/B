
rm(list = ls());gc()
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



my36colors <-c("BLCA"='#E5D2DD', "BRCA"='#53A85F', "CESC"='#F1BB72',"CHOL"= '#F3B1A0', "COAD"='#D6E7A3',"CTCL"= '#57C3F3',"DLBCL"= '#476D87',
               "ESCA"='#E95C59',"GBC"='#E59CC4', "GIST"='#AB3282', "HNSC"='#23452F',"LC"= '#BD956A',"LIHC"= '#8C549C',"NB"= '#585658',
               "OV"='#9FA3A8', "PAAD"='#E0D4CA', "RCC"='#5F3D69', "STAD"='#C5DEBA', "THCA"='#58A4C3',"THYM"= '#E4C755', "UCEC"='#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

load("./color.B.Rdata")

type=c( "Blood"="#4daf4a", "LN_Met"="#984ea3","Adjacent"="#377eb8","Cancer"="#e41a1c")
class=c("GC"="#377eb8" ,"EF"="#e41a1c")

##load sc-BCR data
load("db_obs230206.Rdata")

db_obs[1:4,1:4]
rownames(db_obs)<-db_obs$cell_id


###############################################################################
#'          Merge paired scRNA data with scBCR data                          '#
###############################################################################
objN <- readRDS("Pancancer_all_remove_use_final_dim18_20230921.rds")

Idents(objN)="dataid";table(Idents(objN))
object=subset(objN,idents = c("PCall","PCall_new"))
Idents(object) <- object$type
object <- subset(object,idents = c("Cancer","Blood","Adjacent","LN_Met"))
table(object$type)


object$patient=str_remove_all(object$patient,"PCall_")
table(object$patient)

####remove HCC62 HCC70 HCC136
Idents(object) <- object$patient
patient <- as.character(unique(Idents(object)))
patient.use <- setdiff(patient,c("HCC62","HCC70","HCC136"))
object=subset(object,idents = patient.use)
table(object$cancer)
Idents(object)="dataid";table(Idents(object))
#object=subset(object,idents = c("PCall","PCall_new"))

object@meta.data$cell_id <- rownames(object@meta.data)
object@meta.data$barcode=str_remove_all(object$cell_id,"PCall_")
object@meta.data$barcode=str_replace_all(object$barcode,"\\.","_")
object@meta.data$barcode=str_replace_all(object$barcode,"_3","")

object@meta.data$barcode=str_replace_all(object$barcode,"_1","")
object@meta.data$barcode=str_replace_all(object$barcode,"_2","")

object@meta.data$barcode=str_replace_all(object$barcode,"_T_","_")
object@meta.data$barcode=str_replace_all(object$barcode,"_P_","_")
object@meta.data$barcode=str_replace_all(object$barcode,"_B_","_")

object@meta.data$barcode= str_split(object$barcode,"_",simplify = T)[,2]

####standardize patient ID
object$patient=ifelse(object$patient %in% "BCRA1","BRCA1",object$patient)
object$patient=ifelse(object$patient %in% "CESC8B","CESC8",object$patient)
object$patient=ifelse(object$patient %in% "COAD14B","COAD14",object$patient)
object$patient=ifelse(object$patient %in% "COAD15B","COAD15",object$patient)
object$patient=ifelse(object$patient %in% "COAD16B","COAD16",object$patient)
object$patient=ifelse(object$patient %in% "GIST1B","GIST1",object$patient)
object$patient=ifelse(object$patient %in% "HCC31","HCC31T",object$patient)
object$patient=ifelse(object$patient %in% "HCC32","HCC32T",object$patient)

object$patient=ifelse(object$patient %in% c("HCC36B"),"HCC36",object$patient)
object$patient=ifelse(object$patient %in% "HCC38B","HCC38",object$patient)
object$patient=ifelse(object$patient %in% c("HCC39B"),"HCC39",object$patient)
object$patient=ifelse(object$patient %in% c("HCC40B"),"HCC40",object$patient)
object$patient=ifelse(object$patient %in% c("HCC41B"),"HCC41",object$patient)

object$patient=ifelse(object$patient %in% c("LC10B"),"LC10T",object$patient)
object$patient=ifelse(object$patient %in% c("LC11B"),"LC11",object$patient)

object$patient=ifelse(object$patient %in% c("PAAD2B"),"PAAD2",object$patient)
object$patient=ifelse(object$patient %in% c("STAD3B"),"STAD3",object$patient)

#object$patient=ifelse(object$patient %in% c("LC2T"),"LC2",object$patient)######直接和LC2T匹配好了

object@meta.data$BCR_id =paste0(object$patient,"_",object$barcode)
head(object@meta.data)

object$BCR_id=ifelse(object$patient %in% c("COAD21","GBC01","STAD10","STAD11","STAD7","THCA13" ) 
                     & object$type %in% c("LN_Met"), paste0(object$patient,"_LN_",object$barcode),
                     object$BCR_id)

object$BCR_id=ifelse(object$patient %in% c("STAD11" ) 
                     & object$type %in% c("Adjacent"), paste0(object$patient,"_LN_",object$barcode),
                     object$BCR_id)

object$BCR_id=ifelse(object$patient %in% c("THCA1","THCA10") 
                     & object$type %in% c("Blood"), paste0(object$patient,"_LN_",object$barcode),
                     object$BCR_id)


#db_obs$cell_id=rownames(db_obs)
head(str_split(db_obs$cell_id,"_",simplify = T))
table(str_split(db_obs$cell_id,"_",simplify = T)[,1])

head(table(str_split(db_obs$cell_id,"_",simplify = T)[,2]) %>% as.data.frame() %>% arrange(desc(Freq)))
db_obs$cell_id=str_replace_all(db_obs$cell_id,"_2","")

str_split(db_obs$cell_id,"_",simplify = T)[str_split(db_obs$cell_id,"_",simplify = T)[,2] == "2" ,][,1]###RCC5
str_split(db_obs$cell_id,"_",simplify = T)[str_split(db_obs$cell_id,"_",simplify = T)[,2] == "LN" ,][,1]###COAD21  GBC01 STAD10 STAD11  STAD7  THCA1 THCA10 THCA13 

###
table(db_obs[duplicated(db_obs$cell_id),]$patient)
table(db_obs[duplicated(db_obs$cell_id),]$cell_id)
db_obs=db_obs[!duplicated(db_obs$cell_id),]

###match paired BCR data
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
#db_obs[b[match(a$BCR_id,rownames(b)),],]
new<-cbind(a,d)


object_BCR<-subset(object, cells = row.names(subset(object@meta.data, object@meta.data$BCR_id %in% y)))
object_BCR=object_BCR[,!(duplicated(object_BCR$BCR_id))]
dim(object_BCR)
object=object_BCR
dim(object)
unique(object$type)
#object$type <- gsub("Lymph|LNN","LN_Met",object$type)

dim(object)
rownames(new)[1:3];colnames(object)[1:3]
object@meta.data<-new  


#load ASCs scRNA data & filter ASCs
ASC <- readRDS("ASC_rmtype_20230922.rds")
ASC_use_cellname <- colnames(ASC)

Idents(object) <- object$celltype_level3_0912
object.ASC.cellname <- WhichCells(object,idents="c15_MZB1+ASC")
object.del.cellname <- setdiff(object.ASC.cellname,ASC_use_cellname)
cellname.use <- setdiff(colnames(object),object.del.cellname)
object3 <- object[,cellname.use]
dim(object3)
dim(object)

###save object
saveRDS(object,file = "./BCR/BCR_add_20230915.rds")
saveRDS(object3,file = "./BCR/BCR_add_ASC_select_20230919.rds")


rownames(new)[1:3];colnames(object)[1:3]

