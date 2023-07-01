setwd("./")

library(Seurat)

anno <- readRDS("./objTN2meta221121.rds")
data <- readRDS("./objTN221115.rds")
dim(data)
colnames(anno)

unique(anno$celltype_l7)

anno <- anno[intersect(rownames(anno),colnames(data)),]

table(anno$celltype_l7)
table(anno$cancer)

data <- data[,rownames(anno)]
data$typeuse <- paste0(anno$celltype_l7,"_",anno$cancer)
table(data$typeuse)

length(unique(data$typeuse))
cellanno <- data@meta.data
celluse <- NULL
set.seed(12345)
for(i in 1:length(unique(data$typeuse))){
  tmp <- cellanno[cellanno$typeuse==unique(data$typeuse)[i],]
  tmp.cell <- rownames(tmp)
  if(length(tmp.cell)>100){
    tmp.cell <- sample(tmp.cell,100)
  }
  else{
    tmp.cell <- tmp.cell
  }
  celluse <- c(celluse,tmp.cell)
}
length(celluse)

data.sample <- data[,celluse]
data.sample <- as.data.frame(GetAssayData(data.sample@assays$RNA,slot="counts"))
write.csv(data.sample,file="./SCENIC/resources/minicluster_counts_19294cell_1122.csv")

anno <- anno[celluse,c("celltype_l7","cancer","type")]
anno$cellID <- rownames(anno)
write.csv(anno,file="./SCENIC/resources/samplecluster_counts_19294cell_1122_anno.csv")





