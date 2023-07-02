

################fig.S3D ROGUE################
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(ggpubr))

objectnPC=readRDS("./PC_B_arranged221031/combined_BCR/mjq/docker_final230201/objectnPC230525.rds")
input.num = 100000
cellid<-sample(1:ncol(objN), input.num, replace=F); length(cellid)
obj_random2w<-objN[,cellid]
dim(obj_random2w)
obj.2_expr=objectnPC@assays$RNA@counts
metadata=objectnPC@meta.data

ent.res <- SE_fun(obj.2_expr)
head(ent.res)
SEplot(ent.res)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

rogue.res <- rogue(obj.2_expr, labels = metadata$groupn, samples = metadata$patient, platform = "UMI", span = 0.8)
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
p=ggboxplot(rogue.res2,x="variable",y="ROGUE",
          color = "variable",palette = class,
          add = "jitter")+ theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ 
  stat_compare_means(comparisons = my_comparisons)+
  coord_cartesian(ylim = c(0, 1));p





