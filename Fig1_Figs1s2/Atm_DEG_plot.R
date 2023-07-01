setwd(".")

library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggrepel)
library(patchwork)
#Atm
file <- list.files("./diffgene/",pattern = "_0317.csv")
file
Atm <- NULL
for(i in 1:length(file)){
  temp <- read.csv(paste0("./diffgene/",file[i]),row.names = 1)
  temp$compare <- colsplit(file[i],"marker",names=c("n1","n2"))$n1
  temp$compare <- colsplit(temp$compare,"_VS_",names = c("n1","n2"))$n2
  temp$gene <- rownames(temp)
  temp <- temp[-grep("^MT-",temp$gene),]
  Atm <- rbind(Atm,temp)
}
Atm$cluster <- Atm$compare
Atm$sig=""
Atm$sig[abs(Atm$avg_log2FC) > 0.25 & Atm$p_val_adj < 0.05] = "sig"
Atm$sig2=paste(Atm$cluster,Atm$sig,sep = "_")
Atm$sig2[str_detect(Atm$sig2,"_$")]="not_sig"
Atm$sig2=str_replace(Atm$sig2,"_sig","")

#控制顺序
Atm$sig2=factor(Atm$sig2,levels = c("not",sort(unique(Atm$cluster))))
Atm$cluster=factor(Atm$cluster,levels = sort(unique(Atm$cluster)))
Atm=Atm%>%arrange(cluster,sig2)
#控制范围
max(Atm$avg_log2FC);min(Atm$avg_log2FC)
Atm$avg_log2FC[Atm$avg_log2FC > 4]=4
Atm$avg_log2FC[Atm$avg_log2FC < c(-4)]= -4

#配色
library(RColorBrewer)
library(scales)
color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1])
names(color_ct)=sort(unique(as.character(Atm$cluster)))

#画图
Atm %>% ggplot(aes(x=cluster,y=avg_log2FC,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
  scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
  scale_y_continuous("average log2FC",expand = c(0.02,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.text.y.left = element_text(size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_text(size = 16)
  )

Atm2=Atm
Atm2$padj_log10_neg= -log10(Atm2$p_val_adj)
Atm2$padj_log10_neg=ifelse(Atm2$avg_log2FC > 0,
                           Atm2$padj_log10_neg,
                           -Atm2$padj_log10_neg)

Atm2$label <- ifelse(abs(Atm2$avg_log2FC)>0.8,Atm2$gene,"")

Atm2 %>% ggplot(aes(x=cluster,y=avg_log2FC,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
  scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
  scale_y_continuous("average log2FC",expand = c(0.02,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.text.y.left = element_text(size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_text(size = 16)
  )+geom_label_repel( aes(label = label),
                    size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE,
                    max.overlaps =50)

###################
plot.list=list()
for (ci in sort(unique(as.character(Atm2$cluster)))) {
  tmpdf=Atm2 %>% filter(cluster == ci)
  minabs=abs(min(tmpdf$padj_log10_neg))
  maxabs=max(tmpdf$padj_log10_neg)
  thre=0
  if(minabs < maxabs) {
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > minabs] = minabs
    thre=minabs
  }
  if(minabs > maxabs) {
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-maxabs)] = -maxabs
    thre=maxabs
  }
  if(minabs == maxabs & maxabs == Inf) {
    thre = min(
      abs(
        range(
          tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < Inf & tmpdf$padj_log10_neg > -Inf]
        )
      )
    )
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-thre)] = -thre
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > thre] = thre
  }

  plotdata = tmpdf
  tmpdf=tmpdf%>%filter(sig2 != "not") 
  tmpdf=tmpdf%>%arrange(desc(avg_log2FC))
  tmpdf.a=head(tmpdf%>%filter(avg_log2FC > 0),15)
  tmpdf.a$d=thre*2*0.05+(-thre)-tmpdf.a$padj_log10_neg
  tmpdf.b=tail(tmpdf%>%filter(avg_log2FC < 0),15)
  tmpdf.b$d=thre*2*0.95-thre  - tmpdf.b$padj_log10_neg
  textdata.down = tmpdf.b
  textdata.up   = tmpdf.a

  ###画图
  tmpplot=plotdata%>%ggplot(aes(x=padj_log10_neg,y=avg_log2FC))+
    geom_point(aes(color=sig2),size=1)+
    geom_hline(yintercept = c(-0.25,0.25),linetype="dashed")+
    geom_text_repel(data = textdata.down,
                    mapping = aes(label=gene),
                    nudge_x=textdata.down$d,
                    direction = "y", hjust = 1,segment.size = 0.2)+
    geom_text_repel(data = textdata.up,
                    mapping = aes(label=gene),
                    nudge_x=textdata.up$d,
                    direction = "y", hjust = 0,segment.size = 0.2)+
    labs(title = ci)+
    scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
    scale_y_continuous("average log2FC",expand = c(0.02,0),limits = c(-2,2))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none",

      axis.ticks.x.bottom = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = 14,color = "black"),
      axis.title.y.left = element_text(size = 16),

      plot.title = element_text(size = 12,hjust = 0.5)
    )

  index=which(ci == sort(unique(as.character(Atm2$cluster))))
  if (index!=1) {
    tmpplot=tmpplot+theme(
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank()
    )
  }
  if (index == length(sort(unique(as.character(Atm2$cluster))))) {
    segment.df=data.frame(x=c(0 - thre / 5,0 + thre / 5),
                          xend=c(-thre,thre),
                          y=c(-3,-3),
                          yend=c(-3,-3))
    tmpplot=tmpplot+geom_segment(data = segment.df,
                                 mapping = aes(x=x,xend=xend,y=y,yend=yend),
                                 arrow = arrow(length=unit(0.3, "cm")))

  }
  plot.list[[get("index")]]=tmpplot
}
wrap_plots(plot.list,ncol = 2)&theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
ggsave("./diffgene/Atm_diffgene_0317.pdf",width = 8,height = 8)
