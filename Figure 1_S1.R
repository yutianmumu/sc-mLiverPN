---
title: "Figure 1__Figure S1"
date: '2023-06-10'
author: "Lin LEI"
---
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(brew)
library(stringr)
library(DT)
library(tidyverse)
library(patchwork)
library(harmony)
library(grid)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(future)
library(ggpubr)
set.seed(123)
setwd("~/Figure1_new")
#load the integrated (harmony) and filtered datasets including all the NCD (n=3), CDAHFD (n=3), and BDL (n=2) samples
integration_new <- readRDS("~/Figure1_new/integration_new.rds")
#UMAP plot
ToUmap <-DimPlot(integration_new, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)+ scale_colour_manual( values = colours)+NoLegend()
ggsave(filename = "ToUmap.pdf", plot = ToUmap, device="pdf", width=5.5, height=5)
#UMAP splitted by treatment
ToUmapsplit <- DimPlot(integration_new, reduction = "umap", label = TRUE,split.by = "treatment", repel = TRUE,pt.size = 0.05)+ scale_colour_manual( values = colours)+NoLegend()
ggsave(filename = "ToUmapsplit.pdf", plot = ToUmapsplit, device="pdf", width=15, height=5)
#hight the different conditions in total UMAP
integration_new<-SetIdent(integration_new,value=integration_new@meta.data$treatment)
hi1<-WhichCells(integration_new,idents = "Control")
hi2<-WhichCells(integration_new,idents = "CDAHFD")
hi3<-WhichCells(integration_new,idents = "BDL")
Control_UMAP<-DimPlot(integration_new,group.by = "treatment", reduction = "umap",cells.highlight= hi1, cols.highlight = c("#737373"),sizes.highlight = 0.05, label = F,repel = F,  pt.size = 0.05,cols = c("#BDBDBD"))
ggsave(filename = "Control_UMAP.pdf", plot = sham_UMAP, device="pdf", width=6.2, height=4.5)

CDAHFD_UMAP<-DimPlot(integration_new,group.by = "treatment", reduction = "umap",cells.highlight= hi2, cols.highlight = c("#2171B5"),sizes.highlight = 0.05, label = F,repel = F,  pt.size = 0.05,cols = c("#BDBDBD"))
ggsave(filename = "CDAHFD_UMAP.pdf", plot = CDAHFD_UMAP, device="pdf", width=6.2, height=4.5)

BDL_UMAP<-DimPlot(integration_new,group.by = "treatment", reduction = "umap",cells.highlight= hi3, cols.highlight = c("#238B45"),sizes.highlight = 0.05, label = F,repel = F,  pt.size = 0.05,cols = c("#BDBDBD"))
ggsave(filename = "BDL_UMAP.pdf", plot = BDL_UMAP, device="pdf", width=6.2, height=4.5)

#Plot cell density onto UMAP
data <- cbind(Embeddings(object=Control[['umap']]),FetchData(Control,'treatment'))
p1 <- ggplot(data = data, mapping = aes(x = UMAP_1,
                                        y = UMAP_2))   + 
  guides(alpha="none") +
  geom_point(colour = "#737373", size = 0.05) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='NA') +
  scale_fill_gradientn(limits=c(0,0.020),breaks=c(0,0.005,0.010,0.015,0.020),colours=c("navyblue",  "white","firebrick"))+theme_classic()
ggsave(filename = "Control_desity_UMAP.pdf", plot = p1, device="pdf", width=5.5, height=4)

CDAHFD<- integration_new[,integration_new@meta.data$orig.ident %in% c("CDAHFD_1","CDAHFD_2")]
data <- cbind(Embeddings(object=CDAHFD[['umap']]),FetchData(CDAHFD,'treatment'))
p2 <- ggplot(data = data, mapping = aes(x = UMAP_1,
                                        y = UMAP_2))   + 
  guides(alpha="none") +
  geom_point(colour = "#2171B5", size = 0.05) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='NA') +
  scale_fill_gradientn(limits=c(0,0.020),breaks=c(0,0.005,0.010,0.015,0.020),colours=c("navyblue",  "white","firebrick"))+theme_classic()
ggsave(filename = "CDAHFD_desity_UMAP.pdf", plot = p1, device="pdf", width=5.5, height=4)

BDL<- integration_new[,integration_new@meta.data$orig.ident %in% c("BDL_1","BDL_2")]
data <- cbind(Embeddings(object=BDL[['umap']]),FetchData(BDL,'treatment'))
p3 <- ggplot(data = data, mapping = aes(x = UMAP_1,
                                        y = UMAP_2))   + 
  guides(alpha="none") +
  geom_point(colour = "#238B45", size = 0.05) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='NA') +
  scale_fill_gradientn(limits=c(0,0.020),breaks=c(0,0.005,0.010,0.015,0.020),colours=c("navyblue",  "white","firebrick"))+theme_classic()
ggsave(filename = "BDL_desity_UMAP.pdf", plot = p1, device="pdf", width=5.5, height=4)

bind<-ggarrange(p1,p2,p3,ncol = 3,nrow = 1,common.legend = TRUE,legend = "bottom")
ggsave(filename = "Bind_desity_UMAP.pdf", plot = bind, device="pdf", width=13, height=4.5)

#Proportions of each cell cluster in each treatment condition
table(integration_new$treatment)
My_levels <- c('Control','CDAHFD','BDL')
integration_new$treatment<-factor(integration_new$treatment, levels= My_levels)

prop.table(table(Idents(integration_new)))
table(Idents(integration_new),integration_new$treatment)
Cellratio<-prop.table(table(Idents(integration_new),integration_new$treatment),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
colourCount=length(unique(Cellratio$Var1))
cells_in_treatment<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  coord_flip()+theme_classic()+labs(x="Sample",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values=colours)

cells_in_treatment<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  theme_classic()+labs(x="Sample",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values=colours)

ggsave(filename = "Total_cluster_Freq.pdf", plot = cells_in_treatment, device="pdf", width=10, height=6)
ggsave(filename = "Total_cluster_Freq_vertical.pdf", plot = cells_in_treatment, device="pdf", width=6, height=6)

#percentage of cell cluster from control, NASH and BDL cells out of total portal niche cells.
table(integration_new$orig.ident)
prop.table(table(Idents(integration_new)))
table(Idents(integration_new),integration_new$orig.ident)
Cellratio<-prop.table(table(Idents(integration_new),integration_new$orig.ident),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
write.table(Cellratio, 
            file="Cellratio.txt", sep="\t",row.names=FALSE, quote=FALSE)

CellRa<-read.table("Cellratio.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
CellRa$group <-Cellratio$Var2
CellRa$group<-gsub('.{2}$', '', CellRa$group)
CellRa$group <- factor(CellRa$group,levels=c("Control","CDAHFD","BDL"))
#Cholangiocytes
p1<-ggplot(data=CellRa[CellRa$Var1=="Cholangiocytes",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,50)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Cholangiocytes") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Cholangiocytes_Freq.pdf", plot = p1, device="pdf", width=3.5, height=2.5)
# LSECs
p2<-ggplot(data=CellRa[CellRa$Var1=="LSECs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,40)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="LSECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "LSECs_Freq.pdf", plot = p2, device="pdf", width=3.5, height=2.5)
#Monocytes & Macrophages
p3<-ggplot(data=CellRa[CellRa$Var1=="Monocytes & Macrophages",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,40)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Monocytes & Macrophages") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Monocytes & Macrophages_Freq.pdf", plot = p3, device="pdf", width=3.5, height=2.5)
##Fibroblasts & Myofibroblasts
p4<-ggplot(data=CellRa[CellRa$Var1=="Fibroblasts & Myofibroblasts",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,30)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Fibroblasts & Myofibroblasts") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Fibroblasts & Myofibroblasts_Freq.pdf", plot = p4, device="pdf", width=3.5, height=2.5)

#cDCs
p5<-ggplot(data=CellRa[CellRa$Var1=="cDCs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,10)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="cDCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "cDCs_Freq.pdf", plot = p5, device="pdf", width=3.5, height=2.5)
#Portal vein ECs
p6<-ggplot(data=CellRa[CellRa$Var1=="Portal vein ECs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,10)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Portal vein ECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Portal vein ECs_Freq.pdf", plot = p6, device="pdf", width=3.5, height=2.5)
#T & NK cells
p7<-ggplot(data=CellRa[CellRa$Var1=="T & NK cells",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,8)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="T & NK cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "T cells_Freq.pdf", plot = p7, device="pdf", width=3.5, height=2.5)

#pDCs
p8<-ggplot(data=CellRa[CellRa$Var1=="pDCs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="pDCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "pDCs_Freq.pdf", plot = p8, device="pdf", width=3.5, height=2.5)
#Lymphatic ECs
p9<-ggplot(data=CellRa[CellRa$Var1=="Lymphatic ECs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,2)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Lymphatic ECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Lymphatic ECs_Freq.pdf", plot = p9, device="pdf", width=3.5, height=2.5)
#VSMCs
p10<-ggplot(data=CellRa[CellRa$Var1=="VSMCs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="VSMCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "VSMCs_Freq.pdf", plot = p10, device="pdf", width=3.5, height=2.5)
#KCs
p11<-ggplot(data=CellRa[CellRa$Var1=="KCs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,7)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="KCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "KCs_Freq.pdf", plot = p11, device="pdf", width=3.5, height=2.5)

#B cells
p12<-ggplot(data=CellRa[CellRa$Var1=="B cells",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,4)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="B cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "B cells_Freq.pdf", plot = p12, device="pdf", width=3.5, height=2.5)

#HSCs
p13<-ggplot(data=CellRa[CellRa$Var1=="HSCs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,3)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="HSCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "HSCs_Freq.pdf", plot = p13, device="pdf", width=3.5, height=2.5)

#Neutrophils
p14<-ggplot(data=CellRa[CellRa$Var1=="Neutrophils",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,4)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Neutrophils") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Neutrophils_Freq.pdf", plot = p14, device="pdf", width=3.5, height=2.5)
#Arterial ECs
p15<-ggplot(data=CellRa[CellRa$Var1=="Arterial ECs",], aes(group,Freq*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,2)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Arterial ECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Arterial ECs_Freq.pdf", plot = p15, device="pdf", width=3.5, height=2.5)

pbind<-p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+ 
  plot_layout(ncol = 5)
ggsave(filename = "combind_Freq.pdf", plot = pbind, device="pdf", width=20, height=9)

#Heatmap

options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)

obj.markers <- FindAllMarkers(integration_new, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
write.table(obj.markers, file="obj.markers.txt", sep="\t", row.names=FALSE, quote=FALSE)

plan("sequential")

top10 = obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, file="Top10DEGs.txt", sep="\t", row.names=FALSE, quote=FALSE)

Top10DEGs<-read.table("Top10DEGs.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
colours=c("#3E8104", "#1F78B4", "#B2DF8A", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#DFC27D", "#E31A1C", "#FDBF6F",
          "#FF7F00", "#CAB2D6", "#654f75", "darkgrey", "#B15928", "#FB8072", "#80B1D3", 
          "darkgreen", "#F4A582", "#01665E")

cols.use <- list(cell_type=colours,
                 treatment=c("#737373","#2171B5","#238B45"))
Heatmap<-DoMultiBarHeatmap(subset(integration_new,downsample=50), features = Top10DEGs$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#or
Heatmap<-DoMultiBarHeatmap(integration_new,features = Top10DEGs$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename = "integration_new_heatmap_addanno.pdf", plot = Heatmap, device="pdf", width=15, height=15)

#vlnplot
features<- c("Epcam",  "C1qa","Ms4a7","Clec4f","Plbd1","Siglech","S100a8",
             "Trbc2","Cd79a","Dpt","Angptl6","Myh11","Upk3b","Dnase1l3", "Adgrg6","Lyve1","Sema3g","Apoa2")

colours2=c("#3E8104", "#1F78B4","#1F78B4",  "#B2DF8A", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#DFC27D", "#E31A1C", "#FDBF6F",
           "#FF7F00", "#CAB2D6",  "#654f75", "darkgrey", "#B15928", "#FB8072", "#80B1D3", 
           "darkgreen", "#F4A582", "#01665E")
VlnP<-VlnPlot(integration_new, features =  features, pt.size =0, cols=colours2,
              stack = T,flip = T )+labs(x="",y="Gene expression")+ theme(axis.title.y.right = element_text(face="italic",color = "#000000"))+
  theme(axis.text.y =element_text(color = "#000000",size=10))+theme(axis.text.x = element_text(angle = 90))+NoLegend()

ggsave(filename = "Vlnplot.pdf", plot = VlnP, device="pdf", width=10, height=10)

##Batch effect & QC
integration_new$orig.ident<- factor(integration_new$orig.ident, levels = c("Control_1","Control_2","Control_3","CDAHFD_1","CDAHFD_2","BDL_1","BDL_2"))

BaUmap <- DimPlot(integration_new, reduction = "umap",group.by = "orig.ident", label = F, repel = TRUE,pt.size = 0.02)+
  scale_colour_manual( values = c(brewer.pal(11,"RdGy")[7:9],brewer.pal(11,"RdYlBu")[9:10],brewer.pal(11,"PRGn")[9:10]))
ggsave(filename = "Batch_effect_Umap.pdf", plot = BaUmap, device="pdf", width=9.3, height=7.2)


QC<-VlnPlot(integration_new, group.by = "orig.ident",pt.size =0.05,features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),cols = c(brewer.pal(11,"PRGn")[9:10],brewer.pal(11,"RdYlBu")[9:10],brewer.pal(11,"RdGy")[7:9]), ncol = 3)&
  theme(axis.text.x = element_text(angle = 90))&labs(x="")
ggsave(filename = "QC.pdf", plot = QC, device="pdf", width=9.3, height=4)

exprs <- data.frame(FetchData(object = integration_new, vars=c("nFeature_RNA", "nCount_RNA", "percent.mito")))
#correlation among cell clusters
integration_new$cell_type<-integration_new@active.ident

table(integration_new$cell_type)
av <-AverageExpression(integration_new,
                       group.by = "cell_type",
                       assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av, 1, sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = 'spearman'))
cor<-pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(15))#,color=colorRampPalette(c("navy","white","firebrick3"))(50)
ggsave(filename = "correlation_heatmap.pdf", plot = cor, device="pdf", width=8.5, height=8)

