---
title: "Figure 2_Figure S2_Figure S3"
date: '2023-07-10'
author: "Lin LEI"
---
setwd("~/2023-03-12-Figure2")
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(DT)
library(tidyverse)
library(patchwork)
library(harmony)
set.seed(123)
#load the integrated datasets object
integration <- readRDS("~/Figure1_new/integration_new.rds")
#subseting Myeloid populations
Myeloid_cell<-subset(integration, idents =c(
  "Monocytes & Macrophages", 
  "cDCs", 
  "KCs", 
  "pDCs",
  "Neutrophils"))
p2 <- DimPlot(Myeloid_cell, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)

a<-ls()
rm(list = a[which(a!='integration')])

scMerge_QC<- NormalizeData(Myeloid_cell) %>% FindVariableFeatures() 
?CaseMatch
library(stringr)
s.genes.m<-str_to_title(cc.genes$s.genes)
g2m.genes.m<-str_to_title(cc.genes$g2m.genes)
CaseMatch(c(s.genes.m,g2m.genes.m),VariableFeatures(scMerge_QC))
g2m_genes=CaseMatch(search = g2m.genes.m,match = rownames(scMerge_QC))
s_genes=CaseMatch(search = s.genes.m,match = rownames(scMerge_QC))
scMerge_QC<-CellCycleScoring(object = scMerge_QC,g2m.features = g2m_genes,s.features = s_genes)

#re-clustering and de-batch effect
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
Myeloid_cell_QC<-ScaleData(scMerge_QC,vars.to.regress=c("S.Score","G2M.Score"),features = rownames(scMerge_QC))
##Explicitly close multisession workers, if they were used
plan("sequential")

Myeloid_cell_QC<- RunPCA(Myeloid_cell_QC,features = VariableFeatures(Myeloid_cell_QC),verbose=FALSE)

ElbowPlot(Myeloid_cell_QC)
pc.num=1:20
#Harmony re-integrating
Myeloid_cell_Har2<- Myeloid_cell_QC%>% RunHarmony("orig.ident", plot_convergence = TRUE, theta = 1)
rm(Myeloid_cell_QC)
Myeloid<- Myeloid_cell_Har2 %>% RunUMAP(reduction = "harmony", dims = 1:40) %>% FindNeighbors(reduction = "harmony", dims = 1:40) %>% FindClusters(resolution = seq(from=0.1, to=0.4, by=0.1)) %>% identity()
#check the data with cell-type specific markers
DimPlot(Myeloid, reduction = "umap",label = T, group.by = 'RNA_snn_res.0.3',pt.size = 0.05)
#DC cells
P1<-FeaturePlot(object = Myeloid, features = c("Runx2","Ccr9","Siglech","Spib","Irf8","Tcf4","Xcr1","Clec9a","Ccr7","H2-Aa"),reduction = "umap",order = F, blend = F, pt.size = 0.05, combine = T,ncol = 4)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
P2<-FeaturePlot(object = Myeloid, features = c("Itgax","Cd209a","Ffar2",
                                               "Klrd1","Itgam","Mgl2","Clec10a","Cd7"),reduction = "umap",order = F, blend = F, pt.size = 0.05, combine = T,ncol = 4)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

P3<-FeaturePlot(object =Myeloid, features = c("Spp1","Trem2","Gpnmb","Ptprc","Mki67","Cavin2","Cxcl9"),
                reduction = "umap",order = F, blend = F, pt.size = 0.05, combine = T,ncol = 3)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
P4<-FeaturePlot(object =Myeloid, features = c("Clec4b1","Timd4","Slc16a9","Cd5l","Slc40a1","Clec4f","Adgre1","Itgam","Ccr2"),reduction = "umap",order = F, blend = F, pt.size = 0.05, combine = T,ncol = 4)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
P5<-FeaturePlot(object = Myeloid, features = c("Pglyrp1","Spn","Trem3"),reduction = "umap",order = T, blend = F, pt.size = 0.05, combine = T,ncol = 4)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#Monocyte cells
P6<-FeaturePlot(object = Myeloid, features = c("Chil3","Ly6c2","S100a4","Sell","Cd177","F13a1","Gm9733"),reduction = "umap",order = T, blend = F, pt.size = 0.05, combine = T,ncol = 4)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#rename the subclusters
Myeloid<-SetIdent(Myeloid,value=Myeloid@meta.data$RNA_snn_res.0.3)
new.cluster.ids <- c("RM1s", 
                     "cCD2s", 
                     "Monocytes", 
                     "RM2s", 
                     "cCD1s",
                     "YS-derived KCs",
                     "BM-derived KCs",
                     "Neutrophils",
                     "BD-LAMs",
                     "pDCs",
                     "cCD1s",
                     "Monocytes", 
                     "Migratory cDCs")
names(new.cluster.ids) <- levels(Myeloid)
Myeloid<- RenameIdents(Myeloid, new.cluster.ids)
levels(Myeloid)
#define the color of cell subclusters
colours2=c("#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388","#4BB394","#46806F","#D95F02","#8FF7D8","#7570B3")
p1 <- DimPlot(Myeloid, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)+ scale_colour_manual(values = colours2)+NoLegend()
ggsave(filename = "Umap.pdf", plot = p1, device="pdf", width=8, height=7.2)

#reordering subclusters
Myeloid@active.ident <- factor(Myeloid@active.ident, 
                               levels=c("RM1s", 
                                        "RM2s", 
                                        "Monocytes", 
                                        "BD-LAMs",
                                        "YS-derived KCs",
                                        "BM-derived KCs",
                                        "cCD1s",
                                        "cCD2s", 
                                        "pDCs",
                                        "Migratory cDCs",
                                        "Neutrophils"
                               ))
levels(Myeloid)
#re-name the subclusters
new.cluster.ids <- c("RM1s", 
                     "RM2s", 
                     "Monocytes", 
                     "BD-LAMs",
                     "YS-derived KCs",
                     "BM-derived KCs",
                     "cDC1s",
                     "cDC2s", 
                     "pDCs",
                     "Migratory cDCs",
                     "Neutrophils"
)
names(new.cluster.ids) <- levels(Myeloid)
Myeloid <- RenameIdents(Myeloid, new.cluster.ids)
levels(Myeloid)
#total Umap
ToUmap <-DimPlot(Myeloid, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1)+ scale_colour_manual( values = colours2)+NoLegend()
ggsave(filename = "Mye_ToUmap.pdf", plot = ToUmap, device="pdf", width=5.5, height=5)

#splitted UMAP
ToUmapsplit <- DimPlot(Myeloid, reduction = "umap", label = TRUE,split.by = "treatment", repel = TRUE,pt.size = 1)+ scale_colour_manual( values = colours2)+NoLegend()
ggsave(filename = "ToUmapsplit.pdf", plot = ToUmapsplit, device="pdf", width=12, height=4.5)

#Featureplot
p1=FeaturePlot(object = Myeloid, features = "Cx3cr1",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p1=p1+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p2=FeaturePlot(object = Myeloid, features = "Ccl8",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p2=p2+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p3=FeaturePlot(object = Myeloid, features = "Ly6c2",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p3=p3+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p4=FeaturePlot(object = Myeloid, features = "Gpnmb",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p4=p4+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p5=FeaturePlot(object = Myeloid, features = "Clec4f",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p5=p5+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p6=FeaturePlot(object = Myeloid, features = "Xcr1",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p6=p6+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p7<-ggarrange(p1,p2,p3,p4,p5,p6,ncol = 3,nrow = 2,common.legend = T, legend = "right")+theme(plot.title = element_text(face="italic",color = 'black'))
ggsave(filename = "totalmaker_UMAP.pdf", plot = p7, device="pdf", width=9, height=6)

ggsave(filename = "totalmaker_UMAP.png", plot = p7, device="png", width=9, height=6)

#Proportions of each cell cluster in each treatment condition
table(Myeloid$treatment)
My_levels <- c('Control','CDAHFD','BDL')
Myeloid$treatment<-factor(Myeloid$treatment, levels= My_levels)

prop.table(table(Idents(Myeloid)))
table(Idents(Myeloid),Myeloid$treatment)
Cellratio<-prop.table(table(Idents(Myeloid),Myeloid$treatment),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1,levels=c("RM1s", 
                                                 "RM2s", 
                                                 "Monocytes", 
                                                 "BD-LAMs",
                                                 "YS-derived KCs",
                                                 "BM-derived KCs",
                                                 "cDC1s",
                                                 "cDC2s", 
                                                 "pDCs",
                                                 "Migratory cDCs",
                                                 "Neutrophils"
))
table(Cellratio$Var1)
colourCount=length(unique(Cellratio$Var1))
cells_in_treatment<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  theme_classic()+labs(x="Sample",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values = c("#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388","#4BB394","#46806F","#D95F02","#8FF7D8","#7570B3"))
ggsave(filename = "Myeloid_cluster_Freq_vertical0801.pdf", plot = cells_in_treatment, device="pdf", width=6, height=6)

#percentage of cell cluster from control, NASH and BDL cells out of total portal niche cells.
#table(integration$orig.ident)
#Control_1   Control_2   Control_3 CDAHFD_1 CDAHFD_2    BDL_1    BDL_2 
#3493     3990     1630     8273     7786     9143     9704 

for (i in 1:nrow(b) ) {
  if(b$Var2[i]=="Control_1"){b$Freq2[i]=b$Freq[i]/3493
  }
  else if (b$Var2[i]=="Control_2"){b$Freq2[i]=b$Freq[i]/3990
  }
  else if (b$Var2[i]=="Control_3"){b$Freq2[i]=b$Freq[i]/1630
  }
  else if (b$Var2[i]=="CDAHFD_1"){b$Freq2[i]=b$Freq[i]/8273
  }
  else if (b$Var2[i]=="CDAHFD_2"){b$Freq2[i]=b$Freq[i]/7786
  }
  else if (b$Var2[i]=="BDL_1"){b$Freq2[i]=b$Freq[i]/9143
  }
  else{b$Freq2[i]=b$Freq[i]/9704
  }
}
options(scipen=200)
Cellratio<-as.data.frame(b)

Cellratio$group <- Cellratio$Var2
Cellratio$group<-gsub('.{2}$', '', Cellratio$group)
write.table(Cellratio, 
            file="Cellratio_Tot.txt", sep="\t",row.names=FALSE, quote=FALSE)

CellRa<-read.table("Cellratio_Tot.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)

CellRa$group <- factor(CellRa$group,levels=c("Control","CDAHFD","BDL"))
##RM1s
p0<-ggplot(data=CellRa[CellRa$Var1=="RM1s",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,20)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="RM1s") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "RM1s_Freq_Total.pdf", plot = p0, device="pdf", width=3.5, height=2.5)

##RM2s
p1<-ggplot(data=CellRa[CellRa$Var1=="RM2s",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,10)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="RM2s") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "RM2s_Freq_Total.pdf", plot = p1, device="pdf", width=3.5, height=2.5)
##Monocytes
p2<-ggplot(data=CellRa[CellRa$Var1=="Monocytes",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,10)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Monocytes") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Monocytes_Freq_Total.pdf", plot = p2, device="pdf", width=3.5, height=2.5)
##BD-LAMs
p3<-ggplot(data=CellRa[CellRa$Var1=="BD-LAMs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,3)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="BD-LAMs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "BD-LAMs_Freq_Total.pdf", plot = p3, device="pdf", width=3.5, height=2.5)
##YS-derived KCs
p4<-ggplot(data=CellRa[CellRa$Var1=="YS-derived KCs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,10)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="YS-derived KCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "YS-derived KCs_Freq_Total.pdf", plot = p4, device="pdf", width=3.5, height=2.5)
##BM-derived KCs
p5<-ggplot(data=CellRa[CellRa$Var1=="BM-derived KCs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="BM-derived KCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "BM-derived KCs_Freq_Total.pdf", plot = p5, device="pdf", width=3.5, height=2.5)
##cCD1s 
p6<-ggplot(data=CellRa[CellRa$Var1=="cDC1s",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,4)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="cDC1s") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "cDC1s_Freq_Total.pdf", plot = p6, device="pdf", width=3.5, height=2.5)
##cCD2s 
p7<-ggplot(data=CellRa[CellRa$Var1=="cDC2s",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,8)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="cDC2s") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "cDC2s_Freq_Total.pdf", plot = p7, device="pdf", width=3.5, height=2.5)
##pDCs
p8<-ggplot(data=CellRa[CellRa$Var1=="pDCs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,2)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="pDCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "pDCs_Freq_Total.pdf", plot = p8, device="pdf", width=3.5, height=2.5)
##Migratory cDCs 
p9<-ggplot(data=CellRa[CellRa$Var1=="Migratory cDCs",], aes(group,Freq2*100+0.01))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Migratory cDCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Migratory cDCs_Freq_Total.pdf", plot = p9, device="pdf", width=3.5, height=2.5)
library(ggpubr)
plotDCs<-p6+p7+p9+ 
  plot_layout(ncol = 3)
ggsave(filename = "plotDCs.pdf", plot = plotDCs, device="pdf", width=12, height=3)

plotDCs<-p6+p7+p9+ 
  plot_layout(ncol = 1)
ggsave(filename = "plotDCs2.pdf", plot = plotDCs, device="pdf", width=4, height=9)
##Neutrophilss 
p10<-ggplot(data=CellRa[CellRa$Var1=="Neutrophils",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.01),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Neutrophils") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
ggsave(filename = "Neutrophils_Freq_Total.pdf", plot = p10, device="pdf", width=3.5, height=2.5)
plotDCs<-ggarrange(p3,p0,p1,p2,p10,ncol = 5,nrow = 1,common.legend = T, legend = "right")+theme(plot.title = element_text(face="italic",color = 'black'))
ggsave(filename = "plotDCs.pdf", plot = plotDCs, device="pdf", width=12, height=3)
pbind<-p3+p0+p1+p2+p10+ 
  plot_layout(ncol = 5)
ggsave(filename = "plot_RM_mono.pdf", plot = pbind, device="pdf", width=20, height=3)

#clean the environment
a<-ls()
rm(list = a[which(a!='integration'&a!="Myeloid")])
gc()

#Cluster DEGs
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
Myeloid_markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
write.table(Myeloid_markers, file="Myeloid_markers.txt", sep="\t", row.names=FALSE, quote=FALSE)
Myeloid_markers  <-read.table(file="Myeloid_markers.txt",header=T,sep="\t",quote="",check.names=F)

plan("sequential")

top10 = Myeloid_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, file="Top10DEGs.txt", sep="\t", row.names=FALSE, quote=FALSE)

top10 <-read.table(file="Top10DEGs.txt",header=T,sep="\t",quote="",check.names=F)
#Heatmap
colours=c("#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388","#4BB394","#46806F","#D95F02","#8FF7D8","#7570B3")
cols.use <- list(cell_type=colours,
                 treatment=c("#737373","#2171B5","#238B45"))
Myeloid@meta.data$treatment<- factor(x = Myeloid@meta.data$treatment, levels=c('Control', 'CDAHFD', 'BDL'))
table(Myeloid$treatment)
My_levels <- c('Control','CDAHFD','BDL')
Myeloid$treatment<-factor(Myeloid$treatment, levels= My_levels)
Heatmap<-DoMultiBarHeatmap(Myeloid, features = top10$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename = "Myeloid_heatmap_addanno.pdf", plot = Heatmap, device="pdf", width=10, height=15)
ggsave(filename = "Myeloid_heatmap_addanno2.pdf", plot = Heatmap, device="pdf", width=8, height=8)

#vlnplot
features<- c("Ms4a7", "C1qa","Ccl8",
             "Chil3","Thbs1", "Gpnmb",
             "Vsig4","Clec4f","Timd4",
             "Xcr1","Cd209a","Siglech","Ccr7","Csf3r")
Myeloid$cell_type<-Myeloid@active.ident
VlnP<-VlnPlot(Myeloid, features =  features, pt.size =0, cols=colours2,fill.by = "ident",
              stack = T,flip = T )+labs(x="",y="Gene expression")+ theme(axis.text.y.right = element_text(face="italic",size=10,color = "#000000"))+
  theme(axis.text.y =element_text(color = "#000000",size=10))+theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+NoLegend()
ggsave(filename = "Mye_Vlnplot.pdf", plot = VlnP, device="pdf", width=4, height=5.5)

##Coefficient of variation analysis，adapted from https://www.nature.com/articles/s43587-022-00246-4#code-availability

#```{r COV, message=FALSE, warning=FALSE}
# CV = SD/Mean. The CV will give you the extent of variability in your gene expression dataset
library(reshape2)
#makes a matrix with cell names as colnames and gene as rownames
dsmyeloid<-subset(Myeloid,downsample=100)
ctx<-GetAssayData(dsmyeloid, slot = "data", assay = "RNA")
#ctx<-GetAssayData(Myeloid, slot = "data", assay = "RNA")
ctx<-t(ctx)
dim(ctx)
ctx<-as_tibble(ctx, rownames = "cell_name")
ctx
id<-FetchData(dsmyeloid, vars = c("cell_type", "treatment"), slot = "data")
id<-as_tibble(id, rownames= "cell_name")
#id
id<-mutate(id, group= paste(cell_type, treatment, sep = "&"))
id
ctx<-full_join(ctx, id, by = "cell_name")
#ctx<-dplyr::select(ctx, cell_name, -ident, - stim)
#gets coefficient of variation for each gene per group
covfun<-function(x) {sd(x)/mean(x)}
covres<-group_by(ctx, group) %>%
  summarise(across(where(is.numeric), ~ covfun(.x)))
covres<-as.data.frame(covres)
meltcovres <- melt(covres,id.vars="group")
meltcovres <-meltcovres  %>%
  separate(group, c("cell_type", "treatment"), "&")
# grouped boxplot
meltcovres2<-meltcovres %>% drop_na(value)
# setting the subcluster order
meltcovres2$treatment<-factor(meltcovres$treatment, levels = c("Control", "CDAHFD","BDL"))
My_levels <- c("RM1s", 
               "RM2s", 
               "Monocytes", 
               "BD-LAMs",
               "YS-derived KCs",
               "BM-derived KCs",
               "cDC1s",
               "cDC2s", 
               "pDCs",
               "Migratory cDCs",
               "Neutrophils"
)
meltcovres2$cell_type<-factor(meltcovres2$cell_type, levels= My_levels)
meltcovres


c.v.plot<-ggplot(meltcovres2, aes(x=cell_type, y=value, fill=treatment)) + 
  geom_boxplot() + theme_bw() + scale_fill_manual(values=c("Control" = "#737373", "CDAHFD" = "#2171B5","BDL" = "#238B45"))+
  theme(panel.grid = element_blank(),  axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+theme(axis.text.x = element_text(angle = 90))


ggsave(filename = "c.v.plot_myeloid.pdf", plot = c.v.plot, device="pdf", width=8, height=3.8)

#Check the significance of C.v value
library(rstatix)
wilcox_res<-meltcovres2 %>%
  group_by(cell_type) %>%
  wilcox_test(value ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="c.v._wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)

#Gene signature calculation

library(dplyr)
library(Seurat)
library(patchwork)
library(macSpectrum)

#load the gene list of gene modules,M1 and M2 polarization analysis
geneset_list <- readxl::read_xlsx("~Myeloid_gene_signature/geneset_list.xlsx")
View(geneset_list)
M1 <- as.character(geneset_list$M1)
M1 <- M1[1:15]
M2 <- as.character(geneset_list$M2)
M2 <- M2[1:31]

#calculate the M1 and M2 gene module score
Myeloid<-AddModuleScore(object = Myeloid,features=list(M1),ctrl = 100,name='M1signature')
Myeloid<-AddModuleScore(object = Myeloid,features=list(M2),ctrl = 100,name='M2signature')

P2<-FeaturePlot(object =Myeloid_mono, features = "M1signature1",reduction = "umap",order = F, blend = F, pt.size = 0.05, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
P3<-FeaturePlot(object =Myeloid_mono,features = "M2signature1", reduction = "umap",order = F,blend = F, pt.size = 0.05, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
P4<-ggarrange(P2,P3,common.legend = TRUE, legend = "right")
ggsave(filename = "M1_M2_signature_Umap_noncutoff2.pdf", plot = P4, device="pdf", width=8.8, height=4)

#load the gene list of gene modules,Profibrogenesis
Profibrogenesis<-as.character(geneset_list$profibrogenesis)
Profibrogenesis<-Profibrogenesis[1:19]
#calculate Profibrogenesis gene module score
Myeloid<-AddModuleScore(object = Myeloid,features=list(Profibrogenesis),ctrl = 100,name='Profibrogenesis')

order<-c("RM1s", 
         "RM2s", 
         "Monocytes", 
         "BD-LAMs",
         "YS-derived KCs",
         "BM-derived KCs",
         "cCD1s",
         "cCD2s", 
         "pDCs",
         "Migratory cDCs",
         "Neutrophils")
boxplot<-ggboxplot(Myeloid@meta.data, x = "cell_type", y = "Profibrogenesis1",order = order,
                   color = "black",fill = 'cell_type',palette=c("#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388","#4BB394","#46806F","#D95F02","#8FF7D8","#7570B3"),bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
ggsave(filename = "Profibrogenesis_core_boxplot_celltype.pdf", plot = boxplot, device="pdf", width=8, height=5)

#MacSpectrum analysis--https://macspectrum.uconn.edu/user_manual/##
My_levels <- c('Control','CDAHFD','BDL')
Myeloid$treatment<-factor(Myeloid$treatment, levels= My_levels)
Myeloid@active.ident <- factor(Myeloid@active.ident, 
                               levels=c("RM1s", 
                                        "RM2s", 
                                        "Monocytes", 
                                        "BD-LAMs",
                                        "YS-derived KCs",
                                        "BM-derived KCs",
                                        "cCD1s",
                                        "cCD2s", 
                                        "pDCs",
                                        "Migratory cDCs",
                                        "Neutrophils"
                               ))
Myeloid$celltype.treatment<- paste(Myeloid@active.ident, Myeloid$treatment, sep = "_")
Myeloid$celltype<- Idents(Myeloid)
My_levels <- c('RM1s_Control',
               'RM1s_CDAHFD',
               'RM1s_BDL',
               'RM2s_Control',
               'RM2s_CDAHFD',
               'RM2s_BDL',
               'Monocytes_Control',
               'Monocytes_CDAHFD',
               'Monocytes_BDL',
               'BD-LAMs_Control',
               'BD-LAMs_CDAHFD',
               'BD-LAMs_BDL',
               'YS-derived KCs_Control',
               'YS-derived KCs_CDAHFD',
               'YS-derived KCs_BDL',
               'BM-derived KCs_Control',
               'BM-derived KCs_CDAHFD',
               'BM-derived KCs_BDL')
Myeloid$celltype.treatment<-factor(Myeloid$celltype.treatment, levels= My_levels)

##RM1
RM1<-subset(x = Myeloid, idents = c("RM1s"))
levels(RM1)
expr_matrix <- GetAssayData(object = RM1, slot = "counts")
dim(expr_matrix)
expr_matrix[1:4,1:4]
class(expr_matrix)
write.csv(expr_matrix,file = 'RM1_expr_matrix.csv')
Matr<-read.csv("RM1_expr_matrix.csv",  
               header = T, sep= ",") 
#convert gene ID
library(clusterProfiler)
ensembl<- bitr(Matr$X,fromType = 'SYMBOL',
               toType = c('ENSEMBL'),
               OrgDb='org.Mm.eg.db')

colnames(Matr)[1] <- 'SYMBOL'
mac_mtx<-dplyr::inner_join(x=ensembl,y=Matr,by="SYMBOL")
mac_mtx<-mac_mtx[,-1]
colnames(mac_mtx)[1] <- 'Ensembl_ID'

p_data <- RM1@meta.data
feature<-p_data$celltype.treatment

colnames(mac_mtx)[1] <- 'geneid'
macoutput<-macspec(mac_mtx, feature, select_hu_mo = "mou")

MPI<-ggplot(macoutput, aes(x=MPI,  fill=Feature))+ 
  geom_density(alpha=0.55,bw=1,colour="black",size=0.25)+theme_classic()+
  scale_fill_manual(values= c("#737373","#2171B5","#238B45"))+
  theme(
    text=element_text(size=15,color="black"),
    plot.title=element_text(size=15,family="myfont",face="bold.italic",hjust=.5,color="black"),
    legend.position=c(0.6,0.8),
    legend.background = element_blank()
  )+  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

ggsave(filename = "RM1_MPI.pdf", plot = MPI, device="pdf", width=7, height=4)

##RM2
RM2<-subset(x = Myeloid, idents = c("RM2s"))
levels(RM2)

expr_matrix <- GetAssayData(object = RM2, slot = "counts")
dim(expr_matrix)
expr_matrix[1:4,1:4]
class(expr_matrix)
write.csv(expr_matrix,file = 'RM2_expr_matrix.csv')
Matr<-read.csv("RM2_expr_matrix.csv",  
               header = T, sep= ",") 
ensembl<- bitr(Matr$X,fromType = 'SYMBOL',
               toType = c('ENSEMBL'),
               OrgDb='org.Mm.eg.db',
)
colnames(Matr)[1] <- 'SYMBOL'
mac_mtx<-dplyr::inner_join(x=ensembl,y=Matr,by="SYMBOL")
mac_mtx<-mac_mtx[,-1]
colnames(mac_mtx)[1] <- 'Ensembl_ID'

p_data <- RM2@meta.data
feature<-p_data$celltype.treatment
colnames(mac_mtx)[1] <- 'geneid'
macoutput<-macspec(mac_mtx, feature, select_hu_mo = "mou")
MPI<-ggplot(macoutput, aes(x=MPI,  fill=Feature))+ 
  geom_density(alpha=0.55,bw=1,colour="black",size=0.25)+theme_classic()+
  scale_fill_manual(values= c("#737373","#2171B5","#238B45"))+
  theme(
    text=element_text(size=15,color="black"),
    plot.title=element_text(size=15,family="myfont",face="bold.italic",hjust=.5,color="black"),
    legend.position=c(0.6,0.8),
    legend.background = element_blank()
  )+  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

ggsave(filename = "RM2_MPI.pdf", plot = MPI, device="pdf", width=7, height=4)

#BD-LAMs
BD_LAMs<-subset(x = Myeloid, idents = c("BD-LAMs"))
expr_matrix <- GetAssayData(object = BD_LAMs, slot = "counts")
dim(expr_matrix)
expr_matrix[1:4,1:4]
class(expr_matrix)
write.csv(expr_matrix,file = 'BD_LAMs_expr_matrix.csv')
Matr<-read.csv("BD_LAMs_expr_matrix.csv",  
               header = T, sep= ",") 
ensembl<- bitr(Matr$X,fromType = 'SYMBOL',
               toType = c('ENSEMBL'),
               OrgDb='org.Mm.eg.db',
)
colnames(Matr)[1] <- 'SYMBOL'
mac_mtx<-dplyr::inner_join(x=ensembl,y=Matr,by="SYMBOL")
mac_mtx<-mac_mtx[,-1]
colnames(mac_mtx)[1] <- 'Ensembl_ID'
write.csv(mac_mtx,file = 'BD_LAMs_mac_mtx.csv',row.names = F)
p_data <- BD_LAMs@meta.data
feature<-p_data$celltype.treatment
write.csv(feature,file = 'BD_LAMs_feature.csv')

colnames(mac_mtx)[1] <- 'geneid'
macoutput<-macspec(mac_mtx, feature, select_hu_mo = "mou")
library(ggplot2)
MPI<-ggplot(macoutput, aes(x=MPI,  fill=Feature))+ 
  geom_density(alpha=0.55,bw=1,colour="black",size=0.25)+theme_classic()+
  scale_fill_manual(values= c("#737373","#2171B5","#238B45"))+
  theme(
    text=element_text(size=15,color="black"),
    plot.title=element_text(size=15,family="myfont",face="bold.italic",hjust=.5,color="black"),
    legend.position=c(0.6,0.8),
    legend.background = element_blank()
  )+  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

ggsave(filename = "BD-LAMs_MPI.pdf", plot = MPI, device="pdf", width=7, height=4)
#Monocytes
Mono<-subset(x = Myeloid, idents = c("Monocytes"))
levels(Mono)
expr_matrix <- GetAssayData(object = Mono, slot = "counts")
dim(expr_matrix)
expr_matrix[1:4,1:4]
class(expr_matrix)
write.csv(expr_matrix,file = 'Mono_expr_matrix.csv')
Matr<-read.csv("Mono_expr_matrix.csv",  
               header = T, sep= ",") 
ensembl<- bitr(Matr$X,fromType = 'SYMBOL',
               toType = c('ENSEMBL'),
               OrgDb='org.Mm.eg.db',
)
colnames(Matr)[1] <- 'SYMBOL'
mac_mtx<-dplyr::inner_join(x=ensembl,y=Matr,by="SYMBOL")
mac_mtx<-mac_mtx[,-1]
colnames(mac_mtx)[1] <- 'Ensembl_ID'

p_data <- Mono@meta.data
feature<-p_data$celltype.treatment
colnames(mac_mtx)[1] <- 'geneid'
macoutput<-macspec(mac_mtx, feature, select_hu_mo = "mou")
MPI<-ggplot(macoutput, aes(x=MPI,  fill=Feature))+ 
  geom_density(alpha=0.55,bw=1,colour="black",size=0.25)+theme_classic()+
  scale_fill_manual(values= c("#737373","#2171B5","#238B45"))+
  theme(
    text=element_text(size=15,color="black"),
    plot.title=element_text(size=15,family="myfont",face="bold.italic",hjust=.5,color="black"),
    legend.position=c(0.6,0.8),
    legend.background = element_blank()
  )+  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

ggsave(filename = "Mono_MPI.pdf", plot = MPI, device="pdf", width=7, height=4)
#Ys-derived KC
YS<-subset(x = Myeloid, idents = c("YS-derived KCs"))
expr_matrix <- GetAssayData(object = YS, slot = "counts")
dim(expr_matrix)
expr_matrix[1:4,1:4]
class(expr_matrix)
write.csv(expr_matrix,file = 'YS_expr_matrix.csv')
Matr<-read.csv("YS_expr_matrix.csv",  
               header = T, sep= ",") 
ensembl<- bitr(Matr$X,fromType = 'SYMBOL',
               toType = c('ENSEMBL'),
               OrgDb='org.Mm.eg.db',
)
colnames(Matr)[1] <- 'SYMBOL'
mac_mtx<-dplyr::inner_join(x=ensembl,y=Matr,by="SYMBOL")
mac_mtx<-mac_mtx[,-1]
colnames(mac_mtx)[1] <- 'Ensembl_ID'
p_data <- YS@meta.data
feature<-p_data$celltype.treatment
colnames(mac_mtx)[1] <- 'geneid'
macoutput<-macspec(mac_mtx, feature, select_hu_mo = "mou")
MPI<-ggplot(macoutput, aes(x=MPI,  fill=Feature))+ 
  geom_density(alpha=0.55,bw=1,colour="black",size=0.25)+theme_classic()+
  scale_fill_manual(values= c("#737373","#2171B5","#238B45"))+
  theme(
    text=element_text(size=15,color="black"),
    plot.title=element_text(size=15,family="myfont",face="bold.italic",hjust=.5,color="black"),
    legend.position=c(0.6,0.8),
    legend.background = element_blank()
  )+  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

ggsave(filename = "YS_MPI.pdf", plot = MPI, device="pdf", width=7, height=4)
##BM-derived KC
BM<-subset(x = Myeloid, idents = c("BM-derived KCs"))
expr_matrix <- GetAssayData(object = BM, slot = "counts")
dim(expr_matrix)
expr_matrix[1:4,1:4]
class(expr_matrix)
write.csv(expr_matrix,file = 'BM_expr_matrix.csv')
Matr<-read.csv("BM_expr_matrix.csv",  
               header = T, sep= ",") 
ensembl<- bitr(Matr$X,fromType = 'SYMBOL',
               toType = c('ENSEMBL'),
               OrgDb='org.Mm.eg.db',
)
colnames(Matr)[1] <- 'SYMBOL'
mac_mtx<-dplyr::inner_join(x=ensembl,y=Matr,by="SYMBOL")
mac_mtx<-mac_mtx[,-1]
colnames(mac_mtx)[1] <- 'Ensembl_ID'
p_data <- BM@meta.data
feature<-p_data$celltype.treatment
colnames(mac_mtx)[1] <- 'geneid'
macoutput<-macspec(mac_mtx, feature, select_hu_mo = "mou")

MPI<-ggplot(macoutput, aes(x=MPI,  fill=Feature))+ 
  geom_density(alpha=0.55,bw=1,colour="black",size=0.25)+theme_classic()+
  scale_fill_manual(values= c("#737373","#2171B5","#238B45"))+
  theme(
    text=element_text(size=15,color="black"),
    plot.title=element_text(size=15,family="myfont",face="bold.italic",hjust=.5,color="black"),
    legend.position=c(0.6,0.8),
    legend.background = element_blank()
  )+  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

ggsave(filename = "BM_MPI.pdf", plot = MPI, device="pdf", width=7, height=4)

#heatmap of pro-inflammatory gene in Myeloid cells
Myeloid_mono<-subset(Myeloid, idents =c( "RM1s", "RM2s","Monocytes", "BD-LAMs","YS-derived KCs",  "BM-derived KCs"))
Myeloid_mono$celltype.treatment<- paste(Myeloid_mono@active.ident, Myeloid_mono$treatment, sep = "_")
Myeloid_mono<-SetIdent(Myeloid_mono,value=Myeloid_mono@meta.data$celltype.treatment)
My_levels <- c('RM1s_Control',
               'RM1s_CDAHFD',
               'RM1s_BDL',
               'RM2s_Control',
               'RM2s_CDAHFD',
               'RM2s_BDL',
               'Monocytes_Control',
               'Monocytes_CDAHFD',
               'Monocytes_BDL',
               'BD-LAMs_Control',
               'BD-LAMs_CDAHFD',
               'BD-LAMs_BDL',
               'YS-derived KCs_Control',
               'YS-derived KCs_CDAHFD',
               'YS-derived KCs_BDL',
               'BM-derived KCs_Control',
               'BM-derived KCs_CDAHFD',
               'BM-derived KCs_BDL')
Myeloid_mono$celltype.treatment<-factor(Myeloid_mono$celltype.treatment, levels= My_levels)
Myeloid_mono@active.ident <- factor(Myeloid_mono@active.ident, 
                                    levels=My_levels
)
levels(Myeloid_mono)
aveMyemono <- AverageExpression(Myeloid_mono, return.seurat = T)
features =  c("Cxcl2","Ccl2", "Il1b","Tnf")

suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = FALSE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE,
                               lines.width = NULL,
                               group.bar.height = 0.04, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())+ 
      theme(axis.text.y = element_text(size = 8,colour = "black"))
  }
  
  return(plot)
}
aveMyemono$group<-factor(aveMyemono@active.ident, levels= My_levels)
library(stringr)
library(grid)
aveMyemono$celltype<-str_remove(string = aveMyemono$group,pattern = c("_Control","_CDAHFD","_BDL"))

My_levels <- c("RM1s", "RM2s","Monocytes",                   
               "BD-LAMs","YS-derived KCs",  "BM-derived KCs")
aveMyemono$celltype<-factor(aveMyemono$celltype, levels= My_levels)
colours=c("#737373","#2171B5","#238B45","#737373","#2171B5","#238B45","#737373","#2171B5","#238B45",
          "#737373","#2171B5","#238B45","#737373","#2171B5","#238B45","#737373","#2171B5","#238B45")
cols.use <- list(group=colours,
                 celltype=c("#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388"))
Proinflam-Myeloid<-DoMultiBarHeatmap(aveMyemono, features = features,draw.lines = F, group.by = "celltype",additional.group.by = "group",additional.group.sort.by = c('group'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename = "Proinflam_heatmap.pdf", plot = Proinflam-Myeloid, device="pdf", width=7, height=2)

#Cx3cr1+ Ccl8+  cells density plot
library("Nebulosa")
library(ggpubr)
#Cx3cr1+ Ccl8+ macrophages density onto UMAP
p1 <- plot_density(chow, c("Cx3cr1", "Ccl8"),joint = TRUE,size = 0.4)&scale_color_distiller(limits=c(-0.0001,0.016),palette = "Spectral")
p2 <- plot_density(Nash, c("Cx3cr1", "Ccl8"),joint = TRUE,size = 0.4)&scale_color_distiller(limits=c(-0.0001,0.016),palette = "Spectral")
p3 <- plot_density(BDL, c("Cx3cr1", "Ccl8"),joint = TRUE,size = 0.4)&scale_color_distiller(limits=c(-0.0001,0.016),palette = "Spectral")
p4<-ggarrange(p1[[3]],p2[[3]],p3[[3]],ncol = 3,nrow = 1,common.legend = TRUE, legend = "bottom")
ggsave(filename = "Cx3cr1+ Ccl8+  cells density.pdf", plot = p4, device="pdf", width=10, height=3.8)

###Strip plot###
#DEGs analysis
Myeloid$celltype.treatment<- paste(Myeloid@active.ident, Myeloid$treatment, sep = "_")
Myeloid$celltype<- Idents(Myeloid)
Idents(Myeloid) <- "celltype.treatment"

cellfordeg<-levels(Myeloid$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(Myeloid,logfc.threshold = 0, ident.1 = paste0(cellfordeg[i],"_BDL"), ident.2 = paste0(cellfordeg[i],"_chow"), verbose = FALSE)
  write.csv(CELLDEG,paste0("BDL_T_",cellfordeg[i],".CSV"))
}

cellfordeg<-levels(Myeloid$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(Myeloid,logfc.threshold = 0, ident.1 = paste0(cellfordeg[i],"_CDAHFD"), ident.2 = paste0(cellfordeg[i],"_chow"), verbose = FALSE)
  write.csv(CELLDEG,paste0("CDAHFD_T_",cellfordeg[i],".CSV"))
}

list.files()

library(ggrepel)
library(RColorBrewer)

# Read the DEGs list_BDL vs Control
RM1 <- read.csv("BDL_T_RM1s.CSV", row.names = 1)
RM2 <- read.csv("BDL_T_RM2s.CSV", row.names = 1)
Mono <- read.csv("BDL_T_Monocytes.CSV", row.names = 1)
LAM <- read.csv("BDL_T_BD-LAMs.CSV", row.names = 1)
YS <- read.csv("BDL_T_YS-derived KCs.CSV", row.names = 1)
BM <- read.csv("BDL_T_BM-derived KCs.CSV", row.names = 1)
cDC1 <- read.csv("BDL_T_cCD1s.CSV", row.names = 1)
cDC2 <- read.csv("BDL_T_cCD2s.CSV", row.names = 1)
pDC<- read.csv("BDL_T_pDCs.CSV", row.names = 1)
Mig <- read.csv("BDL_T_Migratory cDCs.CSV", row.names = 1)
Neu <- read.csv("BDL_T_Neutrophils.CSV", row.names = 1)
# Read the DEGs list_CDAHFD vs Control
RM1 <- read.csv("CDAHFD_T_RM1s.CSV", row.names = 1)
RM2 <- read.csv("CDAHFD_T_RM2s.CSV", row.names = 1)
Mono <- read.csv("CDAHFD_T_Monocytes.CSV", row.names = 1)
LAM <- read.csv("CDAHFD_T_BD-LAMs.CSV", row.names = 1)
YS <- read.csv("CDAHFD_T_YS-derived KCs.CSV", row.names = 1)
BM <- read.csv("CDAHFD_T_BM-derived KCs.CSV", row.names = 1)
cDC1 <- read.csv("CDAHFD_T_cCD1s.CSV", row.names = 1)
cDC2 <- read.csv("CDAHFD_T_cCD2s.CSV", row.names = 1)
pDC<- read.csv("CDAHFD_T_pDCs.CSV", row.names = 1)
Mig <- read.csv("CDAHFD_T_Migratory cDCs.CSV", row.names = 1)
Neu <- read.csv("CDAHFD_T_Neutrophils.CSV", row.names = 1)

padj <- 0.05
Typelist<-list(RM1,RM2,Mono,LAM, YS,BM,cDC1,cDC2,pDC,Mig,Neu)
Typelist<-list(RM1,RM2,Mono,LAM, YS,BM,cDC1,cDC2,pDC,Mig,Neu)
for(i in 1:length(Typelist)){
  Typelist[[i]]$change <- ifelse(Typelist[[i]]$avg_log2FC >= 0.25 & Typelist[[i]]$p_val_adj < padj, 
                                 "Up regulate",
                                 ifelse(Typelist[[i]]$avg_log2FC <  -0.25 & Typelist[[i]]$p_val_adj < padj,
                                        "Down regulate", 
                                        "Not significant"))
}

names(Typelist)[1:11] <- c("RM1","RM2","Mono","LAM", "YS","BM","cDC1","cDC2","pDC","Mig","Neu")
for (i in 1:length(Typelist)) {
  Typelist[[i]]$group <- names(Typelist[i])
}

plot_dat <- dplyr::bind_rows(Typelist)

colours=c("#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388","#4BB394","#46806F","#D95F02","#8FF7D8","#7570B3")
My_levels <- c("RM1","RM2","Mono","LAM", "YS","BM","cDC1","cDC2","pDC","Mig","Neu")
plot_dat$group<-factor(plot_dat$group, levels= My_levels)

#CDAHFD
ggplot(plot_dat)+ 
  geom_jitter(data =plot_dat[plot_dat$change == "Not significant",], 
              aes(group, avg_log2FC), color = "lightgrey",
              size=0.85, width = 0.4, alpha= .8)+
  geom_jitter(data = plot_dat[plot_dat$change != "Not significant",], 
              aes(group, avg_log2FC, color = group),
              size=0.85, width = 0.4, alpha= .8)+
  geom_tile(aes(group, 0, fill = group),
            height=0.8,
            color = "black",
            alpha = 0.5,
            show.legend = F,
            width=0.85) +
  geom_text_repel(data = top_label, aes(x =group, y =avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 3,
                  box.padding = unit(0.3, "lines")) 
  geom_text(data = plot_dat[!duplicated(plot_dat$group), ], 
            aes(group, 0, label = group),
            size =2,
            color ="black") +
  scale_fill_manual(values = colours)+
  scale_color_manual(name = "Cell subtype", values = colours2)+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+NoLegend()
ggsave("CDAHFD_myeloid_label.pdf", height = 7, width = 8)

# BDL
ggplot(plot_dat)+ 
  geom_jitter(data =plot_dat[plot_dat$change == "Not significant",], 
              aes(group, avg_log2FC), color = "lightgrey",
              size=0.85, width = 0.4, alpha= .8)+
  geom_jitter(data = plot_dat[plot_dat$change != "Not significant",], 
              aes(group, avg_log2FC, color = group),
              size=0.85, width = 0.4, alpha= .8)+
  geom_tile(aes(group, 0, fill = group),
            height=0.8,
            color = "black",
            alpha = 0.5,
            show.legend = F,
            width=0.85) +
  geom_text_repel(data = top_label, aes(x =group, y =avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 3,
                  box.padding = unit(0.3, "lines")) +
  geom_text(data = plot_dat[!duplicated(plot_dat$group), ], 
            aes(group, 0, label = group),
            size =2,
            color ="black") +
  scale_fill_manual(values = colours)+
  scale_color_manual(name = "Cell subtype", values = colours)+
  theme_bw()+
  theme(axis.text.x = element_text( vjust = 0.5),
        legend.position = "top")+coord_flip()+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+NoLegend()
ggsave("BDL_myeloid_label.pdf", height = 7, width = 8)

#radar charts
####For example BD-LAM####################################################
coord_radar <- function (theta = "x", start = 0, direction = 1) 
{  theta <- match.arg(theta, c("x", "y"))
r <- if (theta == "x") 
  "y"
else "x"
ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
        direction = sign(direction),
        is_linear = function(coord) TRUE)}

NASH<-read.table("~/LC-Myeloid-GSEA/NASH_LAM_common.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
BDL<-read.table("~/LC-Myeloid-GSEA/BDL_LAM_common.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
NASH<-NASH[order(NASH$pvalue),]
idlist<-NASH$ID
BDL<-BDL[match(idlist,BDL$ID),]
mydata<-rbind(NASH,BDL)

mydata$Description=as.factor(mydata$Description)
color=c( "#2171B5","red","red", "red","red","red", 
         "red","red", "red","red","red","red","red",
         "#238B45","red", "red","red","red", "red",
         "red","red", "red","red","red","red", "red")
color2=c( "black", "black","black","black","black","black",
          "black","black","black","black", 
          "black","black", "black","black",
          "black","black","black","black", 
          "black","black", "black","black",
          "black","black","black","black")
P<-ggplot(data=mydata,aes(x=forcats::fct_inorder(ID), y=-log10(pvalue))) + 
  geom_polygon(aes(x=forcats::fct_inorder(ID), y=-log10(pvalue),fill=NES,color = Description), size=2,mydata,group=mydata$Description,fill=color,color=color,alpha=0.2)+
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdGy")))+
  
  geom_point(aes(x=forcats::fct_inorder(ID), y=-log10(pvalue),fill=NES), mydata,color="white",shape=21,size=4)+
  
  coord_radar()+
  theme_bw()+
  ylim(0,4)+
  theme(axis.text.x=element_text(size = 11,colour="black"),
        axis.title=element_text(size=15,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        panel.grid.major = element_line(color="grey80"),
        axis.line = element_blank(),
        axis.ticks =  element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x  = element_text(size = 8,colour = "black"))
P+geom_polygon(aes(x=forcats::fct_inorder(ID), y=-log10(0.05),fill="black",color = "black"), size=1,mydata,group=mydata$Description,fill=color2,color=color2,alpha=0.3)
ggsave("LAM-common-gsea-7-7.pdf", height = 7, width = 7)

##LAM## unique wikipathway
coord_radar <- function (theta = "x", start = 0, direction = 1) 
{  theta <- match.arg(theta, c("x", "y"))
r <- if (theta == "x") 
  "y"
else "x"
ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
        direction = sign(direction),
        is_linear = function(coord) TRUE)}

NASH<-read.table("~/NASH_LAM_uni.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
BDL<-read.table("~/BDL_LAM_uni.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
NASH<-NASH[order(NASH$pvalue),]
idlist<-NASH$ID
BDL<-BDL[match(idlist,BDL$ID),]
mydata<-rbind(NASH,BDL)

mydata$Description=as.factor(mydata$Description)
color=c( "#2171B5","red","red", "red","red","red", "red", 
         "red","red", "red","red","red","red", "red","red","red", 
         "#238B45","red", "red","red","red", "red","red","red","red", "red","red", "red",
         "red","red", "red","red")
color2=c( "black", "black","black","black", "black","black","black", "black","black","black", "black","black","black", "black",
          "black","black", "black","black",
          "black","black","black","black", 
          "black","black", "black","black",
          "black","black","black","black", "black","black")
P<-ggplot(data=mydata,aes(x=forcats::fct_inorder(ID), y=-log10(pvalue))) + 
  geom_polygon(aes(x=forcats::fct_inorder(ID), y=-log10(pvalue),fill=NES,color = Description), size=2,mydata,group=mydata$Description,fill=color,color=color,alpha=0.2)+
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdGy")))+
  
  geom_point(aes(x=forcats::fct_inorder(ID), y=-log10(pvalue),fill=NES), mydata,color="white",shape=21,size=4)+
  
  coord_radar()+
  theme_bw()+
  ylim(0,4)+
  theme(axis.text.x=element_text(size = 11,colour="black"),
        axis.title=element_text(size=15,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        panel.grid.major = element_line(color="grey80"),
        axis.line = element_blank(),
        axis.ticks =  element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x  = element_text(size = 8,colour = "black"))
P+geom_polygon(aes(x=forcats::fct_inorder(ID), y=-log10(0.05),fill="black",color = "black"), size=1,mydata,group=mydata$Description,fill=color2,color=color2,alpha=0.3)

ggsave("LAM-uni-gsea-7-7.pdf", height = 7, width = 7)


