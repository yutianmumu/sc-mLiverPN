---
title: "Figure S5"
author: "Lin LEI"
---
setwd("~/FigureS5")
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(DT)
library(tidyverse)
library(patchwork)
library(harmony)
library(patchwork)
library(data.table)
library(tidyverse)
library(ggsci)
library(reshape2)
library(scRNAtoolVis)
library(ggunchull)
set.seed(123)

#load the integrated datasets object
integration <- readRDS("~/Figure1_new/integration_new.rds")
#subseting Myeloid populations
Endothelial<-subset(integration, idents =c("LSECs",
                                    "Portal vein ECs",
                                    "Lymphatic ECs",
                                    "Arterial ECs"))
p2 <- DimPlot(Endothelial, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)

a<-ls()
rm(list = a[which(a!='integration')])

scMerge_QC<- NormalizeData(Endothelial) %>% FindVariableFeatures() 
?CaseMatch
library(stringr)
s.genes.m<-str_to_title(cc.genes$s.genes)
g2m.genes.m<-str_to_title(cc.genes$g2m.genes)
CaseMatch(c(s.genes.m,g2m.genes.m),VariableFeatures(scMerge_QC))
#cell cycle evaluation
g2m_genes=CaseMatch(search = g2m.genes.m,match = rownames(scMerge_QC))
s_genes=CaseMatch(search = s.genes.m,match = rownames(scMerge_QC))
scMerge_QC<-CellCycleScoring(object = scMerge_QC,g2m.features = g2m_genes,s.features = s_genes)

#reclustering and de-batch effect
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
## Explicitly close multisession workers, if they were used

Endothelial_QC<-ScaleData(scMerge_QC,vars.to.regress=c("S.Score","G2M.Score"),features = rownames(scMerge_QC))

plan("sequential")

Endothelial_QC<- RunPCA(Endothelial_QC,features = VariableFeatures(Endothelial_QC),verbose=FALSE)

ElbowPlot(Endothelial_QC)
pc.num=1:20

Endothelial_Har2<- Endothelial_QC%>% RunHarmony("orig.ident", plot_convergence = TRUE, theta = 1)

Endothelial_Har2<- Endothelial_Har2 %>% RunUMAP(reduction = "harmony", dims = 1:40) %>% FindNeighbors(reduction = "harmony", dims = 1:40) %>% FindClusters(resolution = seq(from=0.1, to=0.4, by=0.1)) %>% identity()

p2 <- DimPlot(Endothelial_Har2, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)+ scale_colour_manual( values = colours)#group.by = "RNA_snn_res.0.3"

endo<-SetIdent(Endothelial_Har2,value=Endothelial_Har2@meta.data$RNA_snn_res.0.4)

#check the data with cell-type specific markers
#Featureplot
P2<-FeatureCornerAxes(object = endo,reduction = 'umap',
                      groupFacet = NULL,
                      legendPos = 'right',
                      features = c("Pecam1","Igfbp3","Adgrg6","Reln","Adam23","Wnt9b","Mki67"),
                      relLength = 0.5,relDist = 0.1,nLayout = 4,
                      lineTextcol = 'black',
                      pSize = 0.5)&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),limits=c(0,4))
ggsave(filename = "endo_marker_Umap_new.pdf", plot = P2, device="pdf", width=12, height=6)

#Define the cell clusters
new.cluster.ids <- c("Peri-portal LSECs", 
                     "Peri-central LSECs", 
                     "Portal vein ECs",
                     "Lymphatic ECs", 
                     "Arterial ECs"
)
names(new.cluster.ids) <- levels(endo)
endo<- RenameIdents(endo, new.cluster.ids)

p2 <- DimPlot(endo, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)

#reorder cluster
endo@active.ident <- factor(endo@active.ident, 
                            levels=c("Arterial ECs",
                                     "Portal vein ECs",
                                     "Lymphatic ECs", 
                                     "Peri-portal LSECs",
                                     "Peri-central LSECs"
                            ))
endo$cell_type<-endo@active.ident

#Proportions of each cell cluster in total endothelial cell for each treatment
table(endo$treatment)
My_levels <- c('Control','CDAHFD','BDL')
endo$treatment<-factor(endo$treatment, levels= My_levels)

prop.table(table(Idents(endo)))
table(Idents(endo),endo$treatment)
Cellratio<-prop.table(table(Idents(endo),endo$treatment),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1,levels=c("Arterial ECs",
                                                 "Portal vein ECs",
                                                 "Lymphatic ECs", 
                                                 "Peri-portal LSECs",
                                                 "Peri-central LSECs"))

table(Cellratio$Var1)
colourCount=length(unique(Cellratio$Var1))
cells_in_treatment<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  theme_classic()+labs(x="Sample",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values = c( "#80B1D3","#B15928",  "#FB8072","#6A3D9A","darkgrey" ))
ggsave(filename = "endo_cluster_Freq_vertical_new.pdf", plot = cells_in_treatment, device="pdf", width=6, height=7.8)

#Proportions of each cell cluster in each treatment condition
table(integration$orig.ident)
b<-as.data.frame(table(Idents(endo),endo$orig.ident)) 
b$Freq2=0
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
            file="Cellratio_Tot_new.txt", sep="\t",row.names=FALSE, quote=FALSE)

CellRa<-read.table("Cellratio_Tot_new.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
CellRa$group <- factor(CellRa$group,levels=c("Control","CDAHFD","BDL"))

p3<-ggplot(data=CellRa[CellRa$Var1=="Peri-portal LSECs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,25)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Peri-portal LSECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

p4<-ggplot(data=CellRa[CellRa$Var1=="Portal vein ECs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,7.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Portal vein ECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

p5<-ggplot(data=CellRa[CellRa$Var1=="Lymphatic ECs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Lymphatic ECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

p6<-ggplot(data=CellRa[CellRa$Var1=="Arterial ECs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Arterial ECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

p7<-ggplot(data=CellRa[CellRa$Var1=="Peri-central LSECs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,12)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Peri-central LSECs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")


library(ggpubr)
p10<-ggarrange(p6,p4,p5,p7,p3,ncol = 4,nrow = 2, legend = "top",common.legend = T)
ggsave(filename = "endo_each_Freq_Total_new.pdf", plot = p10, device="pdf", width=15, height=8)


#Cluster DEGs
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
#Explicitly close multisession workers, if they were used
obj.markers <- FindAllMarkers(endo, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
write.table(obj.markers, file="obj.markers.txt", sep="\t", row.names=FALSE, quote=FALSE)
plan("sequential")

top10 = obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, file="Top10DEGs.txt", sep="\t", row.names=FALSE, quote=FALSE)
#heatmap
colours=c( "#80B1D3","#B15928",  "#FB8072", "#6A3D9A","darkgrey" )
cols.use <- list(cell_type=colours,
                 treatment=c("#737373","#2171B5","#238B45"))
endo@meta.data$treatment<- factor(x = endo@meta.data$treatment, levels=c('Control', 'CDAHFD', 'BDL'))

endo$cell_type<-endo@active.ident
table(endo$treatment)
My_levels <- c('Control','CDAHFD','BDL')
endo$treatment<-factor(endo$treatment, levels= My_levels)
library(grid)
Heatmap<-DoMultiBarHeatmap(endo, features = top10$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme(axis.text.y = element_text(face="italic",color = 'black'))
ggsave(filename = "endo_heatmap_addanno_new.pdf", plot = Heatmap, device="pdf", width=6, height=6.6)

#load genesets
endo_geneset <- readxl::read_xlsx("I:~/endo_geneset.xlsx")
Fibrinolysis<-as.character(endo_geneset$fibrinolysis)
Fibrinolysis<-Fibrinolysis[1:5]
Profibrogenesis<-as.character(endo_geneset$profibrogenesis)
Profibrogenesis<-Profibrogenesis[1:19]
#geneset score caculation_Fibrinolysis
endo<-AddModuleScore(object = endo,features=list(Fibrinolysis),ctrl = 100,name='Fibrinolysis')
FeaturePlot(object =endo, features = "Fibrinolysis1",split.by = "treatment", reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

library(ggpubr)
boxplot<-ggboxplot(endo@meta.data, x = "cell_type", y = "Fibrinolysis1",
                   color = "black",fill = 'treatment',palette=c("#737373","#2171B5","#238B45"),
                   bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))+ labs(x = '')+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype = "dashed", color = "black")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
ggsave(filename = "Fibrinolysis_score_boxplot.pdf", plot = boxplot, device="pdf", width=8.5, height=4)
#determine P value
data<-endo@meta.data
library(rstatix)
wilcox_res<-data%>%
  group_by(cell_type) %>%
  wilcox_test(Fibrinolysis1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="Fibrinolysis_wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)
#geneset score caculation_Profibrogenesis
endo<-AddModuleScore(object = endo,features=list(Profibrogenesis),ctrl = 100,name='Profibrogenesis')
FeaturePlot(object =endo, features = "Profibrogenesis1",split.by = "treatment", reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

boxplot<-ggboxplot(endo@meta.data, x = "cell_type", y = "Profibrogenesis1",
                   color = "black",fill = 'treatment',palette=c("#737373","#2171B5","#238B45"),
                   bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))+ labs(x = '')+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype = "dashed", color = "black")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
ggsave(filename = "Profibrogenesis_score_boxplot.pdf", plot = boxplot, device="pdf", width=8.5, height=4)
data<-endo@meta.data
library(rstatix)
wilcox_res<-data%>%
  group_by(cell_type) %>%
  wilcox_test(Profibrogenesis1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="Profibrogenesis_wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)

#Go function analysis
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(viridis)
#load DEGs list
markers<-read.table("obj.markers_Control_new.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
View(markers)
#select DEGs for Go analysis
sigposDEG.all<- subset(markers, p_val_adj<0.01&abs(avg_log2FC)>1) 

Arter_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Arterial ECs') 
Por_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Portal vein ECs') 
Lym_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Lymphatic ECs')
Peripor_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Peri-portal LSECs') 
Pericen_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Peri-central LSECs') 

list<-list(Arter_up,Por_up,Lym_up,Peripor_up,Pericen_up)
names(list)[1:5] <- c("Arterial ECs",
                      "Portal vein ECs",
                      "Lymphatic ECs", 
                      "Peri-portal LSECs",
                      "Peri-central LSECs")

for(i in 1:length(list)){GO <- enrichGO(gene =list[[i]]$gene,
                                        #universe = row.names(dge.celltype),
                                        OrgDb = 'org.Mm.eg.db',
                                        keyType = 'SYMBOL',
                                        ont = "BP", 
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.01,
                                        qvalueCutoff = 0.01)
GO<- data.frame(GO)
write.csv(GO,paste0("GO_new_", names(list[i]),".CSV"))
} 

list2<-list(Arter_up$gene,Por_up$gene,Lym_up$gene,Pro_up$gene,Peripor_up$gene,Pericen_up$gene)
names(list2)[1:6] <- c("Arterial ECs","Portal vein ECs","Lymphatic ECs", "Prolif.ECs" ,"Peri-portal LSECs","Peri-central LSECs")
ck <- compareCluster(geneCluster = list2[1:6], fun = "enrichGO", OrgDb = 'org.Mm.eg.db',keyType = 'SYMBOL')

ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="SYMBOL")#keyType="SYMBOL"
table(ck@compareClusterResult$Cluster)
t(ck@compareClusterResult[1,])
dotplot(ck)

#doublets checking
#devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)
#PK value identification
sweep.res.list <- paramSweep_v3(endo, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() 
DoubletRate = 0.03 
homotypic.prop <- modelHomotypic(endo$cell_type) 
nExp_poi <- round(DoubletRate*ncol(endo))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gc()
endo_wo_doublets <- doubletFinder_v3(endo, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
doubletsendo<-DimPlot(endo_wo_doublets,size=1.5, reduction = "umap",split.by = "treatment", group.by = "DF.classifications_0.25_0.005_146")
ggsave(filename = "doubletsendo.pdf", plot = doubletsendo, device="pdf", width=10, height=4)


library("Nebulosa")
Nash<- endo[,endo@meta.data$orig.ident %in% c("CDAHFD_1","CDAHFD_2")]
Control<- endo[,endo@meta.data$orig.ident %in% c("Control_1","Control_2","Control_3")]
BDL<- endo[,endo@meta.data$orig.ident %in% c("BDL_1","BDL_2")]
#Ackr1
p1 <- plot_density(Control, c("Cd34", "Ackr1"),joint = TRUE,size = 0.5)&scale_color_distiller(limits=c(-0.0001,0.011),palette = "Spectral")
p2 <- plot_density(Nash, c("Cd34", "Ackr1"),joint = TRUE,size = 0.5)&scale_color_distiller(limits=c(-0.0001,0.011),palette = "Spectral")
p3 <- plot_density(BDL, c("Cd34","Ackr1"),joint = TRUE,size = 0.5)&scale_color_distiller(limits=c(-0.0001,0.011),palette = "Spectral")
p1[[3]]
library(ggpubr)
Ackr1_density<-ggarrange(p1[[3]],p2[[3]],p3[[3]],ncol = 3,nrow = 1,common.legend = TRUE,legend = "bottom")
ggsave(filename = "Cd34+Ackr1+_density.pdf", plot = Ackr1_density, device="pdf", width=8.5, height=3.5)

