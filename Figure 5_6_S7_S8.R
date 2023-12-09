library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(DT)
library(tidyverse)
library(patchwork)
library(harmony)
library(ggpubr)
library(reshape2)
library(scRNAtoolVis)
library(ggunchull)
library(viridis)
set.seed(123)
integration <- readRDS("~/Figure1_new/integration_new.rds")
#subseting Myeloid populations
Mesenchymal<-subset(integration, idents =c(
  "Fibrolasts & myofibroblasts", 
  "HSCs"))
p2 <- DimPlot(Myeloid_cell, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)
levels(Mesenchymal)
control<- Mesenchymal[,Mesenchymal@meta.data$orig.ident %in% c("control_1","control_2","control_3")]
control<-subset(control, idents =c( "2", "4","5","6","7","9"))
levels(control)
#Load the dataset from Lei et al 
load("~/Mesenchymal_cells/Lei et al.RData")
levels(lei_health)
lei<-subset(x = lei_health, idents = c("Fib-1",
                                      "Fib-2",
                                      "Fib-3",
                                      "Fib-4",
                                      "Fib-5",
                                      "HSC"))
DefaultAssay(lei) <- "RNA"

#Integrating datasets of control and lei et al
lei$cell_type<-lei@active.ident
control$cell_type<-control@active.ident
scMerge<- merge(lei, y=c(control))
scMerge_QC<- NormalizeData(scMerge) %>% FindVariableFeatures() 

library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
scMerge_QC<-ScaleData(scMerge_QC,features = rownames(scMerge_QC))

plan("sequential")
scMerge_QC<- RunPCA(scMerge_QC,features = VariableFeatures(scMerge_QC),verbose=FALSE)

scMerge_Har2<- scMerge_QC%>% RunHarmony("orig.ident", plot_convergence = TRUE, theta = 1)
rm(scMerge_QC)

scMerge_Har2<- scMerge_Har2 %>% RunUMAP(reduction = "harmony", dims = 1:40) %>% FindNeighbors(reduction = "harmony", dims = 1:40) %>% FindClusters(resolution = seq(from=0.1, to=0.4, by=0.1)) %>% identity()
scMerge_Har2$cell_type<-factor(scMerge_Har2$cell_type)
scMerge_Har2@active.ident<-scMerge_Har2$cell_type
DimPlot(scMerge_Har2, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1.5)
colours = c(brewer.pal(12,"Paired")[1:3],brewer.pal(12,"Paired")[5:12],brewer.pal(8,"Dark2")[1:4],brewer.pal(8,"Dark2")[6:8],brewer.pal(11,"RdBu")[4:5])
p2 <- DimPlot(p2, reduction = "umap", label = TRUE, label.size = 6,repel = TRUE,pt.size = 1)+ scale_colour_manual( values = colours)
ggsave(filename = "control_lei_etal.pdf", plot = p2, device="pdf", width=6, height=4.5)

#######correlation analysis##########
table(scMerge_Har2$cell_type)
av <-AverageExpression(scMerge_Har2,
                       group.by = "cell_type",
                       assays = "integrated")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av, 1, sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = 'spearman'))
cor<-cor(av[cg,])
cor2<-pheatmap::pheatmap(cor(av[cg,],method = 'spearman',), display_numbers = ifelse(cor > 0.94, "*"," "),    
                        number_color = "black",        
                        fontsize=12,                    
                        color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(15))
ggsave(filename = "celltype_decision_correlation_heatmap.pdf", plot = cor2, device="pdf", width=6, height=5)
#cell cluster annotation
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
mesen_markers <- FindAllMarkers(Mesenchymal, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0)
write.table(mesen_markers, file="Mesen_markers.txt", sep="\t", row.names=FALSE, quote=FALSE)
plan("sequential")

top10 = mesen_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, file="Top10DEGs.txt", sep="\t", row.names=FALSE, quote=FALSE)
FeaturePlot(object = Mesenchymal, features = "Col6a5",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(object = Mesenchymal, features = "Col15a1",min.cutoff =0, max.cutoff = 2,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(object = Mesenchymal, features = "Pla1a",min.cutoff =0, max.cutoff = 2,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(object = Mesenchymal, features = "Meox1",min.cutoff =0, max.cutoff = 2,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(object = Mesenchymal, features = "Aldh1a2",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(object = Mesenchymal, features = "Pbk",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(object = Mesenchymal, features = "H1f5",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

#Rename the cell clusters
levels(Mesenchymal)
new.cluster.ids <- c("Col6a5+ MFB", 
                     "Meox1+ MFB", 
                     "Fib-1",
                     "Aldh1a2+ MFB",
                     "Fib-4", 
                     "Fib-5",
                     "Fib-2",
                     "HSCs",
                     "Robo2+ MFB",
                     "Fib-3",
                     "Prolif.MFB-1",
                     "Prolif.MFB-2")
names(new.cluster.ids) <- levels(Mesenchymal)
Mesenchymal<- RenameIdents(Mesenchymal, new.cluster.ids)
DimPlot(Mesenchymal, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)
Mesenchymal@active.ident <- factor(Mesenchymal@active.ident, 
                                   levels=c("Fib-1",
                                            "Fib-2",
                                            "Fib-3",
                                            "Fib-4", 
                                            "Fib-5",
                                            "HSCs",
                                            "Robo2+ MFB",
                                            "Col6a5+ MFB", 
                                            "Meox1+ MFB", 
                                            "Aldh1a2+ MFB",
                                            "Prolif.MFB-1",
                                            "Prolif.MFB-2"))

colour=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "navy", "#6A3D9A",  "#CAB2D6", "#E31A1C", "#FDBF6F", "#FF7F00", "#FFFF99","#FB9A99")
p3 <- DimPlot(Mesenchymal, reduction = "umap", label = T, repel = TRUE,pt.size = 1) +scale_colour_manual(values =colour ) 
ggsave(filename = "Mesenchymal_Umap.pdf", plot = p3, device="pdf", width=7, height=4.5)

colour2=c("#E31A1C",
          "#FDBF6F",
          "#A6CEE3", 
          "#FF7F00", 
          "#33A02C",
          "navy",
          "#1F78B4", 
          "#6A3D9A",
          "#CAB2D6",
          "#B2DF8A",
          "#FFFF99",
          "#FB9A99")
p3 <- DimPlot(Mesenchymal, reduction = "umap",group.by = "RNA_snn_res.0.8", label = T, repel = TRUE,pt.size = 1) +scale_colour_manual(values =colour2 ) 
ggsave(filename = "Mesenchymal_notdefined_Umap.pdf", plot = p3, device="pdf", width=5.5, height=4.5)
Mesenchymal$cell_type<-Mesenchymal@active.ident

saveRDS(Mesenchymal, file = "~/Mesenchymal_cells/Mesenchymal_0630.rds")
###########################################load RDS#################################################
###########################################################################################################
#plot UMAP of marker expression
Mesenchymal <- readRDS("~/Mesenchymal_cells/Mesenchymal_0630.rds")
TOP10 <- read.table("Top10DEGs.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
marker <- read.table("Mesen_markers.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)

p1=FeaturePlot(object = Mesenchymal, features = "Pdgfra",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p1=p1+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p2=FeaturePlot(object = Mesenchymal, features = "Dcn",min.cutoff =0, max.cutoff = 5.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p2=p2+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p3=FeaturePlot(object = Mesenchymal, features = "Pdgfrb",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p3=p3+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p4=FeaturePlot(object = Mesenchymal, features = "Dpt",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p4=p4+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p5=FeaturePlot(object = Mesenchymal, features = "Col1a1",min.cutoff =0, max.cutoff = 4,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p5=p5+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p6=FeaturePlot(object = Mesenchymal, features = "Col15a1",min.cutoff =0, max.cutoff = 2,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p6=p6+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p7=FeaturePlot(object = Mesenchymal, features = "Acta2",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p7=p7+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p8=FeaturePlot(object = Mesenchymal, features = "Mki67",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 0.5, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p8=p8+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p9<-ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow = 2,common.legend = T, legend = "right")+theme(plot.title = element_text(face="italic",color = 'black'))
ggsave(filename = "mesenchymal_UMAP.pdf", plot = p9, device="pdf", width=16, height=8)

p1=FeaturePlot(object = Mesenchymal, features = "Robo2",min.cutoff =0, max.cutoff = 3,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p1=p1+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()

p2=FeaturePlot(object = Mesenchymal, features = "Col6a5",min.cutoff =0, max.cutoff = 4,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p2=p2+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()

p3=FeaturePlot(object = Mesenchymal, features = "Meox1",min.cutoff =0, max.cutoff = 2,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p3=p3+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()

p4=FeaturePlot(object = Mesenchymal, features = "Aldh1a2",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p4=p4+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()

p5=FeaturePlot(object = Mesenchymal, features = "Col1a1",min.cutoff =0, max.cutoff = 5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p5=p5+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()


p7=FeaturePlot(object = Mesenchymal, features = "Acta2",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p7=p7+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()

p8=FeaturePlot(object = Mesenchymal, features = "Mki67",min.cutoff =0, max.cutoff = 2.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p8=p8+theme(plot.title = element_text(face="italic",color = 'black', size=24))+NoAxes()


p9<-ggarrange(p5,p7,p1,p2,p3,p4,p8,ncol = 4,nrow = 2,common.legend = T, legend = "right")+theme(plot.title = element_text(face="italic",color = 'black'))
ggsave(filename = "mesenchymal_marker_UMAP0907.pdf", plot = p9, device="pdf", width=16, height=8)

p10=FeaturePlot(object = Mesenchymal, features = "Slit2",min.cutoff =0, max.cutoff = 1,reduction = "umap",order = T, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p10=p10+theme(plot.title = element_text(face="italic",color = 'black'))
ggsave(filename = "Slit2.pdf", plot = p10, device="pdf", width=5.5, height=4.5)

#signature of PMSC/PMSC-MFB
sig.PMSCs=c("Col1a2","Col15a1","Igfbp6","Loxl1","Mgp","Slit2","Thy1")
sig.HSCs=c("Bmp10","Hgf","Masp1","Megf9","Nt5e","Plac8")
Mesenchymal<-AddModuleScore(object = Mesenchymal,features=list(c("Col1a2","Col15a1","Igfbp6","Loxl1","Mgp","Slit2","Thy1")),ctrl = 100,name='signature_PMSC')
sigPMSC<-FeaturePlot(object =Mesenchymal, features = "signature_PMSC1",reduction = "umap",order = F, blend = F, pt.size = 1,min.cutoff =0, max.cutoff = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
sigPMSC<-sigPMSC+labs(title="PMSC 7-gene signature")
#signature of HSC/HSC-MFB
Mesenchymal<-AddModuleScore(object = Mesenchymal,features=list(c("Bmp10","Hgf","Masp1","Megf9","Nt5e","Plac8")),ctrl = 100,name='signature_HSC')
sigHSC<-FeaturePlot(object =Mesenchymal, features = "signature_HSC1",reduction = "umap",order = F, blend = F, pt.size = 1,min.cutoff =0, max.cutoff = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
sigHSC<-sigHSC+labs(title="HSC 6-gene signature")
com<-sigPMSC|sigHSC
ggsave(filename = "signature_PMSC_HSC.pdf", plot = com, device="pdf", width=11, height=4.5)
#Heatmap
top10 <- read.table("Top10DEGs.txt",sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
colour=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "navy", "#6A3D9A",  "#CAB2D6", "#E31A1C", "#FDBF6F", "#FF7F00", "#FFFF99","#FB9A99")
cols.use <- list(cell_type=colour,
                 treatment=c("#737373","#2171B5","#238B45"))
Mesenchymal@meta.data$treatment<- factor(x = Mesenchymal@meta.data$treatment, levels=c('control', 'CDAHFD', 'BDL'))

Mesenchymal$cell_type<-Mesenchymal@active.ident
table(Mesenchymal$treatment)
My_levels <- c('control','CDAHFD','BDL')
Mesenchymal$treatment<-factor(Mesenchymal$treatment, levels= My_levels)
library(grid)
Heatmap<-DoMultiBarHeatmap(Mesenchymal, features = top10$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme(axis.text.y = element_text(face="italic",color = 'black'))+
  guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(filename = "Mesenchymal_heatmap_addanno_new2.pdf", plot = Heatmap, device="pdf", width=6, height=6.6)

#composition of clusters in mesenchymal cells
table(Mesenchymal$treatment)
My_levels <- c('control','CDAHFD','BDL')
Mesenchymal$treatment<-factor(Mesenchymal$treatment, levels= My_levels)

prop.table(table(Idents(Mesenchymal)))
table(Idents(Mesenchymal),Mesenchymal$treatment)
Cellratio<-prop.table(table(Idents(Mesenchymal),Mesenchymal$treatment),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1,levels=c("Fib-1",
                                                 "Fib-2",
                                                 "Fib-3",
                                                 "Fib-4", 
                                                 "Fib-5",
                                                 "HSCs",
                                                 "Robo2+ MFB",
                                                 "Col6a5+ MFB", 
                                                 "Meox1+ MFB", 
                                                 "Aldh1a2+ MFB",
                                                 "Prolif.MFB-1",
                                                 "Prolif.MFB-2"))
table(Cellratio$Var1)
colourCount=length(unique(Cellratio$Var1))
colour=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "navy", "#6A3D9A",  "#CAB2D6", "#E31A1C", "#FDBF6F", "#FF7F00", "#FFFF99","#FB9A99")
cells_in_treatment<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  theme_classic()+labs(x="Sample",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values = colour)
ggsave(filename = "Mesenchymal_cluster_Freq_vertical_new.pdf", plot = cells_in_treatment, device="pdf", width=6, height=7.8)

#highlight cells in each condition
Mesenchymal2<-SetIdent(Mesenchymal,value=Mesenchymal@meta.data$treatment)
hi1<-WhichCells(Mesenchymal2,idents = "control")
hi2<-WhichCells(Mesenchymal2,idents = "CDAHFD")
hi3<-WhichCells(Mesenchymal2,idents = "BDL")
sham_UMAP<-DimPlot(Mesenchymal,group.by = "treatment", reduction = "umap",cells.highlight= hi1, cols.highlight = c("#737373"),sizes.highlight = 0.5, label = F,repel = F,  pt.size = 1,cols = c("#BDBDBD"))
ggsave(filename = "shamUMAP_Mesenchymal_distrbution_new.pdf", plot = sham_UMAP, device="pdf", width=6.2, height=4.5)

CDAHFD_UMAP<-DimPlot(Mesenchymal,group.by = "treatment", reduction = "umap",cells.highlight= hi2, cols.highlight = c("#2171B5"),sizes.highlight = 0.5, label = F,repel = F,  pt.size = 1,cols = c("#BDBDBD"))
ggsave(filename = "CDAHFD_UMAP_Mesenchymal_distrbution_new.pdf", plot = CDAHFD_UMAP, device="pdf", width=6.2, height=4.5)

BDL_UMAP<-DimPlot(Mesenchymal,group.by = "treatment", reduction = "umap",cells.highlight= hi3, cols.highlight = c("#238B45"),sizes.highlight = 0.5, label = F,repel = F,  pt.size = 1,cols = c("#BDBDBD"))
ggsave(filename = "BDL_UMAP_Mesenchymal_distrbution_new.pdf", plot = BDL_UMAP, device="pdf", width=6.2, height=4.5)
library(ggpubr)
P4<-ggarrange(sham_UMAP,CDAHFD_UMAP,BDL_UMAP,ncol = 3,nrow = 1, legend = "top")
ggsave(filename = "Mesenchymal_distrbution_UMAP_new2.pdf", plot = P4, device="pdf", width=12, height=4.5)

#Proportions of each cell cluster in each treatment condition

integration <- readRDS("~/Figure1_new/integration_new.rds")
table(integration$orig.ident)

b<-as.data.frame(table(Idents(Mesenchymal),Mesenchymal$orig.ident)) 
b$Freq2=0
for (i in 1:nrow(b) ) {
  if(b$Var2[i]=="control_1"){b$Freq2[i]=b$Freq[i]/3493
  }
  else if (b$Var2[i]=="control_2"){b$Freq2[i]=b$Freq[i]/3990
  }
  else if (b$Var2[i]=="control_3"){b$Freq2[i]=b$Freq[i]/1630
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
CellRa$group <- factor(CellRa$group,levels=c("control","CDAHFD","BDL"))

#Fib-1
p3<-ggplot(data=CellRa[CellRa$Var1=="Fib-1",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,6)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Fib-1") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
# Fib-2
p4<-ggplot(data=CellRa[CellRa$Var1=="Fib-2",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,3)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Fib-2") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Fib-3
p5<-ggplot(data=CellRa[CellRa$Var1=="Fib-3",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Fib-3") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Fib-4
p6<-ggplot(data=CellRa[CellRa$Var1=="Fib-4",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,2.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Fib-4") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Fib-5
p7<-ggplot(data=CellRa[CellRa$Var1=="Fib-5",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,2)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Fib-5") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

#HSC
p8<-ggplot(data=CellRa[CellRa$Var1=="HSCs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,2)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="HSCs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Robo2+ MFB
p9<-ggplot(data=CellRa[CellRa$Var1=="Robo2+ MFB",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Robo2+ MFB") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

#Col6a5+ MFB
p10<-ggplot(data=CellRa[CellRa$Var1=="Col6a5+ MFB",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,8)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Col6a5+ MFB") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Aldh1a2+ MFB
p11<-ggplot(data=CellRa[CellRa$Var1=="Aldh1a2+ MFB",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,6)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Aldh1a2+ MFB") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Prolif.MFB-1
p12<-ggplot(data=CellRa[CellRa$Var1=="Prolif.MFB-1",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,0.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Prolif.MFB-1") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Prolif.MFB-2
p13<-ggplot(data=CellRa[CellRa$Var1=="Prolif.MFB-2",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,0.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Prolif.MFB-2") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")


library(ggpubr)
p14<-ggarrange(p3,p4,p5,p6,p7,p9,p10,p11,p12,p13,ncol = 5,nrow = 2, legend = "top",common.legend = T)

ggsave(filename = "Mesenchymal_each_Freq_Total_new.pdf", plot = p14, device="pdf", width=18.75, height=8)

#go function analysis
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(viridis)
Mesenchymal <- readRDS("~/Mesenchymal_cells/Mesenchymal_0630.rds")
control<- Mesenchymal[,Mesenchymal@meta.data$orig.ident %in% c("control_1","control_2","control_3")]
control<-subset(x =control, idents = c("Fib-1","Fib-2","Fib-3","Fib-4","Fib-5","HSCs"))
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
control_markers <- FindAllMarkers(control, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0)
write.table(control_markers, file="control_markers.txt", sep="\t", row.names=FALSE, quote=FALSE)
plan("sequential")
markers<-read.table("control_markers.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
View(markers)
sigposDEG.all<- subset(markers, p_val_adj<=1&abs(avg_log2FC)>0.25)
Fib1_up <- subset(sigposDEG.all, p_val_adj<0.05&avg_log2FC>0.25&cluster=='Fib-1') 
Fib2_up <- subset(sigposDEG.all, cluster=='Fib-2') 
Fib3_up <- subset(sigposDEG.all, p_val_adj<0.05&avg_log2FC>0.25&cluster=='Fib-3') 
Fib4_up <- subset(sigposDEG.all, p_val_adj<0.05&avg_log2FC>0.25&cluster=='Fib-4') 
Fib5_up <- subset(sigposDEG.all, p_val_adj<0.05&avg_log2FC>0.25&cluster=='Fib-5') 
HSC_up <- subset(sigposDEG.all, p_val_adj<0.05&avg_log2FC>0.25&cluster=='HSCs') 

list<-list(Fib1_up,Fib2_up,Fib3_up,Fib4_up,Fib5_up,HSC_up)
names(list)[1:6] <- c("Fib-1",
                       "Fib-2",
                       "Fib-3",
                       "Fib-4", 
                       "Fib-5",
                       "HSCs"
                       )

for(i in 1:length(list)){GO <- enrichGO(gene =list[[i]]$gene,
                                        OrgDb = 'org.Mm.eg.db',
                                        keyType = 'SYMBOL',
                                        ont = "BP", 
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.01,
                                        qvalueCutoff = 0.01)
GO<- data.frame(GO)
write.csv(GO,paste0("GO_new_", names(list[i]),".CSV"))
} 

Fib1 <- read.csv("GO_new_Fib-1.CSV", row.names = 1)
Fib2 <- read.csv("GO_new_Fib-2.CSV", row.names = 1)
Fib3 <- read.csv("GO_new_Fib-3.CSV", row.names = 1)
Fib4 <- read.csv("GO_new_Fib-4.CSV", row.names = 1)
Fib5 <- read.csv("GO_new_Fib-5.CSV", row.names = 1)
HSC <- read.csv("GO_new_HSCs.CSV", row.names = 1)
Fib1$group<-"Fib-1"
Fib2$group<-"Fib-2"
Fib3$group<-"Fib-3"
Fib4$group<-"Fib-4"
Fib5$group<-"Fib-5"
HSC$group<-"HSCs"
Fib1<-Fib1[1:20,]
Fib2<-Fib2[1:20,]
Fib3<-Fib3[1:20,]
Fib4<-Fib4[1:20,]
Fib5<-Fib5[1:20,]
HSC<-HSC[1:20,]
all <- rbind(Fib1,Fib2,Fib3,Fib4,Fib5,HSC)
all$Description<- gsub("-"," ",all$Description)
library(forcats)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
dat2<-str_split(all$GeneRatio,"/",simplify = T)[,1]
all$Num<-as.numeric(dat2) 
dat3<-str_split(all$GeneRatio,"/",simplify = T)[,2]
all$Tot<-as.numeric(dat3) 
all$Generatio<-as.numeric(all$Num/all$Tot)
My_levels <- c("Fib-1",
               "Fib-2",
               "Fib-3",
               "Fib-4", 
               "Fib-5",
               "HSCs")
all$group<-factor(all$group, levels= My_levels)
ggplot(all, aes(group, Description)) +
  geom_point(aes(color=qvalue, size=Generatio))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
selct<-list(
  "extracellular matrix organization",
  "external encapsulating structure organization",
  "ossification",
  "regulation of epithelial cell proliferation",
  "branching morphogenesis of an epithelial tube",
  "morphogenesis of a branching epithelium",
  "morphogenesis of a branching structure",
  "response to fibroblast growth factor",
  "maternal placenta development",
  "reproductive system development",
  "respiratory tube development",
  "ureteric bud development",
  "mesonephros development",
  "sprouting angiogenesis",
  "bone mineralization",
  "gland development",
  "nitric oxide biosynthetic process",
  "regulation of nitric oxide metabolic process",
  "euron projection extension",
  "T cell mediated cytotoxicity",
  "leukocyte mediated cytotoxicity",
  "lung development",
  "fat cell differentiation"
)

selctterm<-as.character(selct)
Fib1 <- read.csv("GO_new_Fib-1.CSV", row.names = 1)
Fib2 <- read.csv("GO_new_Fib-2.CSV", row.names = 1)
Fib3 <- read.csv("GO_new_Fib-3.CSV", row.names = 1)
Fib4 <- read.csv("GO_new_Fib-4.CSV", row.names = 1)
Fib5 <- read.csv("GO_new_Fib-5.CSV", row.names = 1)
HSC <- read.csv("GO_new_HSCs.CSV", row.names = 1)
Fib1$group<-"Fib-1"
Fib2$group<-"Fib-2"
Fib3$group<-"Fib-3"
Fib4$group<-"Fib-4"
Fib5$group<-"Fib-5"
HSC$group<-"HSCs"
Fib1$Description<- gsub("-"," ",Fib1$Description)
Fib2$Description<- gsub("-"," ",Fib2$Description)
Fib3$Description<- gsub("-"," ",Fib3$Description)
Fib4$Description<- gsub("-"," ",Fib4$Description)
Fib5$Description<- gsub("-"," ",Fib5$Description)
HSC$Description<- gsub("-"," ",HSC$Description)

Fib1=Fib1[Fib1$Description%in%selctterm,]
Fib2=Fib2[Fib2$Description%in%selctterm,]
Fib3=Fib3[Fib3$Description%in%selctterm,]
Fib4=Fib4[Fib4$Description%in%selctterm,]
Fib5=Fib5[Fib5$Description%in%selctterm,]
HSC=HSC[HSC$Description%in%selctterm,]

all <- rbind(Fib1,Fib2,Fib3,Fib4,Fib5,HSC)
library(forcats)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
dat2<-str_split(all$GeneRatio,"/",simplify = T)[,1]
all$Num<-as.numeric(dat2) 
dat3<-str_split(all$GeneRatio,"/",simplify = T)[,2]
all$Tot<-as.numeric(dat3) 
all$Generatio<-as.numeric(all$Num/all$Tot)
My_levels <- c("Fib-1",
               "Fib-2",
               "Fib-3",
               "Fib-4",
               "Fib-5",
               "HSCs")
My_levels2<-c("extracellular matrix organization",
              "external encapsulating structure organization",
              "regulation of epithelial cell proliferation",
              "ossification", 
              "morphogenesis of a branching structure",
              "gland development",
              "bone mineralization",
              "branching morphogenesis of an epithelial tube",
              "morphogenesis of a branching epithelium",
              "response to fibroblast growth factor",
              "ureteric bud development",
              "mesonephros development",
              "sprouting angiogenesis",
              "fat cell differentiation",
              "lung development",
              "maternal placenta development",
              "reproductive system development",
              "respiratory tube development",
              "nitric oxide biosynthetic process",
              "regulation of nitric oxide metabolic process",
              "euron projection extension",
              "T cell mediated cytotoxicity",
              "leukocyte mediated cytotoxicity")
all$group<-factor(all$group, levels= My_levels)
all$Description<-factor(all$Description, levels= My_levels2)
library(viridis)
MyeGO<-ggplot(all, aes(group, Description)) +
  theme_bw()+
  geom_point(aes(fill=qvalue, size=Generatio),shape=21,colour="black",alpha=0.8)+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_viridis(option = "D")

ggsave(filename = "mesenchymal_control_GO0906.pdf", plot = MyeGO, device="pdf", width=5.5, height=4.5)

##################################-----Monocle-----###################################
library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
set.seed(123) 

expr_matrix <- as(as.matrix(Mesenchymal@assays$RNA@counts), 'sparseMatrix') 
p_data <- Mesenchymal@meta.data
f_data <- data.frame(gene_short_name = row.names(Mesenchymal),row.names = row.names(Mesenchymal))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),
                                     num_cells_expressed >= 10)) 

######DEG TEST##########
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1)
head(diff)
deg <- subset(diff, qval < 0.05) 
deg <- deg1[order(deg$qval,decreasing=F),]

head(deg)
write.table(deg,file="train.monocle.DEG_mesenchymal.xls",col.names=T,row.names=F,sep="\t",quote=F)

ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
pdf("train.ordergenes_mesenchymal.pdf")
plot_ordering_genes(cds)
dev.off()

cds <- reduceDimension(cds, max_components = 2,
                        method = 'DDRTree')

cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 7)
library(RColorBrewer)
P1<-plot_cell_trajectory(cds, show_tree = T,show_backbone = F, show_branch_points  = FALSE,color_by = "cell_type")+scale_colour_manual(values=c("Fib-1" ="#A6CEE3",
                                                                                                                                               "Fib-2" ="#1F78B4",
                                                                                                                                               "Fib-3" ="#B2DF8A",
                                                                                                                                               "Fib-4" = "#33A02C",
                                                                                                                                               "Fib-5" = "navy",
                                                                                                                                               "HSCs" = "#6A3D9A",
                                                                                                                                               "Robo2+ MFB" = "#CAB2D6",
                                                                                                                                               "Col6a5+ MFB"="#E31A1C",
                                                                                                                                               "Meox1+ MFB"="#FDBF6F",
                                                                                                                                               "Aldh1a2+ MFB"="#FF7F00",
                                                                                                                                               "Prolif.MFB-1"="#FFFF99",
                                                                                                                                               "Prolif.MFB-2"="#FB9A99"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(title="Peri-portal mesenchymal cells")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Peri-portal mesenchymal cells.pdf", plot = P1, device="pdf", width=8, height=8)

P2<-plot_cell_trajectory(cds, cell_size = 1.5, show_tree = T,show_backbone = F,show_branch_points  = FALSE,color_by = "Pseudotime")+scale_colour_distiller(palette =  "Greys",direction = 1)+ 
  labs(title="Peri-portal mesenchymal cells")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "Peri-portal mesenchymal cells_Pseudotime.pdf", plot = P2, device="pdf", width=8, height=7)

P3<-plot_cell_trajectory(cds, cell_size = 1.5, show_tree = T, show_backbone = F,show_branch_points  = FALSE,color_by = "treatment")+scale_colour_manual(values=c("#737373","#2171B5","#238B45"))+ 
  labs(title="Peri-portal mesenchymal cells")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "Peri-portal mesenchymal cells_treatment.pdf", plot = P3, device="pdf", width=8, height=7)


library(viridis)

P4<-plot_cell_trajectory(cds, color_by = "cell_type",show_branch_points = F)+scale_colour_manual(values=c("Fib-1" ="#A6CEE3",
                                                                                "Fib-2" ="#1F78B4",
                                                                                "Fib-3" ="#B2DF8A",
                                                                                "Fib-4" = "#33A02C",
                                                                                "Fib-5" = "navy",
                                                                                "HSCs" = "#6A3D9A",
                                                                                "Robo2+ MFB" = "#CAB2D6",
                                                                                "Col6a5+ MFB"="#E31A1C",
                                                                                "Meox1+ MFB"="#FDBF6F",
                                                                                "Aldh1a2+ MFB"="#FF7F00",
                                                                                "Prolif.MFB-1"="#FFFF99",
                                                                                "Prolif.MFB-2"="#FB9A99"))+
  facet_wrap(~cell_type, nrow = 2)
ggsave(filename = "Peri-portal mesenchymal_splited_celltype.pdf", plot = P4, device="pdf", width=15, height=8)

#Cell density along pseudotime--celltype
library("ggridges")
plotdf=pData(cds)
P6<-ggplot(plotdf, aes(x=Pseudotime,y=cell_type,fill=cell_type,))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(7,22.5),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = c("#E31A1C",
                                 "#FDBF6F",
                                 "#A6CEE3", 
                                 "#FF7F00", 
                                 "#33A02C",
                                 "navy",
                                 "#1F78B4", 
                                 "#6A3D9A",
                                 "#CAB2D6",
                                 "#B2DF8A",
                                 "#FFFF99",
                                 "#FB9A99"))+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.y = element_text(color = 'black'))+
  theme(axis.ticks = element_line(color = 'black')) 

#Cell density along pseudotime----treatment
plotdf=pData(cds)
P7<-ggplot(plotdf, aes(x=Pseudotime,y=cell_type,fill=treatment,))+
  geom_density_ridges(scale=1, alpha=0.5) +
  geom_vline(xintercept = c(7,22.5),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(    panel.grid = element_blank()
  )+scale_fill_manual(values = c("#737373","#2171B5","#238B45"))+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.y = element_text(color = 'black'))+
  theme(axis.ticks = element_line(color = 'black')) 
ggsave(filename = "cell_type_山脊图.pdf", plot = P7, device="pdf", width=8, height=4)
topN <- c("Timp1", "Acta2", "Col1a1", "Cxcl12", "Clec3b", "Klf4")
topN_cds <- cds[topN,]
plot_genes_in_pseudotime(topN_cds, color_by = "Pseudotime")

P5<-plot_genes_in_pseudotime(topN_cds,ncol = 3, color_by = "cell_type", cell_size = 0.5)+scale_colour_manual(values=c("Fib-1" ="#A6CEE3",
                                                                                        "Fib-2" ="#1F78B4",
                                                                                        "Fib-3" ="#B2DF8A",
                                                                                        "Fib-4" = "#33A02C",
                                                                                        "Fib-5" = "navy",
                                                                                        "HSCs" = "#6A3D9A",
                                                                                        "Robo2+ MFB" = "#CAB2D6",
                                                                                        "Col6a5+ MFB"="#E31A1C",
                                                                                        "Meox1+ MFB"="#FDBF6F",
                                                                                        "Aldh1a2+ MFB"="#FF7F00",
                                                                                        "Prolif.MFB-1"="#FFFF99",
                                                                                        "Prolif.MFB-2"="#FB9A99"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "Timp1_Tagln_Klf4_Acta2_pseudotime0910.pdf", plot = P5, device="pdf", width=12, height=8)


library(monocle)
diff_test_res <- differentialGeneTest(cds[expressed_genes[1:20000],],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
pht <- plot_pseudotime_heatmap(cds[sig_gene_names,],
                               num_clusters = 6,
                               cores = 1,
                               show_rownames = T,
                               return_heatmap=T)
pht
pht$tree_row

clusters <- cutree(pht$tree_row, k = 6)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "GeneClusters"
table(clustering)
head(clustering)

library(magrittr)
library(tidyverse)
topNGenes <- top_n(diff_test_res, n = 600, desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

pht <- plot_pseudotime_heatmap2(
  cds[topNGenes,],
  num_clusters = 4,
  show_rownames = T
)
pht
library(ClusterGVis)
visCluster(object = pht,plot.type = "line")
modu<-pht$wide.res$gene
cluster<-pht$long.res$cluster
modu<-as.data.frame(modu)
colnames(modu)[1] <- 'gene'
modu$gene
pdf(file = "Mesenchymal_monocle_heatmap.pdf",height = 7,width = 6,onefile = F)
visCluster(object = pht,plot.type = "heatmap",pseudotime_col = c("#FFFFFF","#000000"), markGenes = c("Col3a1","Smoc2" ,"Bmp4","Ctgf", "Col5a2", "Pid1","Ackr3",
                                                                                                     "Klf4","Acta2","Pi16","Clec3b","Atf3","Uck2",
                                                                                                     "Il1r1","Il1b","Entpd2","Slc25a5","Vegfd",
                                                                                                     "Timp1","Gucy1a1","Gucy1b1","Ntrk1","Stmn1",
                                                                                                     "Jun","Ccl19","Tgfbr3","Cxcl1","Cxcl2","Cxcl10",
                                                                                                     "Cxcl13","Col1a2","Cxcl12","Wnt11","Ccl2","Atox1",
                                                                                                     "Eif5a","Runx1","Dpt","Clec3b","Pi16"))
dev.off()





#GESA analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)
library(ReactomePA)
library(enrichplot)
data <- read.csv("Mesenchymal-control_vs_BDL_new.csv",sep = ",",row.names = 1,  na.strings = "NA",stringsAsFactors = FALSE)
GF=data[,c(1,4)]#extract gene symbel and logFC value
#not nessisery, genePMP=transform(genePMP,logFC=-1*logFC)#reverse the negative value to positive value of logFC
names(GF) <- c( "SYMBOL","logFC")
GFID <- bitr(GF[,1], fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")#revert symbolID as entrezID and ensemble and creat a new datafram
geneList<-dplyr::inner_join(x=GFID,y=GF,by="SYMBOL")
a=geneList[!duplicated(geneList$ENTREZID), ]
b=a[,c(2,4)]

#start to creat own geneList
#feature 1: numeric vector
geneListready <- b[,2]
## feature 2: named vector
names(geneListready) <- as.character(b[,1])
## feature 3: decreasing order
geneListready <- sort(geneListready, decreasing = TRUE)
head(geneListready)
#Pathway Enrichment Analysis
y <- gsePathway(geneListready, nPerm=10000,
                pvalueCutoff=1,
                pAdjustMethod="BH", verbose=FALSE,organism = "mouse")
C=y@result
write.table(y@result, file="Reactome_gseresult_BDL_control.txt", sep="\t", row.names=T, quote=FALSE)

which(C == "R-MMU-9020702")
which(C == "R-MMU-5668541")
which(C == "R-MMU-380108")
which(C == "R-MMU-2559582")
which(C == "R-MMU-5358346")
which(C == "R-MMU-2132295")
which(C == "R-MMU-1236975")
which(C == "R-MMU-1500931")
which(C == "R-MMU-1483257")
which(C == "R-MMU-983712")
which(C == "R-MMU-382551")
which(C == "R-MMU-425407")
which(C == "R-MMU-196854")
which(C == "R-MMU-8957322")
gseaplot<-gseaplot2(y,c(11,15,14,74,113,118,162),base_size = 7.5,
                    rel_heights = c(1, 0.4, 0.3),
                    subplots = 1:2,
                    pvalue_table = FALSE,
                    color = brewer.pal(n=7,"Dark2"),
                    ES_geom = "line")
ggsave(filename = "GSEA_Mesenchymal_BDL_vs_control_new.pdf", plot = gseaplot, device="pdf", width=7, height=5)

gseadata=read.table("BDL_control_Gsea_selected.TXT",sep = "\t", header =T,quote = "", na.strings = "NA",stringsAsFactors = FALSE,comment.char ="")
library(ggpubr)
CB=ggdotchart(gseadata, x = "Description", y = "NES",
              color = "group",
              palette = c("#000000","#000000"), 
              sorting = "none", 
              add = "segments", 
              rotate = TRUE,
              add.params = list(color = "lightgray", size = 2), 
              group = "group",  
              ggtheme = theme_pubr() 
)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")


CB=CB + geom_point(aes(size=-1*log10(pvalue),color=group),shape=21,alpha=0.8,fill=c("#737373","#737373",
                                                                                    "#737373","#737373",
                                                                                    "#737373","#737373",
                                                                                    "#737373","#737373",
                                                                                    "#737373","#238B45",
                                                                                    "#238B45","#238B45",
                                                                                    "#238B45","#238B45",
                                                                                    "#238B45","#238B45",
                                                                                    "#238B45","#238B45",
                                                                                    "#238B45","#238B45"))+scale_size(range=c(4,8))+
theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
CB

ggsave(filename = "GSEA_Mesenchymal_BDL_vs_control_0902.pdf", plot = CB, device="pdf", width=6, height=6)

P2<-FeaturePlot(object = Mesenchymal, features ='Klf4',reduction = "umap",order = F,max.cutoff = 3.5, blend = F, pt.size = 1, combine = T,ncol = 1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
P2<-P2+theme(plot.title = element_text(face="italic",size=18, color = 'black'))+NoAxes()

P3<-FeaturePlot(object = Mesenchymal, features ='Acta2',reduction = "umap",order = F,max.cutoff = 3, blend = F, pt.size = 1, combine = T,ncol = 1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
P3<-P3+theme(plot.title = element_text(face="italic",size=18, color = 'black'))+NoAxes()


P5<-FeaturePlot(object = Mesenchymal, features ='Mfap4',reduction = "umap",order = F,max.cutoff = 3.5, blend = F, pt.size = 1, combine = T,ncol = 1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
P5<-P5+theme(plot.title = element_text(face="italic",size=18, color = 'black'))+NoAxes()

#Define the fib/MFBs at different stages of activation by Acta2 expression level
Acta2_mid=colnames(subset(x = Mesenchymal2, subset = Acta2 >=2&Acta2<10, slot = 'counts'))

Acta2_highORActa2_mid=ifelse(colnames(subset(Mesenchymal2)) %in% Acta2_high,'Acta2-high',
                             ifelse(colnames(subset(Mesenchymal2)) %in% Acta2_mid,'Acta2-mid','Acta2-low'))
table(Acta2_highORActa2_mid)
Mesenchymal2@meta.data$Acta2_highORActa2_mid=Acta2_highORActa2_mid
Mesenchymal2$Acta2_highORActa2_mid=factor(Mesenchymal2$Acta2_highORActa2_mid, levels=c('Acta2-low','Acta2-mid','Acta2-high'))

p8=DimPlot(Mesenchymal2, reduction = "umap", label = FALSE, split.by = "treatment", repel = TRUE,pt.size =1,group.by = 'Acta2_highORActa2_mid')+ 
  scale_colour_manual( values = c('lightblue',"purple","red"))

p8=DimPlot(Mesenchymal2, reduction = "umap", label = FALSE, repel = TRUE,pt.size =1,group.by = 'Acta2_highORActa2_mid')+ 
  scale_colour_manual( values = c('lightblue',"purple","red"))
ggsave(filename = "UMAP_acta2low_high0910.pdf", plot = p8, device="pdf", width=5.5, height=4)

P4<-VlnPlot(object = Mesenchymal2, features = c("Klf4"), group.by= 'Acta2_highORActa2_mid',cols =c('lightblue',"purple","red"))+ 
  scale_colour_manual(values = c('blue',"purple","red"))+ 
  geom_boxplot(width=.2,col="black",fill="white")+
  theme(plot.title = element_text(face="italic",color = 'black'))+
  NoLegend()
P4
ggsave(filename = "Klf4_acta2low_high.pdf", plot = P4, device="pdf", width=4, height=4)
