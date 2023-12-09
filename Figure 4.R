#Cholangiocytes analysis_0601
#Set work dictionary and load workspace.
setwd("I:/Postdoc/Cholangiocytes")
Cholangiocytes <- readRDS("~/0213_Cholangiocytes_clean_Har2_sublustering.rds")
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(DT)
library(tidyverse)
library(patchwork)
library(harmony)
library(scRNAtoolVis)
library(ggunchull)
set.seed(123)
#Cell annotation
p2 <- DimPlot(Cholangiocytes, reduction = "umap", label = TRUE, repel = TRUE,group.by = 'RNA_snn_res.0.4',pt.size = 0.05)
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
obj.markers <- FindAllMarkers(Cholangiocytes, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
write.table(obj.markers, file="obj.markers_new.txt", sep="\t", row.names=FALSE, quote=FALSE)
plan("sequential")

top10 = obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, file="Top10DEGs_new.txt", sep="\t", row.names=FALSE, quote=FALSE)
#check the DEGs genes expression
p3<-FeaturePlot(object =Cholangiocytes, features = c("Epcam","Sox9","Pecam1","Alb","Col1a1","Acta2","Mki67","Ccl2","Il1b","Tnf","Tgfb1","Gli1"),reduction = "umap",order =F, blend = F, pt.size = 0.05, combine = T,ncol = 4)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p4<-FeaturePlot(object =Cholangiocytes, features = c("Pigr","Sox9","Pecam1","Alb","Col1a1","Acta2","Mki67","Ccl2","Il1b","Tnf","Tgfb1","Gli1"),reduction = "umap",order =F, blend = F, pt.size = 0.05, combine = T,ncol = 4)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p5<-FeaturePlot(object =Cholangiocytes, features = c("Pigr","Tff2"),split.by ='treatment',   reduction = "umap",order =F, blend = F, pt.size = 0.05, combine = T,ncol = 4)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
levels(Cholangiocytes)
#Rename the subclusters
new.cluster.ids <- c("chol-1", 
                     "chol-2", 
                     "chol-3",
                     "Prolif.chol-1", 
                     "Prolif.chol-2")
names(new.cluster.ids) <- levels(Cholangiocytes)
Cholangiocytes<- RenameIdents(Cholangiocytes, new.cluster.ids)

DimPlot(Cholangiocytes, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)

#Reorder cluster
Cholangiocytes@active.ident <- factor(Cholangiocytes@active.ident, 
                                      levels=c("chol-1", 
                                               "chol-2", 
                                               "chol-3",
                                               "Prolif.chol-1", 
                                               "Prolif.chol-2"
                                      ))
Cholangiocytes$cell_type<-Cholangiocytes@active.ident
#---------------------------------------------------------------#
Cholangiocytes <- readRDS("~/Cholangiocytes.rds")
VlnPlot(Cholangiocytes, features = "Ccl2", pt.size =0, split.by = 'treatment',group.by = "cell_type",
        stack = F,flip = T,cols=c("#737373","#2171B5","#238B45") )+labs(x="",y="Gene expression")+ theme(axis.title.y.right = element_text(face="italic",color = "#000000"))+
  theme(axis.text.y =element_text(color = "#000000",size=10))+theme(axis.text.x = element_text(angle = 90))+geom_vline(xintercept=c(1.5,2.5,3.5,4.5), colour='black',linetype='dotted',size=0.5)+
  theme(plot.tag = element_text(face="italic",color = 'black'))
#UMAP
colours=c("#006837", "#c2e699","#66c2a4","#a6611a","#dfc27d") 
p1<-DimPlot(Cholangiocytes, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1)+scale_colour_manual(values =colours )+
  guides(colour = guide_legend(override.aes = list(size=4)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))
ggsave(filename = "Umap_w_extra_Chol.pdf", plot = p1, device="pdf", width=7.5, height=5.4)
Cholangiocytes<-subset(x = Cholangiocytes, idents = c("chol-1", 
                                                               "chol-2",
                                                               "chol-3",
                                                               "Prolif.chol-1", 
                                                               "Prolif.chol-2"))
p1<-DimPlot(Cholangiocytes, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1)+scale_colour_manual(values =colours )+
  guides(colour = guide_legend(override.aes = list(size=4)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))
ggsave(filename = "Umap_Chol.pdf", plot = p1, device="pdf", width=7.5, height=5.4)
#Cell cycle analysis
p11<-DimPlot(Cholangiocytes, reduction = "umap",group.by = "Phase", label = F, repel = TRUE,pt.size = 1)+scale_colour_manual(values =c("#a6611a","red","blue"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))
ggsave(filename = "cellcycle_Chol.pdf", plot = p11, device="pdf", width=7, height=5.4)

table(Cholangiocytes$orig.ident)
My_levels <- c('Control_1','Control_2','Control_3','CDAHFD_1','CDAHFD_2','BDL_1','BDL_2')
Cholangiocytes$orig.ident<-factor(Cholangiocytes$orig.ident, levels= My_levels)

prop.table(table(Cholangiocytes$Phase))
table(Cholangiocytes$Phase,Cholangiocytes$orig.ident)
Cellratio<-prop.table(table(Cholangiocytes$Phase,Cholangiocytes$orig.ident),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1,levels=c("S",
                                                 "G2M",
                                                 "G1"))

table(Cellratio$Var1)
colourCount=length(unique(Cellratio$Var1))
cellcycle<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  theme_classic()+labs(x="",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values = c( "blue","red","#a6611a"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = "cellcycle_Freq_vertical_new.pdf", plot = cellcycle, device="pdf", width=5, height=5)

#Find DEGs of subcluster once again
library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
obj.markers <- FindAllMarkers(Cholangiocytes, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
write.table(obj.markers, file="obj.markers_new.txt", sep="\t", row.names=FALSE, quote=FALSE)
plan("sequential")
top10 = obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, file="Top10DEGs_new.txt", sep="\t", row.names=FALSE, quote=FALSE)
#Heatmap
colours=c("#006837", "#c2e699","#66c2a4","#a6611a","#dfc27d") 
cols.use <- list(cell_type=colours,
                 treatment=c("#737373","#2171B5","#238B45"))
#Add metadata
Cholangiocytes@meta.data$treatment<- factor(x = Cholangiocytes@meta.data$treatment, levels=c('Control', 'CDAHFD', 'BDL'))
Cholangiocytes$cell_type<-Cholangiocytes@active.ident
table(Cholangiocytes$treatment)
My_levels <- c('Control','CDAHFD','BDL')
Cholangiocytes$treatment<-factor(Cholangiocytes$treatment, levels= My_levels)
#DoMultiBarHeatmap function
library(grid)
Heatmap<-DoMultiBarHeatmap(Cholangiocytes, features = top10$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme(axis.text.y = element_text(face="italic",color = 'black'))+
  guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(filename = "Cholangiocytes_heatmap.pdf", plot = Heatmap, device="pdf", width=6, height=6.6)

#Composition proportion analysis
table(Cholangiocytes$treatment)
My_levels <- c('Control','CDAHFD','BDL')
Cholangiocytes$treatment<-factor(Cholangiocytes$treatment, levels= My_levels)

prop.table(table(Idents(Cholangiocytes)))
table(Idents(Cholangiocytes),Cholangiocytes$treatment)
Cellratio<-prop.table(table(Idents(Cholangiocytes),Cholangiocytes$treatment),margin=2)
options(scipen=200)
Cellratio
Cellratio<-as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1,levels=c("chol-1", 
                                                 "chol-2", 
                                                 "chol-3",
                                                 "Prolif.chol-1", 
                                                 "Prolif.chol-2"
))

table(Cellratio$Var1)
colourCount=length(unique(Cellratio$Var1))
cells_in_treatment<-ggplot(Cellratio)+geom_bar(aes(x=Var2,y=Freq*100,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour=NA)+
  theme_classic()+labs(x="Sample",y="% cells in sample")+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+
  scale_fill_manual(values = c("#006837", "#c2e699","#66c2a4","#a6611a","#dfc27d"))
ggsave(filename = "Cholangiocytes_cluster_Freq_vertical_new.pdf", plot = cells_in_treatment, device="pdf", width=6, height=7.8)

##highlight the cholangiocytes in each condition, respectively.
Cholangiocytes2<-SetIdent(Cholangiocytes,value=Cholangiocytes@meta.data$treatment)
hi1<-WhichCells(Cholangiocytes2,idents = "Control")
hi2<-WhichCells(Cholangiocytes2,idents = "CDAHFD")
hi3<-WhichCells(Cholangiocytes2,idents = "BDL")
sham_UMAP<-DimPlot(Cholangiocytes,group.by = "treatment", reduction = "umap",cells.highlight= hi1, cols.highlight = c("#737373"),sizes.highlight = 0.5, label = F,repel = F,  pt.size = 1,cols = c("#BDBDBD"))
ggsave(filename = "shamUMAP_Cholangiocytes_distrbution_new.pdf", plot = sham_UMAP, device="pdf", width=6.2, height=4.5)

CDAHFD_UMAP<-DimPlot(Cholangiocytes,group.by = "treatment", reduction = "umap",cells.highlight= hi2, cols.highlight = c("#2171B5"),sizes.highlight = 0.5, label = F,repel = F,  pt.size = 1,cols = c("#BDBDBD"))
ggsave(filename = "CDAHFD_UMAP_Cholangiocytes_distrbution_new.pdf", plot = CDAHFD_UMAP, device="pdf", width=6.2, height=4.5)

BDL_UMAP<-DimPlot(Cholangiocytes,group.by = "treatment", reduction = "umap",cells.highlight= hi3, cols.highlight = c("#238B45"),sizes.highlight = 0.5, label = F,repel = F,  pt.size = 1,cols = c("#BDBDBD"))
ggsave(filename = "BDL_UMAP_Cholangiocytes_distrbution_new.pdf", plot = BDL_UMAP, device="pdf", width=6.2, height=4.5)
library(ggpubr)
P4<-ggarrange(sham_UMAP,CDAHFD_UMAP,BDL_UMAP,ncol = 3,nrow = 1, legend = "top")
ggsave(filename = "Cholangiocytes_distrbution_UMAP_new2.pdf", plot = P4, device="pdf", width=12, height=4.5)

#Proportions of each cell cluster in each treatment condition
b<-as.data.frame(table(Idents(Cholangiocytes),Cholangiocytes$orig.ident)) 
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

#chol-1
p3<-ggplot(data=CellRa[CellRa$Var1=="chol-1",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,40)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="chol-1") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
# Portal vein ECs
p4<-ggplot(data=CellRa[CellRa$Var1=="chol-2",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,15)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="chol-2") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Lymphatic ECs
p5<-ggplot(data=CellRa[CellRa$Var1=="chol-3",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,10)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="chol-3") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Arterial ECs
p6<-ggplot(data=CellRa[CellRa$Var1=="Prolif.chol-1",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Prolif.chol-1") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Prolif.ECs
p7<-ggplot(data=CellRa[CellRa$Var1=="Prolif.chol-2",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,0.6)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Prolif.chol-2") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

library(ggpubr)
p10<-ggarrange(p3,p4,p5,p6,p7,ncol = 3,nrow = 2, legend = "top",common.legend = T)

ggsave(filename = "Cholangiocytes_each_Freq_Total_new.pdf", plot = p10, device="pdf", width=11.25, height=8)

#Go function analysis
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(tidyverse)
library(viridis)
markers<-read.table("obj.markers_new.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
View(markers)
#select top DEGs for Go function analysis
sigposDEG.all<- subset(markers, p_val_adj<0.01&abs(avg_log2FC)>1) 

Chol1_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='chol-1') 
Chol2_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='chol-2') 
Chol3_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='chol-3') 
Prolif.Chol1_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Prolif.chol-1')
Prolif.Chol2_up <- subset(sigposDEG.all, avg_log2FC>1&cluster=='Prolif.chol-2') 

list<-list(Chol1_up,Chol2_up,Chol3_up,Prolif.Chol1_up,Prolif.Chol2_up)
names(list)[1:5] <- c("chol-1", 
                      "chol-2", 
                      "chol-3",
                      "Prolif.chol-1", 
                      "Prolif.chol-2")

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

#read the function analysis results
Chol1 <- read.csv("GO_new_chol-1.CSV", row.names = 1)
Chol2 <- read.csv("GO_new_chol-2.CSV", row.names = 1)
Chol3 <- read.csv("GO_new_chol-3.CSV", row.names = 1)
Prolif.Chol1 <- read.csv("GO_new_Prolif.chol-1.CSV",row.names = 1)
Prolif.Chol2<- read.csv("GO_new_Prolif.chol-2.CSV",row.names = 1)


Chol1$group<-"chol-1"
Chol2$group<-"chol-2"
Chol3$group<-"chol-3"
Prolif.Chol1$group<-"Prolif.chol-1"
Prolif.Chol2$group<-"Prolif.chol-2"

Chol1<-Chol1[1:20,]
Chol2<-Chol2[1:20,]
Chol3<-Chol3[1:15,]
Prolif.Chol1<-Prolif.Chol1[1:20,]
Prolif.Chol2<-Prolif.Chol2[1:20,]

all <- rbind(Chol1,Chol2,Chol3,Prolif.Chol1,Prolif.Chol2)
all$Description<- gsub("-"," ",all$Description)

library(forcats)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
dat2<-str_split(all$GeneRatio,"/",simplify = T)[,1]
all$Num<-as.numeric(dat2) 
dat3<-str_split(all$GeneRatio,"/",simplify = T)[,2]
all$Tot<-as.numeric(dat3) 
all$Generatio<-as.numeric(all$Num/all$Tot)

My_levels <- c("chol-1", 
               "chol-2", 
               "chol-3",
               "Prolif.chol-1", 
               "Prolif.chol-2")
all$group<-factor(all$group, levels= My_levels)
ggplot(all, aes(group, Description)) +
  geom_point(aes(color=qvalue, size=Generatio))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
selct<-list(
  "DNA replication",
  "chromosome condensation",
  "mitotic cell cycle phase transition",
  "positive regulation of cell cycle",
  "mitotic nuclear division",
  
  "apoptotic mitochondrial changes",
  "release of cytochrome c from mitochondria",
  "intrinsic apoptotic signaling pathway",
  "regulation of body fluid levels",
  "water homeostasis",
  
  "response to interferon gamma",
  "leukocyte migration",
  "neutrophil migration",
  "neutrophil chemotaxis",
  "granulocyte chemotaxis",
  "tumor necrosis factor production",
  "regulation of leukocyte migration",
  "cell chemotaxis",
  "acute inflammatory response",
  
  "reverse cholesterol transport",
  "negative regulation of cholesterol transport",
  "cholesterol transport",
  "phospholipid efflux",
  "receptor mediated endocytosis"
)

selctterm<-as.character(selct)

Chol1 <- read.csv("GO_new_chol-1.CSV", row.names = 1)
Chol2 <- read.csv("GO_new_chol-2.CSV", row.names = 1)
Chol3 <- read.csv("GO_new_chol-3.CSV", row.names = 1)
Prolif.Chol1 <- read.csv("GO_new_Prolif.chol-1.CSV",row.names = 1)
Prolif.Chol2<- read.csv("GO_new_Prolif.chol-2.CSV",row.names = 1)

Chol1$group<-"chol-1"
Chol2$group<-"chol-2"
Chol3$group<-"chol-3"
Prolif.Chol1$group<-"Prolif.chol-1"
Prolif.Chol2$group<-"Prolif.chol-2"

Chol1$Description<- gsub("-"," ",Chol1$Description)
Chol2$Description<- gsub("-"," ",Chol2$Description)
Chol3$Description<- gsub("-"," ",Chol3$Description)
Prolif.Chol1$Description<- gsub("-"," ",Prolif.Chol1$Description)
Prolif.Chol2$Description<- gsub("-"," ",Prolif.Chol2$Description)

Chol1=Chol1[Chol1$Description%in%selctterm,]
Chol2=Chol2[Chol2$Description%in%selctterm,]
Chol3=Chol3[Chol3$Description%in%selctterm,]
Prolif.Chol1=Prolif.Chol1[Prolif.Chol1$Description%in%selctterm,]
Prolif.Chol2=Prolif.Chol2[Prolif.Chol2$Description%in%selctterm,]

all <- rbind(Chol1,Chol2,Chol3,Prolif.Chol1,Prolif.Chol2)

library(forcats)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
dat2<-str_split(all$GeneRatio,"/",simplify = T)[,1]
all$Num<-as.numeric(dat2) 
dat3<-str_split(all$GeneRatio,"/",simplify = T)[,2]
all$Tot<-as.numeric(dat3) 
all$Generatio<-as.numeric(all$Num/all$Tot)

My_levels <- c("chol-1", 
               "chol-2", 
               "chol-3",
               "Prolif.chol-1", 
               "Prolif.chol-2")

all$group<-factor(all$group, levels= My_levels)
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

ggsave(filename = "cholangiocyte_GO_new.pdf", plot = MyeGO, device="pdf", width=5, height=4.2)
#######Correlation analysis among cholangiocytes subclusters##########
library (pheatmap)
Cholangiocytes$cell_type<-Cholangiocytes@active.ident

table(Cholangiocytes$cell_type)
av <-AverageExpression(Cholangiocytes,
                       group.by = "cell_type",
                       assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av, 1, sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = 'spearman'))
cor<-pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(15))
ggsave(filename = "cholangiocyte_correlation.pdf", plot = cor, device="pdf", width=6, height=5)
#Load the gene sets
Geneset <- readxl::read_xlsx("~/cholan_signature.xlsx")
#Convert into list
Profibrogenesis<-as.character(Geneset$profibrogenesis)
Profibrogenesis<-Profibrogenesis[1:19]

library(ggpubr)
boxplot<-ggboxplot(Cholangiocytes@meta.data, x = "cell_type", y = "Profibrogenesis1",
                   color = "black",fill = 'treatment',palette=c("#737373","#2171B5","#238B45"),
                   bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))+ labs(x = '')
ggsave(filename = "Profibrogenesis_score_boxplot.pdf", plot = boxplot, device="pdf", width=8.5, height=4.5)

#P value caculation
data<-Cholangiocytes@meta.data
library(rstatix)
wilcox_res<-data%>%
  group_by(cell_type) %>%
  wilcox_test(Profibrogenesis1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="Profibrogenesis_wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)

#Dot heatmap for function related marker genes
Cholangiocytes$celltype.treatment<- paste(Cholangiocytes$cell_type, Cholangiocytes$treatment, sep = "_")
My_levels <- c('chol-1_Control',
               'chol-1_CDAHFD',
               'chol-1_BDL',
               'chol-2_Control',
               'chol-2_CDAHFD',
               'chol-2_BDL',
               'chol-3_Control',
               'chol-3_CDAHFD',
               'chol-3_BDL',
               'Prolif.chol-1_Control',
               'Prolif.chol-1_CDAHFD',
               'Prolif.chol-1_BDL',
               "Prolif.chol-2_Control",
               "Prolif.chol-2_CDAHFD",
               "Prolif.chol-2_BDL"
)

Cholangiocytes$celltype.treatment<-factor(Cholangiocytes$celltype.treatment, levels= My_levels)
levels(Cholangiocytes)
Idents(Cholangiocytes) <- "celltype.treatment"
cholan_geneset <- readxl::read_xlsx("~/cholan_signature.xlsx")

bind<-as.character(cholan_geneset)

plot<-DotPlot(Cholangiocytes,group= "celltype.treatment",cols = c("#1E3163", "#00C1D4", "#FFED99","#FF7600"),   
              features =bind , dot.scale = 8)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()

plot<-plot+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+theme(axis.text.x = element_text(angle = 90))+
  scale_y_discrete(labels = c('','chol-1','','','chol-2','',
                              '','chol-3','',"",'Prolif.chol-1','',
                              '','Prolif.chol-2',''))+
  theme(axis.text.x = element_text(angle=90,
                                   vjust=1, 
                                   size=11,
                                   hjust=1,
                                   color = 'black')) +
  theme(axis.text.y = element_text(face="italic",
                                   color = 'black'))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.justification=c(1,2))+
  geom_hline(yintercept=c(3.5,6.5,9.5,12.5), linetype='dotted', col = 'black',size=0.5)+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
ggsave(filename = "cholangiocytes_dotheatmap.pdf", plot = plot, device="pdf", width=5, height=8)

#Interaction between cholangiocytes and myeloid cells
library(harmony)
library(CellChat)
library(ComplexHeatmap)
library(ggpubr)
set.seed(123)
Cholangiocytes <- readRDS("~/Cholangiocytes.rds")
colours=c("#006837", "#c2e699","#66c2a4","#a6611a","#dfc27d") 
p1<-DimPlot(Cholangiocytes, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1)+scale_colour_manual(values =colours )+
  guides(colour = guide_legend(override.aes = list(size=4)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

Myeloid <- readRDS("~/Myeloid.rds")
colours2=c("#258DD2", "#46806F","#0764A2","#6CC5FF","#4BB394","#8AC15A","#BAF388","#7570B3","#054169",
           "#D95F02","#8FF7D8")
DimPlot(Myeloid, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)+ scale_colour_manual(values = colours2)+NoLegend()
Myeloid$cell_type<-Myeloid@active.ident
table(Myeloid$cell_type)

scMerge<- merge(Myeloid, y=c(Cholangiocytes))
plot<-VlnPlot(scMerge, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
scMerge_QC<- NormalizeData(scMerge) %>% FindVariableFeatures() 
?CaseMatch
library(stringr)
s.genes.m<-str_to_title(cc.genes$s.genes)
g2m.genes.m<-str_to_title(cc.genes$g2m.genes)
CaseMatch(c(s.genes.m,g2m.genes.m),VariableFeatures(scMerge_QC))
g2m_genes=CaseMatch(search = g2m.genes.m,match = rownames(scMerge_QC))
s_genes=CaseMatch(search = s.genes.m,match = rownames(scMerge_QC))
scMerge_QC<-CellCycleScoring(object = scMerge_QC,g2m.features = g2m_genes,s.features = s_genes)
scMerge_QC<-ScaleData(scMerge_QC)
scRNAa=RunPCA(scMerge_QC,features = c(s_genes,g2m_genes))
p<-DimPlot(scRNAa,reduction = "pca",group.by = "Phase")

library(future)
options(future.globals.maxSize=800000000000)
plan("multisession",workers=4)
scMerge_QC<-ScaleData(scMerge_QC,features = rownames(scMerge_QC))
plan("sequential")
scMerge_QC<- RunPCA(scMerge_QC,features = VariableFeatures(scMerge_QC),verbose=FALSE)

scMerge_Har2<- scMerge_QC%>% RunHarmony("orig.ident", plot_convergence = TRUE, theta = 1)
rm(scMerge_QC)
scMerge_Har2<- scMerge_Har2 %>% RunUMAP(reduction = "harmony", dims = 1:40) %>% FindNeighbors(reduction = "harmony", dims = 1:40) %>% FindClusters(resolution = seq(from=0.1, to=0.4, by=0.1)) %>% identity()
p2 <- DimPlot(scMerge_Har2, reduction = "umap", label = TRUE, repel = TRUE,group.by = 'RNA_snn_res.0.3',pt.size = 0.05)

scMerge_Har2$cell_type<-factor(scMerge_Har2$cell_type)
scMerge_Har2$treatment<-factor(scMerge_Har2$treatment)
scMerge_Har2$orig.ident<-factor(scMerge_Har2$orig.ident)

table(scMerge_Har2$cell_type)
scMerge_Har2<-SetIdent(scMerge_Har2,value=scMerge_Har2@meta.data$cell_type)

colours=c("#006837", "#c2e699","#66c2a4","#a6611a","#dfc27d","#258DD2", "#6CC5FF","#0764A2","#054169","#8AC15A","#BAF388","#4BB394","#46806F","#D95F02","#8FF7D8","#7570B3") 
scMerge_Har2@active.ident <- factor(scMerge_Har2@active.ident, 
                                    levels=c("Chol-1", 
                                             "Chol-2", 
                                             "Chol-3",
                                             "Prolif.Chol-1", 
                                             "Prolif.Chol-2",
                                             "RM1s", 
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
p1<-DimPlot(scMerge_Har2, reduction = "umap", label = F, repel = TRUE,pt.size = 0.3)+scale_colour_manual( values = colours)+
  guides(colour = guide_legend(override.aes = list(size=4)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))

setwd("I:/Postdoc/Cholangiocytes/interaction")
ggsave(filename = "integration_cholangiocytes_myeloid.pdf", plot =p1, device="pdf", width=8, height=4.5)

#############################Cell chat##############################################
control<- scMerge_Har2[,scMerge_Har2@meta.data$orig.ident %in% c("control_1","control_2","control_3")]
control$cell_type = Idents(control)
table(control$cell_type)

control_datainput <- control[["RNA"]]@data
control_meta <- control@meta.data
control_meta <- control_meta[,c(1,26)]
colnames(control_meta) <- c("group","labels")


control_cellchat <- createCellChat(object = control_datainput, 
                                meta = control_meta, 
                                group.by = "labels")
control_cellchat <- addMeta(control_cellchat, meta = control_meta)
control_cellchat <- setIdent(control_cellchat, ident.use = "labels") 
levels(control_cellchat@idents)
groupSize <- as.numeric(table(control_cellchat@idents)) 
CellChatDB <- CellChatDB.mouse

CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) 
control_cellchat@DB <- CellChatDB.use
control_cellchat <- subsetData(control_cellchat, features = NULL)
control_cellchat <- identifyOverExpressedGenes(control_cellchat)
control_cellchat <- identifyOverExpressedInteractions(control_cellchat)
control_cellchat <- projectData(control_cellchat, PPI.mouse)
#--------------computeCommunProb-----
control_cellchat <- computeCommunProb(control_cellchat, raw.use=T)
control_cellchat <- filterCommunication(control_cellchat, min.cells = 10)
control.net <- subsetCommunication(control_cellchat)  
write.csv(control.net, file = "control_net_inter_all.csv")

control_cellchat <- computeCommunProbPathway(control_cellchat)
control_cellchat <- aggregateNet(control_cellchat)
save(control_cellchat,control_datainput, control_meta, control.net, groupSize, file = "control_cellchat.RData")

par(mfrow = c(1,2))
netVisual_circle(control_cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")


netVisual_circle(control_cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

mat <- control_cellchat@net$weight
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
  
}

pathway.show <- control.net$pathway_name
vertex.receiver = seq(1,10)

netVisual_aggregate(control_cellchat, 
                    signaling = pathway.show[1],
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy") 


#circle plot
netVisual_aggregate(control_cellchat, 
                    signaling = pathway.show[1], 
                    layout = "circle")
par(mfrow = c(10,10), xpd=TRUE)
netVisual_aggregate(control_cellchat, 
                    signaling = pathway.show[1], 
                    layout = "chord")



netVisual_heatmap(control_cellchat, signaling =  pathway.show[1], color.heatmap = "Reds")

netVisual_bubble(control_cellchat, 
                 sources.use = "Cholangiocytes",
                 remove.isolate = FALSE,
                 font.size=14)

#CDAHFD
CDAHFD<- scMerge_Har2[,scMerge_Har2@meta.data$orig.ident %in% c("CDAHFD_1","CDAHFD_2")]

table(Idents(CDAHFD))
CDAHFD$cell_type = Idents(CDAHFD)
table(CDAHFD$cell_type)
CDAHFD_datainput <- CDAHFD[["RNA"]]@data
CDAHFD_meta <- CDAHFD@meta.data
CDAHFD_meta <- CDAHFD_meta[,c(1,26)]
colnames(CDAHFD_meta) <- c("group","labels")

CDAHFD_cellchat <- createCellChat(object = CDAHFD_datainput, 
                                  meta = CDAHFD_meta, 
                                  group.by = "labels")
CDAHFD_cellchat <- addMeta(CDAHFD_cellchat, meta = CDAHFD_meta)
CDAHFD_cellchat <- setIdent(CDAHFD_cellchat, ident.use = "labels")
levels(CDAHFD_cellchat@idents)
groupSize <- as.numeric(table(CDAHFD_cellchat@idents)) 

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
CDAHFD_cellchat@DB <- CellChatDB.use
CDAHFD_cellchat <- subsetData(CDAHFD_cellchat, features = NULL)
CDAHFD_cellchat <- identifyOverExpressedGenes(CDAHFD_cellchat)
CDAHFD_cellchat <- identifyOverExpressedInteractions(CDAHFD_cellchat)
CDAHFD_cellchat <- projectData(CDAHFD_cellchat, PPI.mouse)
#--------------computeCommunProb-----
CDAHFD_cellchat <- computeCommunProb(CDAHFD_cellchat, raw.use=T)
CDAHFD_cellchat <- filterCommunication(CDAHFD_cellchat, min.cells = 15)
CDAHFD.net <- subsetCommunication(CDAHFD_cellchat) 
write.csv(CDAHFD.net, file = "CDAHFD_net_inter.csv")
CDAHFD_cellchat <- computeCommunProbPathway(CDAHFD_cellchat)
CDAHFD_cellchat <- aggregateNet(CDAHFD_cellchat)
save(CDAHFD_cellchat,CDAHFD_datainput, CDAHFD_meta, CDAHFD.net, groupSize, file = "CDAHFD_cellchat.RData")


par(mfrow = c(1,2))
netVisual_circle(CDAHFD_cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(CDAHFD_cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

mat <- CDAHFD_cellchat@net$weight
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

pathway.show <- CDAHFD.net$pathway_name
vertex.receiver = seq(1,10)
netVisual_aggregate(CDAHFD_cellchat, 
                    signaling = pathway.show[1],
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy") 


#circle plot
netVisual_aggregate(CDAHFD_cellchat, 
                    signaling = pathway.show[1], 
                    layout = "circle")
par(mfrow = c(10,10), xpd=TRUE)
netVisual_aggregate(CDAHFD_cellchat, 
                    signaling = pathway.show[1], 
                    layout = "chord")

netVisual_heatmap(CDAHFD_cellchat, signaling =  pathway.show[1], color.heatmap = "Reds")

BDL<- scMerge_Har2[,scMerge_Har2@meta.data$orig.ident %in% c("BDL_1","BDL_2")]
BDL<-subset(BDL,idents=unique(Idents(BDL))[!unique(Idents(BDL))=="Mesothelial cells"])
BDL$cell_type = Idents(BDL)
table(Idents(BDL))
table(BDL$cell_type)
BDL_datainput <- BDL[["RNA"]]@data
BDL_meta <- BDL@meta.data
BDL_meta <- BDL_meta[,c(1,26)]
colnames(BDL_meta) <- c("group","labels")
BDL_cellchat <- createCellChat(object = BDL_datainput, 
                               meta = BDL_meta, 
                               group.by = "labels")
BDL_cellchat <- addMeta(BDL_cellchat, meta = BDL_meta)
BDL_cellchat <- setIdent(BDL_cellchat, ident.use = "labels") 
levels(BDL_cellchat@idents)
groupSize <- as.numeric(table(BDL_cellchat@idents)) 
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
BDL_cellchat@DB <- CellChatDB.use
BDL_cellchat <- subsetData(BDL_cellchat, features = NULL)
BDL_cellchat <- identifyOverExpressedGenes(BDL_cellchat)
BDL_cellchat <- identifyOverExpressedInteractions(BDL_cellchat)
BDL_cellchat <- projectData(BDL_cellchat, PPI.mouse)

#--------------computeCommunProb-----
BDL_cellchat <- computeCommunProb(BDL_cellchat, raw.use=T)
BDL_cellchat <- filterCommunication(BDL_cellchat, min.cells = 10)
BDL.net <- subsetCommunication(BDL_cellchat) 
write.csv(BDL.net, file = "BDL_net_inter.csv")
BDL_cellchat <- computeCommunProbPathway(BDL_cellchat)
BDL_cellchat <- aggregateNet(BDL_cellchat)
save(BDL_cellchat,BDL_datainput, BDL_meta, BDL.net, groupSize, file = "BDL_cellchat.RData")

par(mfrow = c(1,2))
netVisual_circle(BDL_cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(BDL_cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
mat <- BDL_cellchat@net$weight
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
pathway.show <- BDL.net$pathway_name
vertex.receiver = seq(1,10)
netVisual_aggregate(BDL_cellchat, 
                    signaling = pathway.show[1],
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy") 
#circle plot
netVisual_aggregate(BDL_cellchat, 
                    signaling = pathway.show[1], 
                    layout = "circle")
par(mfrow = c(10,10), xpd=TRUE)
netVisual_aggregate(BDL_cellchat, 
                    signaling = pathway.show[1], 
                    layout = "chord")

netVisual_heatmap(BDL_cellchat, signaling =  pathway.show[1], color.heatmap = "Reds")
netVisual_bubble(BDL_cellchat, 
                 sources.use = "Cholangiocytes",
                 remove.isolate = FALSE,
                 font.size=14)

#cell-cell interaction comparative analysis

cellchat.list <- list(control=control_cellchat,
                      CDAHFD=CDAHFD_cellchat,
                      BDL=BDL_cellchat
)
control_cellchat <- netAnalysis_computeCentrality(control_cellchat)
CDAHFD_cellchat <- netAnalysis_computeCentrality(CDAHFD_cellchat)
BDL_cellchat <- netAnalysis_computeCentrality(BDL_cellchat)
cellchat <- mergeCellChat(cellchat.list,
                          add.names = names(cellchat.list))
library(CellChat)
par(mfrow = c(1,2))
p1 <- compareInteractions(cellchat, show.legend = F,color.use =c("#737373","#2171B5","#238B45"), group = c(1,2,3),width = 0.6,size.text=12)
p2 <- compareInteractions(cellchat, show.legend = F, color.use =c("#737373","#2171B5","#238B45"), group = c(1,2,3),width = 0.6, measure = "weight",size.text=12)
num_int_stren<-p1+p2
ggsave(filename = "num_int_stren.pdf", plot =num_int_stren, device="pdf", width=5, height=3)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat,color.use = colours, weight.scale = T,edge.width.max = 5,comparison = c(1,2))
netVisual_diffInteraction(cellchat, color.use = colours, weight.scale = T, edge.width.max = 5,measure = "weight",comparison = c(1,2))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, color.use = colours,weight.scale = T,edge.width.max = 5,comparison = c(1,3))
netVisual_diffInteraction(cellchat, weight.scale = T,  color.use = colours,measure = "weight",edge.width.max = 5,comparison = c(1,3))

gg1 <- netVisual_heatmap(cellchat, color.use=colours,
                         comparison = c(1,2), color.heatmap = c("#11355a","#b2182b"))


gg2 <- netVisual_heatmap(cellchat, color.use=colours,measure = "weight",
                         comparison = c(1,2),color.heatmap = c("#11355a","#b2182b"))
gg1 + gg2

weight.max <- getMaxWeight(cellchat.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(cellchat.list)) {
  netVisual_circle(cellchat.list[[i]]@net$count, 
                   weight.scale = T, label.edge= F,
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 8, color.use = colours,
                   title.name = paste0("Number of interactions - ", names(cellchat.list)[i]))
}
cellchat.list = lapply(cellchat.list, function(x){x = netAnalysis_computeCentrality(x)})
num.link <- sapply(cellchat.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # group2 the dot size in the different datasets
gg <- list()
for (i in 1:length(cellchat.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat.list[[i]],color.use = colours, title = names(cellchat.list)[i], weight.MinMax = weight.MinMax)
}
pdot=patchwork::wrap_plots(plots = gg)&ylim(0,30)&xlim(0,65)
ggsave(filename = "dot_inc_out2.pdf", plot = pdot, device="pdf", width=12, height=3.5)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
rankSimilarity(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

pathway.show=c("ANNEXIN","CD40","NGR","TRAIL", "IL6","BAFF", "IL1",
               "MHC-II","NCAM","WNT","FGF","SEMA5","NOTCH",
               "BMP","BMP10","GDF","MK","PTN","SPP1","XCR",
               "BRADYKININ","TNF","ESAM","PARs","VWF","PVR") 
gg1 <- rankNet(cellchat, mode = "comparison", cutoff.pvalue = 0.01,slot.name='netP',
               tol = 0.05,thresh = 0.01,stacked = T,do.flip = F,color.use = c("#737373","#2171B5"), do.stat = T,comparison = c(1,2)) + 
  theme(legend.position = "top")+theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))
gg11 <- rankNet(cellchat, mode = "comparison",  cutoff.pvalue = 0.01,
                tol = 0.05,color.use = c("#737373","#238B45"),do.flip = F,slot.name='netP', 
                thresh = 0.01, stacked = T, do.stat = T,comparison = c(1,3))+ 
  theme(legend.position = "top")+theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))
gg1 <- rankNet(cellchat, mode = "comparison", cutoff.pvalue = 1,slot.name='netP',
               tol = 0.05,thresh = 0.05,stacked = T,do.flip = F,color.use = c("#737373","#2171B5"), do.stat = T,comparison = c(1,2)) + 
  theme(legend.position = "top")+theme(axis.text.y = element_text( color = "black"))

gg11 <- rankNet(cellchat, mode = "comparison",  cutoff.pvalue = 1,
                tol = 0.05,color.use = c("#737373","#238B45"),do.flip = F,slot.name='netP', 
                thresh = 0.05, stacked = T, do.stat = T,comparison = c(1,3))+ 
  theme(legend.position = "top")+theme(axis.text.y = element_text( color = "black"))

library(ComplexHeatmap)
library(CellChat)
i = 1
pathway.union <- union(cellchat.list[[i]]@netP$pathways, cellchat.list[[i+1]]@netP$pathways)
pathway.union <- union(pathway.union, cellchat.list[[i+2]]@netP$pathways)
#outgoing
ht1 = netAnalysis_signalingRole_heatmap(cellchat.list[[i]], color.use =colours,pattern = "outgoing", signaling = pathway.union, title = names(cellchat.list)[i], width = 4, height = 20,color.heatmap = "Reds")
ht2 = netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]],color.use =colours, pattern = "outgoing", signaling = pathway.union, title = names(cellchat.list)[i+1], width = 4, height = 20,color.heatmap = "Reds")
ht3 = netAnalysis_signalingRole_heatmap(cellchat.list[[i+2]],color.use =colours, pattern = "outgoing", signaling = pathway.union, title = names(cellchat.list)[i+2], width = 4, height =20,color.heatmap = "Reds")
draw(ht1 + ht2+ ht3, ht_gap = unit(0.5, "cm"))
#incoming
ht4 = netAnalysis_signalingRole_heatmap(cellchat.list[[i]],color.use =colours, pattern = "incoming", signaling = pathway.union, title = names(cellchat.list)[i], width = 4, height = 25, color.heatmap = "Reds")
ht5 = netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]], color.use =colours,pattern = "incoming", signaling = pathway.union, title = names(cellchat.list)[i+1], width = 4, height = 25, color.heatmap = "Reds")
ht6 = netAnalysis_signalingRole_heatmap(cellchat.list[[i+2]], color.use =colours,pattern = "incoming", signaling = pathway.union, title = names(cellchat.list)[i+2], width = 4, height = 25, color.heatmap = "Reds")

draw(ht1 + ht2+ ht3+ht4 + ht5+ht6, ht_gap = unit(0.5, "cm"))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "BD-LAMs",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#2171B5"),comparison = c(1, 2))
ggsave(filename = "BD-LAMs_control_CDAHFD.pdf", plot = gg1, device="pdf", width=5, height=3)

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "BD-LAMs",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#238B45"),comparison = c(1, 3),signaling.exclude = "COLLAGEN")
ggsave(filename = "BD-LAMs_control_BDL.pdf", plot = gg2, device="pdf", width=5, height=3)

gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Chol-1",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#238B45"),comparison = c(1, 3),signaling.exclude = "COLLAGEN")
ggsave(filename = "Chol-1_control_BDL.pdf", plot = gg3, device="pdf", width=5, height=3)

gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Chol-2",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#238B45"),comparison = c(1, 3),signaling.exclude = "COLLAGEN")
ggsave(filename = "Chol-2_control_BDL.pdf", plot = gg4, device="pdf", width=5, height=3)

gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Chol-3",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#238B45"),comparison = c(1, 3),signaling.exclude = "COLLAGEN")
ggsave(filename = "Chol-3_control_BDL.pdf", plot = gg5, device="pdf", width=5, height=3)

gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Prolif.Chol-1",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#238B45"),comparison = c(1, 3),signaling.exclude = "COLLAGEN")
ggsave(filename = "Prolif.Chol-1_control_BDL.pdf", plot = gg6, device="pdf", width=5, height=3)

gg7 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Prolif.Chol-2",  dot.alpha = 1,dot.size = 2, color.use = c("black", "#737373","#238B45"),comparison = c(1, 3),signaling.exclude = "COLLAGEN")
ggsave(filename = "Prolif.Chol-2_control_BDL.pdf", plot = gg7, device="pdf", width=5, height=3)

#cholangiocytes as targets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_aggregate(cellchat.list[[1]],color.use = colours, targets.use =c("Chol-1","Chol-2","Chol-3","Prolif.Chol-1","Prolif.Chol-2") ,thresh = 0.001, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[1]))
netVisual_aggregate(cellchat.list[[3]],color.use = colours, targets.use =c("Chol-1","Chol-2","Chol-3","Prolif.Chol-1","Prolif.Chol-2") ,thresh = 0.001, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[3]))
#cholangiocytes as sources
par(mfrow = c(1,2), xpd=TRUE)
netVisual_aggregate(cellchat.list[[1]],color.use = colours, sources.use =c("Chol-1","Chol-2","Chol-3","Prolif.Chol-1","Prolif.Chol-2") ,thresh = 0.001, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[1]))
netVisual_aggregate(cellchat.list[[3]],color.use = colours, sources.use =c("Chol-1","Chol-2","Chol-3","Prolif.Chol-1","Prolif.Chol-2") ,thresh = 0.001, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[3]))


cellchat.list <- list(control=control_cellchat,BDL=BDL_cellchat
)
control_cellchat <- netAnalysis_computeCentrality(control_cellchat)
#CDAHFD_cellchat <- netAnalysis_computeCentrality(CDAHFD_cellchat)
BDL_cellchat <- netAnalysis_computeCentrality(BDL_cellchat)

cellchat_control_BDL <- mergeCellChat(cellchat.list,
                                   add.names = names(cellchat.list))
pos.dataset = "BDL"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat2 <- identifyOverExpressedGenes(cellchat_control_BDL, 
                                        group.dataset = "datasets", 
                                        pos.dataset = pos.dataset,
                                        features.name = features.name,
                                        only.pos = FALSE,
                                        thresh.pc = 0.1,
                                        thresh.fc = 0.1,
                                        thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat2, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat2, net = net,
                              datasets = "BDL",ligand.logFC = 0.05,
                              receptor.logFC = 0.6)

# Chord diagram
tmp.pairLR.use <- list()
tmp.pairLR.use[["Cell-Cell Contact"]] <- data.frame( 
  interaction_name=CellChatDB$interaction$interaction_name[CellChatDB$interaction$annotation=="Cell-Cell Contact"]
)
tmp.pairLR.use[["Secreted Signaling"]] <- data.frame( 
  interaction_name=CellChatDB$interaction$interaction_name[CellChatDB$interaction$annotation=="Secreted Signaling"]
)
tmp.pairLR.use[["ECM-Receptor"]] <- data.frame( 
  interaction_name=CellChatDB$interaction$interaction_name[CellChatDB$interaction$annotation=="ECM-Receptor"]
)


netVisual_chord_gene(cellchat.list[[2]],
                     targets.use = c("Chol-1","Chol-2","Chol-3","Prolif.Chol-1","Prolif.Chol-2"),
                     slot.name = 'net',net = net.up,color.use = colours, pairLR.use = tmp.pairLR.use$`Secreted Signaling`,
                     lab.cex = 0.8,small.gap = 0.5, thresh = 0.05,
                     title.name = paste0("Up-regulated signaling in ",names(cellchat.list)[2]))
ggsave(filename = "Up-regulated_Secreted_sinaling_BDL.pdf", plot = gg8, device="pdf", width=8, height=4)

net.up2 <- subsetCommunication(cellchat2, net = net,
                               datasets = "BDL",ligand.logFC = 0.4,
                               receptor.logFC = NULL)

netVisual_chord_gene(cellchat.list[[2]],
                     sources.use = c("Chol-1","Chol-2","Chol-3","Prolif.Chol-1","Prolif.Chol-2"),
                     targets.use = c("BD-LAMs"),pairLR.use = tmp.pairLR.use$`Secreted Signaling`,
                     slot.name = 'net',net = net.up2, thresh = 0.05,
                     big.gap = 15,color.use = colours,
                     lab.cex = 0.8,small.gap = 3,
                     title.name = paste0("Up-regulated signaling in ",names(cellchat.list)[2]))

cellchat.list <- list(control=control_cellchat,CDAHFD=CDAHFD_cellchat)
control_cellchat <- netAnalysis_computeCentrality(control_cellchat)
CDAHFD_cellchat <- netAnalysis_computeCentrality(CDAHFD_cellchat)
#BDL_cellchat <- netAnalysis_computeCentrality(BDL_cellchat)

#cellchat.list = lapply(cellchat.list, function(x){x = netAnalysis_computeCentrality(x)})

cellchat_control_CDAHFD <- mergeCellChat(cellchat.list,
                                      add.names = names(cellchat.list))
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "CDAHFD"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat3 <- identifyOverExpressedGenes(cellchat_control_CDAHFD, 
                                        group.dataset = "datasets", 
                                        pos.dataset = pos.dataset,
                                        features.name = features.name,
                                        only.pos = FALSE,
                                        thresh.pc = 0.1,
                                        thresh.fc = 0.1,
                                        thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat3, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS

net.up3 <- subsetCommunication(cellchat3, net = net,
                               datasets = "CDAHFD",ligand.logFC = NULL,
                               receptor.logFC = 0.2)
netVisual_chord_gene(cellchat.list[[2]],
                     sources.use = c("Chol-1","Chol-2","Chol-3"),
                     targets.use = c("BD-LAMs"),pairLR.use = tmp.pairLR.use$`Secreted Signaling`,
                     slot.name = 'net',net = net.up3, thresh = 0.05,
                     big.gap = 15,color.use = colours,
                     lab.cex = 0.8,small.gap = 3,
                     title.name = paste0("Up-regulated signaling in ",names(cellchat.list)[2]))


net.up4<- subsetCommunication(cellchat3, net = net,
                              datasets = "CDAHFD",ligand.logFC = 0.05,
                              receptor.logFC = NULL)
netVisual_chord_gene(cellchat.list[[2]],
                     targets.use = c("Chol-1","Chol-2","Chol-3"),
                     sources.use = c("BD-LAMs"),
                     pairLR.use = tmp.pairLR.use$`Secreted Signaling`,
                     slot.name = 'net',net = net.up4, thresh = 0.05,
                     big.gap = 18,color.use = colours,
                     lab.cex = 0.8,small.gap = 4,
                     title.name = paste0("Up-regulated signaling in ",names(cellchat.list)[2]))

