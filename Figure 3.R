---
title: "Figure 3"
author: "Lin LEI"
---
setwd("~/Figure3")
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
T_NK<-subset(integration, idents =c("T & NK cells"))
p2 <- DimPlot(T_NK, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.05)

a<-ls()
rm(list = a[which(a!='integration')])

scMerge_QC<- NormalizeData(T_NK) %>% FindVariableFeatures() 
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
T_NK_QC<-ScaleData(scMerge_QC,vars.to.regress=c("S.Score","G2M.Score"),features = rownames(scMerge_QC))
##Explicitly close multisession workers, if they were used
plan("sequential")

T_NK_QC<- RunPCA(T_NK_QC,features = VariableFeatures(T_NK_QC),verbose=FALSE)

ElbowPlot(T_NK_QC)
pc.num=1:20
#Harmony re-integrating
T_NK_Har2<- T_NK_QC%>% RunHarmony("orig.ident", plot_convergence = TRUE, theta = 1)
rm(T_NK_QC)
T_NK<- T_NK_Har2 %>% RunUMAP(reduction = "harmony", dims = 1:40) %>% FindNeighbors(reduction = "harmony", dims = 1:40) %>% FindClusters(resolution = seq(from=0.1, to=0.4, by=0.1)) %>% identity()
#check the data with cell-type specific markers
##Featureplot
p1=FeaturePlot(object = T_NK, features = "Cd3e",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p1=p1+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p2=FeaturePlot(object = T_NK, features = "Cd8a",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p2=p2+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p3=FeaturePlot(object = T_NK, features = "Cd4",min.cutoff =0, max.cutoff = 2,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p3=p3+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p4=FeaturePlot(object = T_NK, features = "Tmem176a",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p4=p4+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p5=FeaturePlot(object = T_NK, features = "Foxp3",min.cutoff =0, max.cutoff = 3.5,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p5=p5+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p6=FeaturePlot(object = T_NK, features = "Scart1",min.cutoff =0, max.cutoff = 3,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p6=p6+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p7=FeaturePlot(object = T_NK, features = "Lyst",min.cutoff =0, max.cutoff = 3,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p7=p7+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p8=FeaturePlot(object = T_NK, features = "Ncr1",min.cutoff =0, max.cutoff = 3,reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
p8=p8+theme(plot.title = element_text(face="italic",color = 'black'))+NoAxes()

p9<-ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow = 2,common.legend = T, legend = "right")+theme(plot.title = element_text(face="italic",color = 'black',size=16))
ggsave(filename = "totalmaker_UMAP.pdf", plot = p9, device="pdf", width=12, height=6)

#Define the cell clusters
new.cluster.ids <- c("CD8+ T cells", 
                     "CD4+ T cells", 
                     "NKT cells", 
                     "NK cells",
                     "Th17 cells", 
                     "γδ T cells",
                     "TRegs"
)
names(new.cluster.ids) <- levels(T_NK)
T_NK <- RenameIdents(T_NK, new.cluster.ids)
levels(T_NK)
#reorder cluster
T_NK@active.ident <- factor(T_NK@active.ident, 
                            levels=c("CD8+ T cells",
                                     "CD4+ T cells",
                                     "Th17 cells",
                                     "TRegs", 
                                     "γδ T cells", 
                                     "NKT cells",
                                     "NK cells"
                            ))
#UMAP plot
colours=c("#E41A1C", "#377EB8","#395c69","#395eea", "#FFFF33","#4DAF4A","#984EA3") 
p1 <- DimPlot(T_NK, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.5) +
  scale_colour_manual(values = colours) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
ggsave(filename = "Umap0812.pdf", plot = p1, device="pdf", width=6, height=4)

#DEGs clusters
T_obj.markers <- FindAllMarkers(T_NK, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
write.table(T_obj.markers, file="T_obj.markers.txt", sep="\t", row.names=FALSE, quote=FALSE)

T_NKobj <-read.table("~/T_obj.markers.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)


#heatmap
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
colours=c("#E41A1C", "#377EB8","#395c69",  "#395eea", "#FFFF33", "#4DAF4A","#984EA3") 
cols.use <- list(cell_type=colours,
                 treatment=c("#737373","#2171B5","#238B45"))
T_NK@meta.data$treatment<- factor(x = T_NK@meta.data$treatment, levels=c('Control', 'CDAHFD', 'BDL'))

T_NK$cell_type<-T_NK@active.ident
table(T_NK$treatment)
My_levels <- c('Control','CDAHFD','BDL')
T_NK$treatment<-factor(T_NK$treatment, levels= My_levels)
library(grid)
Heatmap<-DoMultiBarHeatmap(T_NK, features = top10$gene, group.by = "cell_type", additional.group.by = "treatment",additional.group.sort.by = c('treatment'),cols.use = cols.use)+ 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename = "T_NK_heatmap_addanno2.pdf", plot = Heatmap, device="pdf", width=6, height=6)

#Proportions of each cell cluster in each treatment condition
b<-as.data.frame(table(Idents(T_NK),T_NK$orig.ident)) 
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
            file="Cellratio_Tot0809.txt", sep="\t",row.names=FALSE, quote=FALSE)

CellRa<-read.table("Cellratio_Tot0809.txt" ,sep = "\t", header = T, na.strings = "NA",stringsAsFactors = FALSE)
CellRa$group <- factor(CellRa$group,levels=c("Control","CDAHFD","BDL"))

#CD8+ T cells
p3<-ggplot(data=CellRa[CellRa$Var1=="CD8+ T cells",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1.8)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="CD8+ T cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
# CD4+ T cells
p4<-ggplot(data=CellRa[CellRa$Var1=="CD4+ T cells",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1.5)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="CD4+ T cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Th17 cells
p5<-ggplot(data=CellRa[CellRa$Var1=="Th17 cells",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="Th17 cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Th17 cells
p6<-ggplot(data=CellRa[CellRa$Var1=="TRegs",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,0.6)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="TRegs") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#Th17 cells
p7<-ggplot(data=CellRa[CellRa$Var1=="γδ T cells",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="γδ T cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#NK cells
p8<-ggplot(data=CellRa[CellRa$Var1=="NK cells",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="NK cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")
#NKT cells
p9<-ggplot(data=CellRa[CellRa$Var1=="NKT cells",], aes(group,Freq2*100))+
  geom_jitter(aes(fill = group),position = position_jitter(0.02),shape=21, size = 6)+
  scale_fill_manual(values=c("#737373","#2171B5","#238B45"))+theme_classic()+ylim(0,1)+
  theme( axis.line=element_line(colour="black",size=0.25),
         axis.title=element_text(size=13,face="plain",color="black"),
         axis.text = element_text(size=12,face="plain",color="black"),
         panel.border = element_blank(),
         legend.position="none")+labs(x="",y="")+labs(title="NKT cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_summary(fun.y=mean, geom="point", shape="—", size=6, color="black")

library(ggpubr)
p10<-ggarrange(p3,p4,p5,p6,p7,p8,p9,ncol = 4,nrow = 2, legend = "top",common.legend = T)
ggsave(filename = "T_NK_each_Freq_Total.pdf", plot = p10, device="pdf", width=16, height=6.5)


#dot heatmap
#Load the gene markers
T_NK_geneset <- readxl::read_xlsx("~/T_NK_signature.xlsx")
View(T_NK_geneset)
#convert into list 
Cytotoxicity<- as.character(T_NK_geneset$Cytotoxicity)[1:6]
Naive<- as.character(T_NK_geneset$Naive)[1:5]
Costimulation<- as.character(T_NK_geneset$Costimulation)[1:3]
suppression<- as.character(T_NK_geneset$suppression)[1:5]
antigen_presentation_processing<- as.character(T_NK_geneset$`antigen presentation and processing`)[1:6]

bind<-c(Naive,Cytotoxicity,Costimulation,Suppression,antigen_presentation_processing)
bind<-bind[!duplicated(bind)]
My_levels <- c('CD8+ T cells_Control',
               'CD8+ T cells_CDAHFD',
               'CD8+ T cells_BDL',
               'CD4+ T cells_Control',
               'CD4+ T cells_CDAHFD',
               'CD4+ T cells_BDL',
               'Th17 cells_Control','Th17 cells_CDAHFD','Th17 cells_BDL',
               'TRegs_Control',
               'TRegs_CDAHFD',
               'TRegs_BDL',
               'γδ T cells_Control',
               'γδ T cells_CDAHFD',
               'γδ T cells_BDL',
               'NK cells_Control',
               'NK cells_CDAHFD',
               'NK cells_BDL',
               "NKT cells_Control",
               "NKT cells_CDAHFD",
               "NKT cells_BDL"
)
T_NK$celltype.treatment<-factor(T_NK$celltype.treatment, levels= My_levels)

plot<-DotPlot(T_NK,group= "celltype.treatment",cols = c("#1E3163", "#00C1D4", "#FFED99","#FF7600"),   
              features =bind , dot.scale = 8)+coord_flip()

plot<-plot+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+theme(axis.text.x = element_text(angle = 90))+#+theme(legend.position = "top")
  scale_y_discrete(labels = c('','CD8+ T cells','','','CD4+ T cells','',
                              '','Th17 cells','',"",'TRegs','',
                              '','γδ T cells','','','NK cells','',
                              '','NKT cells',''))+
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
  geom_hline(yintercept=c(3.5,6.5,9.5,12.5,15.5,18.5,21.5), linetype='dotted', col = 'black',size=0.5)
ggsave(filename = "T_NK_genesets_dotheatmap.pdf", plot = plot, device="pdf", width=8, height=9)

#Residency markers Cd69 and Rgs1
p10<-FeaturePlot(object = T_NK, features = c("Cd69","Rgs1"),reduction = "umap",order = T,max.cutoff = 4,  blend = F, pt.size = 1, combine = T,ncol = 2)&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))&NoAxes()

library(ggpubr)
p10<-ggarrange(p10[[1]],p10[[2]],ncol = 2,nrow = 1,common.legend = TRUE,legend = "right")
ggsave(filename = "residence_marker_Umap.pdf", plot = p10, device="pdf", width=6, height=2.8)

##Cxcr6+ CD8T cell density analysis
library("Nebulosa")
Nash<- T_NK[,T_NK@meta.data$orig.ident %in% c("CDAHFD_1","CDAHFD_2")]
Control<- T_NK[,T_NK@meta.data$orig.ident %in% c("Control_1","Control_2","Control_3")]
BDL<- T_NK[,T_NK@meta.data$orig.ident %in% c("BDL_1","BDL_2")]

p1<-plot_density(Control, c("Cxcr6", "Cd8a"),reduction = "umap",  joint = TRUE,size = 1)&scale_color_distiller(limits=c(-0.0001,0.003),palette = "Spectral")
p2 <- plot_density(Nash, c("Cxcr6", "Cd8a"),reduction = "umap",  joint = TRUE,size = 1)&scale_color_distiller(limits=c(-0.0001,0.003),palette = "Spectral")
p3 <- plot_density(BDL, c("Cxcr6", "Cd8a"),reduction = "umap",joint = TRUE,size = 1)&scale_color_distiller(limits=c(-0.0001,0.003),palette = "Spectral")
p4<-ggarrange(p1[[3]],p2[[3]],p3[[3]],ncol = 3,nrow = 1,common.legend = TRUE, legend = "right")
ggsave(filename = "Cxcr6+Cd8a+ cells density.pdf", plot = p4, device="pdf", width=12, height=3.5)

#Gene signature calculation
T_NK<-AddModuleScore(object = T_NK,features=list(c("Prf1","Ctsw","Cst7", "Gzmb","Gzma","Ccl5","Klf2", "Nkg7","Gzmk","Fasl","Cd27","Cd5")),ctrl = 100,name='cytotoxicity')
FeaturePlot(object =T_NK, features = "cytotoxicity1",reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

T_NK<-AddModuleScore(object = T_NK,features=list(c("Lef1", "Tcf7", "Sell", "Ccr7", "Foxp1")),ctrl = 100,name='Naive')
FeaturePlot(object =T_NK, features = "Naive1",reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

T_NK<-AddModuleScore(object =T_NK,features=list(c("Cxcr6","Prf1","Stat3","Foxo1","Ikzf1","Gata3","Sh2d2a","Gzmk","Gzma","Gzmb","Pdcd1","Bhlhe40")),ctrl = 100,name='autoaggressive')
P12<-FeaturePlot(object =T_NK, features = "autoaggressive1",reduction = "umap",order = F, blend = F, pt.size = 1, combine = T)&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

#Boxplot_Cytotoxicity
boxplot<-ggboxplot(T_NK@meta.data, x = "cell_type", y = "cytotoxicity1",
                   color = "black",fill = 'treatment',palette=c("#737373","#2171B5","#238B45"),
                   bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))+ labs(x = '')+NoLegend()
ggsave(filename = "cytotoxicity_score_boxplot.pdf", plot = boxplot, device="pdf", width=7.8, height=3.5)

data<-T_NK@meta.data
library(rstatix)
wilcox_res<-data%>%
  group_by(cell_type) %>%
  wilcox_test(cytotoxicity1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="cytotoxicity_score_wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)

#Boxplot_autoaggressive score
boxplot<-ggboxplot(T_NK@meta.data, x = "cell_type", y = "autoaggressive1",
                   color = "black",fill = 'treatment',palette=c("#737373","#2171B5","#238B45"),
                   bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))+ labs(x = '')+NoLegend()
ggsave(filename = "autoaggressive_score_boxplot.pdf", plot = boxplot, device="pdf", width=8, height=3.5)

data<-T_NK@meta.data
library(rstatix)
wilcox_res<-data%>%
  group_by(cell_type) %>%
  wilcox_test(autoaggressive1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="autoaggressive_score_wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)

#Boxplot_Naive
boxplot<-ggboxplot(T_NK@meta.data, x = "cell_type", y = "Naive1",
                   color = "black",fill = 'treatment',palette=c("#737373","#2171B5","#238B45"),
                   bxp.errorbar = T,bxp.errorbar.width = 0.5,size = 0.5,outlier.shape = NA,legend="right")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))+ labs(x = '')+NoLegend()
ggsave(filename = "Naive_score_boxplot.pdf", plot = boxplot, device="pdf", width=7.8, height=3.5)

data<-T_NK@meta.data
library(rstatix)
wilcox_res<-data%>%
  group_by(cell_type) %>%
  wilcox_test(Naive1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

write.table(wilcox_res, 
            file="Naive_score_wilcox_res.txt", sep="\t",row.names=FALSE, quote=FALSE)

#Cxcr6+Cd8a+ and Cxcr6-Cd8a+_GSVA analysis
T_NK<-SetIdent(T_NK,value=T_NK@meta.data$cell_type)
CD8cells<-subset(x = T_NK, idents = c("CD8+ T cells"))
Othercells<-subset(x = T_NK, idents = c("CD4+ T cells",
                                        "Th17 cells",
                                        "TRegs", 
                                        "γδ T cells", 
                                        "NK cells",
                                        "NKT cells"))
DimPlot(CD8cells)
#Define Cxcr6+Cd8a+ positive cells
CD8cells$pos.neg = ifelse(CD8cells@assays$RNA@counts['Cxcr6',]>0,'pos','neg')
table(CD8cells$pos.neg)
P4<-DimPlot(CD8cells, group.by ="pos.neg",  reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1)+scale_colour_manual(values =c("#9ecae1","#E0367A"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+NoLegend()+scale_x_continuous(limits = c(-8, 10)) +
  scale_y_continuous(limits = c(-6, 6))

P5<-DimPlot(Othercells, group.by ="treatment",  reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1)+scale_colour_manual(values =c("lightgrey","lightgrey","lightgrey"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text( color = "black"),axis.text.y = element_text( color = "black"))+NoLegend()+scale_x_continuous(limits = c(-8, 10)) +
  scale_y_continuous(limits = c(-6, 6))
library(ggpubr)
P6<-ggarrange(P4,P5,ncol = 2,nrow = 1,align = "hv",common.legend = T,legend = "right")
ggsave(filename = "CXCR6+_CD8T_higlighted.pdf", plot = P6, device="pdf", width=8, height=3.5)
#GSVA 
BiocManager::install("GSEABase")
BiocManager::install("GSVA", version = "3.14") 
library('GSEABase')
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

CD8Tmat <- as.matrix(CD8cells@assays$RNA@counts)
meta <- CD8cells@meta.data[,c("orig.ident", "treatment", "cell_type", "pos.neg")]
library(msigdbr)
mouse <- msigdbr(species = "Mus musculus")
mouse[1:5,1:5]
table(mouse$gs_cat) 
table(mouse$gs_subcat)
mouse_GO_bp = msigdbr(species = "Mus musculus",
                      category = "C5", 
                      subcategory = "GO:BP") %>% 
  dplyr::select(gs_name,gene_symbol)
mouse_GO_bp_Set = mouse_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
mouse_REACTOME_cp = msigdbr(species = "Mus musculus",
                            category = "C2", 
                            subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name,gene_symbol)
mouse_REACTOME_cp_Set = mouse_REACTOME_cp %>% split(x = .$gene_symbol, f = .$gs_name)

T_gsva <- gsva(expr = CD8Tmat, 
               gset.idx.list = mouse_GO_bp_Set,
               kcdf="Poisson", 
               parallel.sz = 5)

T_gsva_REACTOME<- gsva(expr = CD8Tmat, 
                       gset.idx.list = mouse_REACTOME_cp_Set,
                       kcdf="Poisson", 
                       parallel.sz = 5)


write.table(T_gsva, 'Cxcr6pos_neg_CD8_gsva.xls', row.names=T, col.names=NA, sep="\t")

group <- CD8cells$pos.neg %>% as.factor()
desigN <- model.matrix(~ 0 + group) 
colnames(desigN) <- levels(group)
fit = lmFit(test_control, desigN)
fit2 <- eBayes(fit)
diff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
write.csv(diff, file = "Diff.csv")


library(limma)
group <- CD8cells@meta.data[,c("pos.neg")]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(T_gsva_REACTOME)
fit <- lmFit(T_gsva_REACTOME,design)
contrast.matrix<-makeContrasts("pos-neg",levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)
test_control_Rectome_gsva <- topTable(fit2,adjust='fdr',coef=1,number=Inf)

write.csv(test_control_Rectome_gsva, file = "Cxcr6pos_Cxcr6neg_Rectome_gsva.csv")
Diff<-test_control_Rectome_gsva 

## barplot
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t,
                       adjp=Diff$adj.P.Val)
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "REACTOME_","")
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
dat_plot$threshold =sapply(1:nrow(dat_plot),function(x){
  if(dat_plot[x,"t"]>2 & dat_plot[x,"adjp"]<0.05){return("up")}
  else if(dat_plot[x,"t"]<-2 & dat_plot[x,"adjp"]<0.05){return("down")}
  else{return("noSig")}
})
df=dat_plot

df <- df[which(abs(df$t)>3),]
dfup = top_n(df, 12, t)
dfdown = top_n(df, -6, t)
df=rbind(dfup,dfdown)
df$id<- gsub('REACTOME_','',df$id)
df$id <- tolower(df$id)
df$id <- gsub('_',' ',df$id)
df$id<-sub("(.)", "\\U\\1",df$id,perl=TRUE)

df$hjust = ifelse(df$t>0,1,0)
df$nudge_y = ifelse(df$t>0,-0.1,0.1)
sortdf <- df[order(df$t),]
sortdf$id <- factor(sortdf$id, levels = sortdf$id)
limt = max(abs(df$t))

biplot<-ggplot(sortdf, aes(id, t,fill=threshold)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","noSig","up"),
                    values = c("#9ecae1","grey","#E0367A"))+
  geom_text(data = df, aes(label = id, y = nudge_y),
            nudge_x =0,nudge_y =0,hjust =df$hjust,
            size = 3.5)+
  labs(x = 'Reactome cp',
       y="t value of GSVA score",
       title = '')+
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 12),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = limt)
ggsave(filename = "CXCR6_pos_neg_REACTOME_GSVA.pdf", plot = biplot, device="pdf", width=6, height=7)







