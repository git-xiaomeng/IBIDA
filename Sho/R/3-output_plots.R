# This script annotates clusters and apply annotations to markers from step.2

library(tidyverse)
library(magrittr)
library(Seurat)
library(RColorBrewer)

hcc <- readRDS("processed/hcc_cluster.RDS")
feature_genes <- read.csv("input/feature_genes.csv", encoding = "utf-8", stringsAsFactors = FALSE)
markers <- read.csv("processed/markers.csv", encoding = "utf-8", stringsAsFactors = FALSE)

# DimPlot
png("plots/DimPlot_celltype.png", width = 800*2.5, height = 600*3, res = 300)
DimPlot(hcc, group.by = "celltype", label = TRUE) + NoLegend() 
dev.off()

png("plots/DimPlot_celltype_sub.png", width = 800*3, height = 600*3.5, res = 200)
DimPlot(hcc, group.by = "celltype_sub", label = TRUE) + NoLegend() 
dev.off()

png("plots/DimPlot_tissue.png", width = 800*3, height = 600*3, res = 300)
DimPlot(hcc, group.by = "tissue") 
dev.off()

png("plots/DimPlot_splitby_tissue.png", width = 800*3, height = 600*2.5, res = 300)
DimPlot(hcc, split.by = "tissue", group.by = "celltype", ncol = 3) + NoAxes()
dev.off()

# FeaturePlot
p <- FeaturePlot(hcc, features = c("CD3D", "NKG7", "CD79A", "MS4A1", "HLA-DRA", "CD74",
                                   "CD68", "KIT", "CTSG", "SELL", "TCF4", "MKI67"),
                 pt.size = 0.01, min.cutoff = 0, combine = FALSE)
p <- lapply(X = p,
            FUN = function(p){
              return(p + NoLegend() + NoAxes() + theme(plot.title = element_text(size = 10)))
            })
png("plots/FeaturePlot.png", width = 800*2, height = 600*2, res = 300)
cowplot::plot_grid(plotlist = p, ncol = 4)
dev.off()


p <- FeaturePlot(hcc, features = c("LYZ", "LAMP3", "SLC40A1", "FCER1A"),
                 pt.size = 0.01, min.cutoff = 0, combine = FALSE)
p <- lapply(X = p,
            FUN = function(p){
              return(p + NoLegend() + NoAxes() + theme(plot.title = element_text(size = 8)))
            })
png("plots/FeaturePlot_targetgenes.png", width = 800*2.5, height = 600*1, res = 300)
cowplot::plot_grid(plotlist = p, ncol = 4)
dev.off()


png("plots/FeaturePlot_tissue_LAMP3.png", width = 800*5, height = 600*1.5, res = 200)
FeaturePlot(hcc, features = "LAMP3", split.by = "tissue", pt.size = 1, min.cutoff = 0)
dev.off()

png("plots/FeaturePlot_tissue_SLC40A1.png", width = 800*5, height = 600*1.5, res = 200)
FeaturePlot(hcc, features = "SLC40A1", split.by = "tissue", pt.size = 1, min.cutoff = 0)
dev.off()


# DotPlot
png("plots/Dotplot_CellSub.png", width = 800*4.5, height = 600*2.5, res = 150)
DotPlot(hcc, features = unique(feature_genes$gene), group.by = "celltype_sub", assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

png("plots/Dotplot_CellSub_targetgenes.png", width = 800*2.5, height = 600*1, res = 150)
DotPlot(hcc, features = c("LYZ", "LAMP3", "SLC40A1", "FCER1A"), group.by = "celltype_sub", assay = "RNA") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

# Propotion of cell subset in each tissue
prop <- prop.table(table(hcc$celltype, hcc$tissue), margin = 2) %>% as.data.frame.table

t <- 0.3
w <- 0.7
color <- hue_pal()

png("plots/frequency_Celltype.png", height = 1200, width = 1500, res = 200)
print(ggplot() + theme_classic() +
        geom_bar(data = prop,
                 aes(x = Var2, y = Freq, fill = Var1),
                 colour = "gray", size = t, width = w, stat="identity", position = "stack") +
        geom_segment(data = tidyr::spread(prop, Var2, Freq),
                     colour = "gray", size = t,
                     aes(x = 1 + w/2,
                         xend = 2 - w/2,
                         y = cumsum(rev(Ascites)),
                         yend = cumsum(rev(Blood)))) +
        geom_segment(data = tidyr::spread(prop, Var2, Freq),
                     colour = "gray", size = t,
                     aes(x = 2 + w/2,
                         xend = 3 - w/2,
                         y = cumsum(rev(Blood)),
                         yend = cumsum(rev(Lymphnode)))) +
        geom_segment(data = tidyr::spread(prop, Var2, Freq),
                     colour = "gray", size = t,
                     aes(x = 3 + w/2,
                         xend = 4 - w/2,
                         y = cumsum(rev(Lymphnode)),
                         yend = cumsum(rev(Normal)))) +
        geom_segment(data = tidyr::spread(prop, Var2, Freq),
                     colour = "gray", size = t,
                     aes(x = 4 + w/2,
                         xend = 5 - w/2,
                         y = cumsum(rev(Normal)),
                         yend = cumsum(rev(Tumor)))) +
        xlab("Tissue") + ylab("Proportion of Cells") + labs(fill = "Subset") +
        theme(axis.title = element_text(size = 14), plot.title = element_text(size = 18)) +
        scale_fill_discrete() + ggtitle("Cell Subsets in each Tissue"))
dev.off()



