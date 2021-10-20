# Caution: Run this script with high spec PC
# This script identifies markers
library(Seurat)
library(tidyverse)
library(future)

plan("multisession", workers = 16)
options(future.globals.maxSize = 4000 * 1024 ** 2)

hcc <- readRDS(file = "processed/hcc_cluster.RDS")

markers <- FindAllMarkers(hcc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")

write.csv(markers, file = "processed/markers.csv")

markers <- read.csv("processed/markers.csv", encoding = "utf-8", stringsAsFactors = FALSE)

# ClusterTree and Heatmap plotting top5 genes from each clusters

hcc <- BuildClusterTree(hcc, reorder = FALSE, reorder.numeric = FALSE)

png("plots/ClusterTree.png", width = 800*4, height = 600*4, res = 200)
PlotClusterTree(hcc)
dev.off()

top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

png("plots/Heatmap.png", width = 800*5, height = 600*3, res = 200)
DoHeatmap(hcc, features = unique(top5$gene)) +
        theme(axis.text.y = element_text(size = 5)) +
        NoLegend()
dev.off()


