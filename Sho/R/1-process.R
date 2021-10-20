# IBIDA Group1 Sho Hanakawa
## The goal is learning general coding of scRNA-seq data. We use 10x genomics data from database
## Zhang Q et al. Cell 2019. PMID: 31675496
## GSE140228: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228

# Tweak file names to be read by Seurat
## GSE140228_UMI_counts_Droplet.mtx.gz > "matrix.mtx.gz"
## GSE140228_UMI_counts_Droplet_barcodes.tsv.gz > "barcodes.tsv.gz"
## GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz - 
### rename MΦ as Macro. then save as "cellinfo.csv"
## GSE140228_UMI_counts_Droplet_genes.tsv.gz -
### remove columns except 1:2 and 1st raw. then compress as "features.tsv.gz"

# This script sets up the Seurat Objects
# Standard pre-processing
library(Seurat)
library(ggplot2)
library(cowplot)

hcc.count <- Read10X(data.dir = "input/GSE140228_UMI_counts_Droplet/")
hcc <- CreateSeuratObject(counts = hcc.count, project = "GSE140228")
rm(hcc.count)

sample <- read.csv("input/GSE140228_UMI_counts_Droplet/cellinfo.csv", encoding = "utf-8", stringsAsFactors = FALSE)

hcc$donor <- sample$Donor
hcc$tissue <- sample$Tissue
hcc$tissue_sub <- sample$Tissue_sub
hcc$celltype <- sample$celltype_global
hcc$celltype_sub <- sample$celltype_sub
hcc$histology <- sample$Histology

hcc[["percent.mt"]] <- PercentageFeatureSet(hcc, pattern = "^MT-")


# QC metrics
p1 <- FeatureScatter(hcc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(hcc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2

hcc_filt <- subset(hcc, nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 80000)

# plot pre&post QC metrics from Seurat object
png("plots/qc_summary.png", width = 800*3, height = 600*3, res = 200)
seurat_pre <- hcc
seurat_post <- hcc_filt
p1 <- VlnPlot(seurat_pre, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0.01)
p2 <- FeatureScatter(seurat_pre, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.001)
p3 <- FeatureScatter(seurat_post, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size = 0.001)
p4 <- FeatureScatter(seurat_post, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.001)
plot_grid(p1,
          p2 + theme(legend.position = "none"),
          p3 + theme(legend.position = "none"),
          p4 + theme(legend.position = "none"),
          labels = c("Pre", NA, "Post", NA), nrow = 2) +
  draw_figure_label(label = "GSE140228_HCC", position = "top.right", fontface = "bold", size = 10)
dev.off()

hcc <- hcc_filt
rm(list=setdiff(ls(), "hcc"))

png("plots/qc_summary_by_conditions.png", width = 800*2, height = 600*2.5, res = 200)
p1 <- VlnPlot(hcc, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by = "donor", pt.size = 0.01)
p2 <- VlnPlot(hcc, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by = "tissue", pt.size = 0.01)
plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          labels = c("Donor", "Tissue"), ncol = 1) +
  draw_figure_label(label = "GSE140228_HCC", position = "top.right", fontface = "bold", size = 10)
dev.off()


# Normalizing the data and identification of variable feaures
hcc <- NormalizeData(hcc)
hcc <- FindVariableFeatures(hcc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(hcc), 10)

png("plots/VariableFeatures.png", width = 800*2, height = 600*2, res = 200)
p1 <- VariableFeaturePlot(hcc)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
p2
dev.off()


rm(list=setdiff(ls(), "hcc"))
saveRDS(hcc, file = "processed/hcc_lognorm.RDS")

# Scaling the data and linear dimensional reduction
n_dims <- 100
hcc <- ScaleData(hcc, verbose = TRUE)
hcc <- RunPCA(hcc, npcs = n_dims)


# NOTE: This process can take a long time for big dataset
hcc <- JackStraw(hcc, num.replicate = 100, dims = n_dims)
hcc <- ScoreJackStraw(hcc, dims = 1:n_dims)


# Visualizing PCA, JacStrawplot and Elbowplot
pdf("plots/PCA.pdf")

print(DimHeatmap(hcc, dims = 1:12, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 13:24, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 25:36, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 37:48, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 49:60, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 61:72, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 73:84, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 85:96, cells = 500, balanced = TRUE))

dev.off()

png("plots/JackStraw.png", width = 800*3, height = 600*2, res = 200)
JackStrawPlot(hcc, dims = 1:n_dims)
dev.off()

png("plots/ElbowPlot.png", width = 800*2, height = 600*2, res = 200)
ElbowPlot(hcc, ndims = n_dims)
dev.off()


saveRDS(hcc, file = "processed/hcc_pca.RDS")


rm(list=ls())
hcc <- readRDS(file = "processed/hcc_pca.RDS")

# Clustering and visualization
## Jackstraw plot indicates 58-60 PCs
n_dim <- 60
res <- 1.2

hcc <- FindNeighbors(hcc, reduction = "pca", dims = 1:n_dim)
hcc <- FindClusters(hcc, resolution = res)

hcc <- RunUMAP(hcc, assay = "features", dims = 1:n_dim)


png(paste0("plots/DimPlot_dim-", n_dim, "_res-", res, ".png"), width = 800*4, height = 600*3, res = 200)

p1 <- DimPlot(hcc, label = TRUE)
p2 <- DimPlot(hcc, group.by = "donor")
p3 <- DimPlot(hcc, group.by = "tissue")
top_row <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1.1, 1, 1))

p4 <- DimPlot(hcc, group.by = "celltype_sub")
p5 <- DimPlot(hcc, group.by = "celltype")
bottom_row <- plot_grid(p4, p5, ncol = 2, rel_widths = c(1.3, 1))

plot_grid(top_row, bottom_row, ncol = 1) +
  draw_figure_label(label = paste0("GSE140228_HCC_dim-", n_dim, "_res-", res),
                    position = "top.right", fontface = "bold", size = 10)

dev.off()

png(paste0("plots/DimPlot_CellSub_dim-", n_dim, "_res-", res, ".png"), width = 800*3, height = 600*2, res = 100)
p1 <- p1 + NoLegend()
p4 <- DimPlot(hcc, group.by = "celltype_sub", label = TRUE) + NoLegend()
plot_grid(p1, p4, ncol = 2) +
  draw_figure_label(label = paste0("GSE140228_HCC_dim-", n_dim, "_res-", res),
                    position = "top.right", fontface = "bold", size = 10)
dev.off()


rm(list=setdiff(ls(), "hcc"))


saveRDS(hcc, file = "processed/hcc_cluster.RDS")

