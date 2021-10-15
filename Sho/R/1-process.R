# IBIDA Group1 Sho Hanakawa
## The goal is learning general coding of scRNA-seq data. We use 10x genomics data from the published article below
## Zhang Q, He Y, Luo N, Patel SJ et al. Landscape and Dynamics of Single Immune Cells in Hepatocellular Carcinoma. Cell 2019 Oct 31;179(4):829-845.e20. PMID: 31675496
## GSE140228: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228

# Tweak file names to be read by Seurat
## GSE140228_UMI_counts_Droplet.mtx.gz > "matrix.mtx.gz"
## GSE140228_UMI_counts_Droplet_barcodes.tsv.gz > "barcodes.tsv.gz"
## GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz > "cellinfo.tsv.gz"
## GSE140228_UMI_counts_Droplet_genes.tsv.gz - remove columns except 1:2 and 1st raw. then compress as "features.tsv.gz"

# This script sets up the Seurat Objects
# Standard pre-processing
library(Seurat)
library(ggplot2)
library(cowplot)

hcc.count <- Read10X(data.dir = "input/GSE140228_UMI_counts_Droplet/")
hcc <- CreateSeuratObject(counts = hcc.count, project = "GSE140228")
rm(hcc.count)

sample <- read.table("input/GSE140228_UMI_counts_Droplet/cellinfo.tsv.gz",sep="\t",header=T,row.names=1)

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
rm(list = c("seurat_pre", "seurat_post", "hcc_filt", "sample", "p1", "p2", "p3", "p4"))

png("plots/qc_summary_by_conditions.png", width = 800*2, height = 600*2.5, res = 200)
p1 <- VlnPlot(hcc, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by = "donor", pt.size = 0.01)
p2 <- VlnPlot(hcc, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by = "tissue", pt.size = 0.01)
plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          labels = c("Donor", "Tissue"), ncol = 1) +
  draw_figure_label(label = "GSE140228_HCC", position = "top.right", fontface = "bold", size = 10)
dev.off()

rm(list = c("p1", "p2"))


# Normalizing the data and identification of variable feaures
hcc <- NormalizeData(hcc)
hcc <- FindVariableFeatures(hcc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(hcc), 10)

png("plots/VariableFeatures.png", width = 800*2, height = 600*2, res = 200)
p1 <- VariableFeaturePlot(hcc)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
p2
dev.off()

rm(list = c("p1", "p2", "top10"))


# Scaling the data and linear dimensional reduction
hcc <- ScaleData(hcc, verbose = TRUE)
hcc <- RunPCA(hcc, npcs = 30, seed.use = 1)

saveRDS(hcc, file = "processed/hcc_scale.RDS")


# NOTE: This process can take a long time for big dataset
hcc <- JackStraw(hcc, num.replicate = 100, dims = 30)
hcc <- ScoreJackStraw(hcc, dims = 1:30)

JackStrawPlot(hcc, dims = 1:30)
ElbowPlot(hcc, ndims = 30)


# Visualizing PCA, JacStrawplot and Elbowplot
pdf("plots/PCA.pdf")

print(DimHeatmap(hcc, dims = 1:9, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 10:20, cells = 500, balanced = TRUE))
print(DimHeatmap(hcc, dims = 20:30, cells = 500, balanced = TRUE))

dev.off()


png("plots/JackStraw.png", width = 800*2, height = 600*2, res = 200)
JackStrawPlot(hcc, dims = 1:30)
dev.off()

png("plots/ElbowPlot.png", width = 800*2, height = 600*2, res = 200)
ElbowPlot(hcc, ndims = 30)
dev.off()


# Clustering and visualization
n_dim <- 28
hcc <- FindNeighbors(hcc, reduction = "pca", dims = 1:n_dim)
hcc <- FindClusters(hcc, random.seed = 1, resolution = 0.9)

hcc <- RunTSNE(hcc, dims = 1:n_dim, seed.use = 1)
hcc <- RunUMAP(hcc, assay = "features", dims = 1:n_dim, seed.use = 1)

png("plots/DimPlot.png", width = 800*3.2, height = 600*3, res = 200)
p1 <- DimPlot(hcc)
p2 <- DimPlot(hcc, group.by = "donor")
p3 <- DimPlot(hcc, group.by = "celltype")
p4 <- DimPlot(hcc, group.by = "tissue")
plot_grid(p1, p2, p3, p4, nrow = 2) +
  draw_figure_label(label = "GSE140228_HCC", position = "top.right", fontface = "bold", size = 10)
dev.off()

saveRDS(hcc, file = "processed/hcc.RDS")

