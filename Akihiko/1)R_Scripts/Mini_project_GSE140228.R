#Daily setup==== 
proj_dir <- "/Users/kichu/Documents/Bioinformatics/Bioinfo_course_AK/Akihiko" 
data_dir <- "/Users/kichu/Documents/Bioinformatics/Bioinfo_course_AK/Akihiko/2)Data"

pacman::p_load(Seurat, magrittr, dplyr, cowplot, patchwork, tidyverse)

#Section 00_Prepare data files for Read10X====
#Put 10X data files in data_dir
#File name should be
  #barcodes.tsv, genes.tsv, and matrix.mtx (CellRanger < 3.0)
  #matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz (CellRanger >= 3.0)
  #GSE140228 is processed using CellRanger v2.3

setwd(data_dir)
barcodesfile <- "GSE140228_UMI_counts_Droplet_barcodes.tsv" 
genesfile <- "GSE140228_UMI_counts_Droplet_genes.tsv"
mtxfile <- "GSE140228_UMI_counts_Droplet.mtx" 

dir.create("10X files for Read10X")
data_dir <- paste0(data_dir, "/10X files for Read10X")

file.copy(from = barcodesfile, to = paste0(data_dir, "/barcodes.tsv"))
file.copy(from = genesfile, to = paste0(data_dir, "/genes.tsv"))
file.copy(from = mtxfile, to = paste0(data_dir, "/matrix.mtx"))

rm(barcodesfile, genesfile, mtxfile)

#genes.tsv should have only 3 columns:Gene Id, Gene name, feature type (from V3).
#genes.tsv doesn't expect header

#As for the genes.tsv of GSE140228,
##Delete the header
##Keep only the left two columns (gene id, gene name)
##Delete unnecessary columns

setwd(data_dir)
genesfile <- read.delim("genes.tsv", header = FALSE, skip = 1)[1:2]
  #Do not recognize the header
  #Skip 1 to remove the first line (previously the header)
  #[1:2] removes columns 3-7 right away
write.table(genesfile, file = "genes.tsv", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
rm(genesfile)

#Read data using Read10X====
expression_matrix <- Read10X(data.dir = data_dir)

#Inspect the data
dim(expression_matrix)

#The data should be in the form of each row being a gene and each column being a cell
expression_matrix[1:5,1:5]

#Initialize the Seurat object with the raw data====
seurat_object <-  CreateSeuratObject(counts = expression_matrix, project = "GSE140228")
rm(expression_matrix)

#Save Data_"00_Initial_Seurat_Object"====
setwd(proj_dir)
saveRDS(object = seurat_object, file = "./3)Results/00_Initial_Seurat_Object.rds") 
seurat_object <- readRDS(file = "./3)Results/00_Initial_Seurat_Object.rds")

#Section 01_Initial QC====
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(seurat_object@meta.data, 5)

#Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

qc_initial <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

setwd(proj_dir)
pdf("./3)Results/01__Initial QC.pdf",width=8.27, height=5.84)
print(qc_initial)
dev.off()
rm(qc_initial)

#Filter out poor quality cells
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

#Standard processing_Data normalization====

#SCTransform
#seurat_object <- SCTransform(seurat_object)  
#Error: vector memory exhausted (limit reached?)
###Question###
  #In what cases should I choose SCTransform instead of LogNormalize?
  #Approximately how much memory is needed for SCTransform?
  #My laptop has 16GB Apple M1 chip.

#LogNormalize
seurat_object <- NormalizeData(seurat_object)

#Save Data_"01_Bef_feature_selection"====
setwd(proj_dir)
saveRDS(object = seurat_object, file = "./3)Results/01_Bef_feature_selection.rds") 
seurat_object <- readRDS(file = "./3)Results/01_Bef_feature_selection.rds")

#Section 02_Standard processing_Feature selection====
seurat_object <- FindVariableFeatures(seurat_object, 
                                      selection.method = "vst", 
                                      nfeatures = 4000) #defaultは2000

# Identify the 10 most highly variable genes
top10 <- VariableFeatures(seurat_object)[1:10]

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
variable_features <- plot_grid(plot1, plot2, ncol = 1) 

setwd(proj_dir)
pdf("./3)Results/02__Variable features.pdf",width=8.27, height=11.69)
print(variable_features)
dev.off()
rm(plot1, plot2, variable_features, top10)

#Section 03_PCA====

#Scaling for all genes
#seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
#Error: vector memory exhausted (limit reached?)
###Question###
#Approximately how much memory is needed for this?

#Scaling for previously determined variable features (4,000 at this time)
seurat_object <- ScaleData(seurat_object)

#PCA on scaled data
seurat_object <- RunPCA(seurat_object)

DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE, combine = FALSE)

###Question###
#How can I export this combined figure as a pfd (or png) file using a script?

#Save Data_"03_Aft_PCA"====
setwd(proj_dir)
saveRDS(object = seurat_object, file = "./3)Results/03_Aft_PCA.rds") 
seurat_object <- readRDS(file = "./3)Results/03_Aft_PCA.rds")

#Section 04_Jackstraw_To estimate the significance of each PC====
seurat_object <- JackStraw(seurat_object, num.replicate = 80) 
  #num.replicate = 100 resulted in Error (session aborted)

seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)

#jackstrawplot and/or elbowplot

jsp <- JackStrawPlot(seurat_object, dims = 1:20)　　
jsp

#Warning message:Removed 56000 rows containing missing values (geom_point). 
###Question###
# Is it okay to ignore missing values?
# p value is zero even for PC20.
# PC1 to PC20 are overlapped on the jack straw plot. 
# Is this true?

#elbow plot
ep <- ElbowPlot(seurat_object)
ep

setwd(proj_dir)
jsp_ep <- plot_grid(jsp, ep, ncol = 1) 
pdf("./3)Results/04__Jackstraw for PC selection.pdf",width=8.27, height=11.69)
print(jsp_ep)
dev.off()
rm(jsp, ep, jsp_ep)

#Save Data_"04_Aft_jackstraw"====
setwd(proj_dir)
saveRDS(object = seurat_object, 
        file = "./3)Results/04_Aft_jackstraw.rds") 
seurat_object <- readRDS(file = "./3)Results/04_Aft_jackstraw.rds")

#Section 05_Select PCs (Dims)====
#The first 30 features/genes that contribute to PC11-PC20 are shown.
pc_related_genes_1 <- VizDimLoadings(seurat_object, dims = 11:15, 
                                   nfeatures = 30, reduction = "pca", 
                                   ncol = 5)
pc_related_genes_2 <- VizDimLoadings(seurat_object, dims = 16:20, 
                                     nfeatures = 30, reduction = "pca", 
                                     ncol = 5)

setwd(proj_dir)
pdf("./3)Results/05__PC_related_genes.pdf",width=11.69, height=8.27)
print(pc_related_genes_1)
print(pc_related_genes_2)
dev.off()
rm(pc_related_genes_1, pc_related_genes_2)

#Determine the dims for clustering
selected_dims = 1:15  

#Run tSNE and UMAP====
seurat_object <- RunTSNE(seurat_object, dims = selected_dims)
seurat_object <- RunUMAP(seurat_object, dims = selected_dims)

#Save Data_"05_Bef_FindClusters"====
setwd(proj_dir)
save(seurat_object, data_dir, proj_dir, selected_dims, file = 
       "./3)Results/05_Bef_FindClusters.Rdata")
load("./3)Results/05_Bef_FindClusters.Rdata")

#Section 06_Clustering at various resolutions====
seurat_object <- FindNeighbors(seurat_object, dims = selected_dims)

res.examined <- c(3.0, 2.0, 1.0, 0.5) #resolutions to be tested
seurat_object <- FindClusters(seurat_object, resolution = res.examined)

#TSNE plot with clustering
colnames(seurat_object@meta.data)
clusters.res <- paste("RNA_snn_res.", res.examined, sep = "")
clusters.res

res_compare_tSNE <- DimPlot(
  seurat_object, reduction = "tsne", group.by = clusters.res
  ) + 
  guides(colour = guide_legend(override.aes = list(size=3), ncol=3))

res_compare_tSNE

###Question###
# guides() is only applied to the last tSNE plot (res0.5).
# How can I apply guides() to the legends of all 4 plots?

#UMAP plot with clustering
res_compare_UMAP <- DimPlot(
  seurat_object, reduction = "umap", group.by = clusters.res
  ) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=3)) 

setwd(proj_dir)
png("./3)Results/06__res_compare_tSNE.png",width=4093, height=2894, res=350)
print(res_compare_tSNE)
dev.off()
png("./3)Results/06__res_compare_UMAP.png",width=4093, height=2894, res=350)
print(res_compare_UMAP)
dev.off()
rm(res_compare_tSNE, res_compare_UMAP)

#See how many cells end up in each cluster
table1 <- table(seurat_object[[clusters.res[1]]])
table2 <- table(seurat_object[[clusters.res[2]]])
table3 <- table(seurat_object[[clusters.res[3]]])
table4 <- table(seurat_object[[clusters.res[4]]])

n <- max(length(table1), length(table2), 
         length(table3), length(table4))

length(table1)=n                     
length(table2)=n
length(table3)=n
length(table4)=n

cell.count_cluster <- cbind(table1, table2, table3, table4)
column_name <- paste("res=", res.examined, sep = "")
colnames(cell.count_cluster) <- column_name
row_name <- paste("cluster", rownames(cell.count_cluster), sep = "")
rownames(cell.count_cluster) <- row_name
write.csv(cell.count_cluster, file = "./3)Results/06__Cell count_cluster.csv")
rm(n, cell.count_cluster, column_name, row_name, res.examined, 
   table1, table2, table3, table4)

###Question###
#Is there any smarter way to compare the clustering at various resolutions?


#Save Data_"06_Aft_FindClusters"====
setwd(proj_dir)
save(seurat_object, data_dir, proj_dir, selected_dims, clusters.res, file = 
       "./3)Results/06_Aft_FindClusters.Rdata")
load("./3)Results/06_Aft_FindClusters.Rdata")

#Section 7_Find markers for the clusters identified at the determined resolution====
#Determine the resolution for clustering
head(seurat_object@meta.data)
clusters.res
Idents(seurat_object) <- seurat_object[[clusters.res[2]]] 
                                    #In the case of res=2.0

#Find markers for clusters====
all.cluster_markers <- FindAllMarkers(seurat_object, test.use = "bimod", 
                               only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
colnames(all.cluster_markers)
head(all.cluster_markers)

#Filter markers with an adjusted p value of less than 0.05
all.cluster_markers <- 
  all.cluster_markers[ all.cluster_markers$p_val_adj < 0.05, ]
 
#pick the top 50 markers
top_50_markers <- 
  all.cluster_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 50, order_by = avg_log2FC) 

top_50_markers <- as.data.frame(top_50_markers)

#save to tab delimited file
write.table(top_50_markers, 
            file="./3)Results/07__top50_markers.txt", 
            quote=F, sep='\t', col.names = NA)

# # Finding all markers of a specific cluster 
# cluster2.markers <- FindMarkers(seurat_object, test.use = "bimod", 
#                                 ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)

# Finding markers between clusters
#To find differentially regulated genes among specific groups
#Here you specify specific clusters to find markers using the ident.1 
#and ident.2 inputs

c1_vs_c2_markers <- FindMarkers(seurat_object, test.use = "bimod", 
                                ident.1 = 1,  ident.2 = 2, min.pct = 0.25)
c1_vs_c2_markers <- c1_vs_c2_markers[ c1_vs_c2_markers$p_val_adj < 0.05 , ]

c5_vs_c0_and_c3 <- FindMarkers(seurat_object, test.use = "bimod", 
                               ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

# To take a look at top DEGs
head(c1_vs_c2_markers, n = 5)

#More details on DEG testing: 
#https://satijalab.org/seurat/v3.0/de_vignette.html

#Visualization of the marker expression====

#Heatmap
top_10_markers <- all.cluster_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
DoHeatmap(seurat_object, features = top_10_markers$gene) + NoLegend()
#scalingされていないgeneは自動的にomitされる

#Looking at the expression of selected genes

##Using violin plots
###Showing normalized gene expression
VlnPlot(seurat_object, features = c("MS4A1", "CD79A"))
###Showing raw counts
VlnPlot(seurat_object, features = c("NKG7", "PF4"), slot="counts", log=TRUE)

##Using intensity plots (UMAP上に色づけ)
FeaturePlot(seurat_object, features = c("CD3E", "PTPRC", "FCER1A", "LYZ"))

#to save to a file directly:
intensity_plot <- FeaturePlot(seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
setwd(proj_dir)
png("./3)Results/marker_intensity.png",width=4093, height=2894, res=350)
print(intensity_plot)
dev.off()
#For png, height and width denote pixels
#You can also save a pdf using the pdf command instead of png 
#command, but for saving pdf, you specify the height and width 
#in inches. Please do not mix them up!

#Method 1: 
tsne_plot <- DimPlot(object = seurat_object, reduction = 'tsne', 
                       group.by = 'seurat_clusters')
png("./pbmc_tsne.png",width=1000, height=800, res=144)
print(tsne_plot)
dev.off()
#You can also save as pdf 
#using pdf() and with paper sizes instead of image resolution.

#Putting the intensity side by side with the clustering
p1 <- DimPlot(object = seurat_object, reduction = 'tsne', group.by = 'seurat_clusters')
p2 <- FeaturePlot(seurat_object, reduction = 'tsne', features = c("MS4A1"))
plot_grid(p1, p2)

#Adding annotation====

new_cluster_ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

names(new_cluster_ids) <- levels(seurat_object)
#levels(x)はobjectのlevelsの値を返す
#names(x)はobject nameをsetあるいはget
seurat_object <- RenameIdents(seurat_object, new_cluster_ids)
#RenameIdent(object, old.ident.name = NULL, new.ident.name = NULL)

#Visualize with the new labels
DimPlot(seurat_object, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
#plot上にcluster name入れるならこっち
DimPlot(seurat_object, reduction = "tsne", pt.size = 0.5) 
#凡例でつけるならこっち

#Saving and loading seurat objects====

#To save a seurat object:
#saveRDS(seurat_object, file = "pbmc3k_final.rds") #ファイル名入力

#Warning: seurat objects can be quite big!
#For this exercise, the saved object is > 274MB
#There is a function called DietSeurat that can be used to remove some of the 
#assays inside the Seurat object to save space.


#To load a seurat object:
#seurat_object <- readRDS(file = "pbmc3k_final.rds")

#####################################
# Very useful list of commands for interacting with the Seurat object:
# Seurat Standard Workflow
# Seurat Object Interaction
# Data Access
# Visualization
# Multi-Assay Feature
# Seurat v2.X vs 3.X

#https://satijalab.org/seurat/essential_commands.html
# 
# If you find severe batch effects in your data, you may want to perform batch integration/correction.
# The main aim of this is to find similar cell types across different batches and “bring them together”
# For the Seurat package, the batch integrated output should not be used for DEG analysis.
# Instead, use the identified cell clusters (composed of cells from different batches) for DEG analysis.
# 
# For details, see the vignette on the Seurat website, or the batch correction script for this course.



