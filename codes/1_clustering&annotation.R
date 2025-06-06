library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "D:/youtube/bioinformatics/Single Cell Analysis/Hands_on/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#####################################################################
########QC and selecting cells for further analysis#################
##################################################################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# This line of code is calculating mitochondrial gene expression percentage per cell and storing it in the Seurat object as metadata.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

##The number of unique genes and total molecules are automatically calculated during - CreateSeuratObject()
##You can find them stored in the object meta data

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

#####################################################################

###################### Visualization #################################

##################################################################


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#####################################################################

###################### NORMALIZATION #################################

##################################################################

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)



#####################################################################

###################### NORMALIZATION #################################

##################################################################
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#####################################################################

###################### SCALING #################################

##################################################################

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = "percent.mt")
#pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

#####################################################################

###################### LINEAR DIMENSIONALITY REDUCTION #################################

##################################################################

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


#####################################################################

###################### Determine the ‘dimensionality’ of the dataset #################################

##################################################################

ElbowPlot(pbmc)



#####################################################################

###################### Clustering #################################

##################################################################

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#####################################################################

###################### NON LINEAR DIMENSIONALITY REDUCTION #################################

##################################################################

pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


# save the object
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


#####################################################################

## Finding differentially expressed features (cluster biomarkers) ###

##################################################################

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

########################################################################
#########Finding Marker Genes for Cluster 0 using ROC Test#############
######################################################################
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc,
        features = c("MS4A1", "CD79A"),
        pt.size = 0.2,             # smaller dots for better clarity
        cols = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
                 "#6a3d9a", "#b15928", "#a6cee3", "#fb9a99", "#fdbf6f"),  # nice color palette
        group.by = "seurat_clusters",  # group by cluster identities
        slot = "data")             # raw expression values

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

pbmc.markers %>%
  group_by(cluster) %>%                     # 1. Group all markers by cluster
  dplyr::filter(avg_log2FC > 1) %>%         # 2. Keep only strong markers (log2FC > 1)
  slice_head(n = 10) %>%                    # 3. Pick top 10 genes per cluster
  ungroup() -> top10                        # 4. Save as `top10`

#This draws a heatmap of expression levels of the top 10 most enriched genes (strong markers) per cluster.
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

################################################################
#Assigning cell type identity to clusters##
#######################################################

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

#saveRDS(pbmc, file = "../output/pbmc3k_final.rds")

# Or you can save the entire workspace

save.image(file = "../output/pbmc_workspace.RData")


