library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#Load the data
#pbmc3k_final <- readRDS("D:/youtube/bioinformatics/Single Cell Analysis/Hands_on/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/output/pbmc3k_final.rds")

load("../output/pbmc_workspace.RData")


######################################################################################
##############Plotting variability in detected molecules across cells#################
######################################################################################

# Plotting a histogram of the total detected molecules (UMIs)
ggplot(pbmc@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Total Molecules Detected per Cell",
       x = "nCount_RNA (Total UMIs)",
       y = "Number of Cells")

# Violin plot of total detected molecules per cell
VlnPlot(pbmc, features = "nCount_RNA", pt.size = 0.1) +
  ggtitle("Variability in Total Molecules Detected per Cell")

# Scatter Plot: Genes vs Molecules
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1) +
  ggtitle("UMI Count vs Number of Genes Detected")

#######################################################################################
#############################Exploring the Seurat object###############################
#######################################################################################

# To view the assays in your Seurat object pbmc, use:
# This will return a character vector listing all assay names (e.g., "RNA", "SCT").
Assays(pbmc)

#To get more details about a specific assay (e.g., "RNA"):
pbmc[["RNA"]]

# Or to see available slots inside an assay object:
slotNames(pbmc[["RNA"]])

# For example, to look at the expression matrix in the default layer:
pbmc[["RNA"]]@layers
pbmc[["RNA"]]@default # is accessing the default assay slot inside the pbmc Seurat object.
pbmc[["RNA"]]@layers[[pbmc[["RNA"]]@default]]

# To check cells and features names:
pbmc[["RNA"]]@cells
head(pbmc[["RNA"]]@features)

pbmc[["RNA"]]@cells
head(pbmc[["RNA"]]@features)

# what are the "3" observations?
pbmc[["RNA"]]@meta.data



#######################################################################################
#######################Exploring through Visualizations################################
#######################################################################################

pbmc$groups <- Idents(pbmc) # Idents() retrieves the current identity class (e.g., clusters like 0, 1, 2...).



############## Explore all genes (features) in the Seurat object #########################
all_genes <- rownames(pbmc)
length(all_genes)  # Total number of genes
head(all_genes, 20)  # Preview first 20 genes


############### Identify top marker genes per cluster ##############################


# Filter markers with strong expression difference, e.g. avg_log2FC > 1
# top10 <- pbmc.markers %>%
#   group_by(cluster) %>%
#   filter(avg_log2FC > 1) %>%
#   slice_max(order_by = avg_log2FC, n = 10) %>%  # Top 10 per cluster
#   ungroup()

# View top genes for all clusters
top_genes <- unique(top10$gene)
print(top_genes)


# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc, features = c("LYZ", "CCL5"), group.by = "groups")

# Visualize top genes for a few clusters (adjust as needed)
RidgePlot(pbmc, features = top_genes[1:10], group.by = "groups")

# What a RidgePlot Shows:
#   X-axis: Expression level of the gene (usually log-normalized counts).
# 
# Y-axis: Different groups or clusters (ordered vertically).
# 
# Density curves: Show the distribution shape of expression values for that gene within each group.


# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc, features = top_genes) # top10_genes is a character vector whereas top10 is a dataframe

# our top10 is a data frame, not a simple list of genes. The VlnPlot() function requires a simple vector of gene names, not a data frame.
# 
# How to fix:
#   
#   Extract the gene names from top10 as a vector and pass that to VlnPlot():

VlnPlot(pbmc, features = top_genes[1:10])
VlnPlot(pbmc, features = top_genes[1:5])


# Get top 10 most variable genes identified by FindVariableFeatures()
top_variable_genes <- head(VariableFeatures(pbmc), 10)
# Plot violin plots of these genes
VlnPlot(pbmc, features = top_variable_genes)


# top10 is a data frame containing marker genes along with other information (like cluster, log fold change, etc.).
# 
# top10$gene accesses the gene column from this data frame.
# 
# unique() removes any duplicate gene names, giving you a vector of distinct gene names.
# 
# This vector is stored in the variable top10_genes.


#### To check how many cell types (clusters or identities) are detected in your Seurat object, run:
levels(pbmc)
length(unique(Idents(pbmc)))


head(pbmc@meta.data)

# To check unique groups in a specific metadata column
unique(pbmc$groups)

unique(Idents(pbmc))

VlnPlot(pbmc, features = c("CD3D", "MS4A1"), group.by = "seurat_clusters")

# Feature plot - visualize feature expression in low-dimensional space
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
FeaturePlot(pbmc, features = features)
#requires understanding how expression is mapped onto the dimensionality-reduced space
# Light purple or white → low or no expression of the gene in those cells.
# 
# Dark purple → high expression in those cells.
# 
# Each dot represents a cell, and its position is based on UMAP (i.e., cells closer together are more transcriptionally similar).


# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(pbmc, features = features) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(pbmc, downsample = 100), features = features, size = 3)


# Randomly selects 100 cells per cluster to prevent overplotting.
# When downsampling is not ideal:
#   Situation	Better Alternative
# Rare cell populations	Don't downsample—retain all cells in rare clusters.
# Analyzing subcluster heterogeneity	Use full data or cluster-specific plots.
# Feature-level publication (e.g., showing exhaustive heatmaps)	Use pseudobulk averages or group medians per cluster.
# Quantitative, not visual analysis	Use expression matrices directly, not heatmaps.
# 


# Plot a legend to map colors to expression levels
FeaturePlot(pbmc, features = "MS4A1")
#gives you spatial and quantitative expression information for the gene MS4A1, projected onto a 2D UMAP (or PCA/tSNE) plot.

# Adjust the contrast in the plot
FeaturePlot(pbmc, features = "MS4A1", min.cutoff = 1, max.cutoff = 3)

# Cells with expression ≥ 3 are capped at 3 on the scale.
# 
# Anything above 3 gets mapped to the brightest color.
# 
# Prevents outlier cells from skewing the color gradient.

# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
FeaturePlot(pbmc, features = c("MS4A1", "PTPRCAP"), min.cutoff = "q10", max.cutoff = "q90")
# This FeaturePlot uses the 10th and 90th percentiles to scale color intensity, reducing the effect of outliers
# and enhancing contrast. It helps visualize the expression of MS4A1 and PTPRCAP more clearly across cells.
# This is a good practice for more interpretable and balanced plots.


# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc, features = c("MS4A1", "CD79A"), blend = TRUE)


# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc, features = features, split.by = "groups", 
        cols = c("blue", "red", "green", "purple", "orange", "brown", "pink", "cyan", "yellow")) + RotatedAxis()

library(RColorBrewer)
DotPlot(pbmc, features = features, split.by = "groups", 
        cols = brewer.pal(9, "Set1")) + RotatedAxis()



###############################################################################################
###############################Applying Themes to plots########################################
###############################################################################################

baseplot <- DimPlot(pbmc, reduction = "umap")
# Add custom labels and titles
baseplot + labs(title = "Clustering of 2,700 PBMCs")

# Seurat also provides several built-in themes
baseplot + DarkTheme()


###############################################################################################
###############################Interative Feature Plots########################################
###############################################################################################

# Include additional data to display alongside cell names by passing in a data frame of
# information.  Works well when using FetchData
plot <- FeaturePlot(pbmc, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))


# Seurat allows manual cell selection from plots like DimPlot() or FeaturePlot() using CellSelector().
# This is valuable for identifying small or subtle clusters that automatic clustering might miss.
# By selecting cells interactively, you can assign them a new identity and perform targeted differential
# expression analysis.

# For example, let’s pretend that DCs had merged with monocytes in the clustering, but we wanted to see
# what was unique about them based on their position in the tSNE plot.

pbmc3k <- RenameIdents(pbmc, DC = "CD14+ Mono")
plot <- DimPlot(pbmc, reduction = "umap")
select.cells <- CellSelector(plot = plot)

# We can then change the identity of these cells to turn them into their own mini-cluster.
head(select.cells)
Idents(pbmc3k, cells = select.cells) <- "NewCells"

# Now, we find markers that are specific to the new cells, and find clear DC markers
newcells.markers <- FindMarkers(pbmc3k, ident.1 = "NewCells", ident.2 = "CD14+ Mono", min.diff.pct = 0.3,
                                only.pos = TRUE)
head(newcells.markers)

# In addition to returning a vector of cell names, CellSelector() can also take the selected
# cells and assign a new identity to them, returning a Seurat object with the identity classes
# already set. This is done by passing the Seurat object used to make the plot into CellSelector(),
# as well as an identity class. As an example, we’re going to select the same set of cells as before,
# and set their identity class to “selected”
pbmc3k1 <- CellSelector(plot = plot, object = pbmc3k, ident = "selected")
levels(pbmc3k1)
