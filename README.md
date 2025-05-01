# scRNA seq data analysis using Seurat V5 package
This repository is a reference material for scRNA analysis course. The reference materials and important links are in my google classroom:
https://classroom.google.com/c/NzcxODY5MzM3OTkx?cjc=usnpoawr
---
---

## 1. scRNA_Analysis_Guided_Clustering_Tutorial-pbmc3k
For this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are 2,700 single cells that were sequenced on the Illumina NextSeq 500.

The **dataset** for this tutorial can be downloaded from here:
https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

📊 **Seurat Object Structure**
The tutorial is built around the Seurat object, which serves as a central container for:

- Raw & normalized data

- Cell-level metadata

- Dimensionality reductions (PCA, UMAP, etc.)

- Clustering results

- Graph-based structures (kNN, SNN)
---------------------------------------------------------------------------------------------------------------
 ### 🧪 Analysis Workflow Covered
- *Data Import & Seurat Object Creation*
  
    CreateSeuratObject()

- *Quality Control (QC)*
  
    Filter cells by nFeature_RNA, nCount_RNA, percent.mt

- *Normalization*
  
    Using NormalizeData() or SCTransform()

- *Highly Variable Feature Selection*
  
    FindVariableFeatures()

- *Scaling & PCA*
  
    ScaleData(), RunPCA()

- *Clustering*
  
    FindNeighbors(), FindClusters()

- *Dimensionality Reduction & Visualization*
  
    RunUMAP(), DimPlot()

- *Marker Gene Identification & Cluster Annotation*
  
    FindAllMarkers()
---
---

## 2. Data visualization

**Title:**  
🧬 **Visualizing Single-Cell Data in Seurat: A Comprehensive Guide**

**Summary:**  
This Seurat visualization vignette demonstrates how to effectively explore and present single-cell RNA-seq data using the Seurat R package. It covers key plotting functions such as `DimPlot`, `FeaturePlot`, `VlnPlot`, `RidgePlot`, `DotPlot`, and `Heatmap`, providing insight into cell clustering, gene expression, and cell identity.

The tutorial showcases how to visualize dimensionality reductions (UMAP/tSNE), display feature-level expression across cells, and create comparative plots across conditions or metadata groups.

The guide emphasizes clarity, aesthetics, and the power of patchwork for combining multiple plots into cohesive figures—making it essential for publication-quality visual analytics in scRNA-seq research.
