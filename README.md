# scRNA seq data analysis using Seurat V5 package
This repository is a reference material for scRNA analysis course. The reference materials and important links are in my google classroom:
https://classroom.google.com/c/NzcxODY5MzM3OTkx?cjc=usnpoawr

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
