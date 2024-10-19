# Joint-scRNA-scATAC-analysis

This project processes multiomic data (RNA-seq and ATAC-seq) from a 10X Genomics PBMC granulocyte sorted dataset. It includes data downloading, preprocessing, peak calling, differential accessibility analysis, and visualization of gene expression and chromatin accessibility correlations.

##Overview
This pipeline performs joint RNA-seq and ATAC-seq analysis using Seurat, ArchR, and other bioinformatics tools. The pipeline consists of the following steps:

- Download raw data from 10X Genomics.
- Preprocess RNA-seq and ATAC-seq data.
- Perform dimensionality reduction and clustering.
- Identify differentially accessible peaks across cell types.
- Visualize gene expression and chromatin accessibility correlations using Seurat and Signac.

P.S. I plan to set up a Docker container to run the analysis at a later time.




