# Spatiotemporal co-transcriptomics uncovers coordinated rhizobia-legume stress programs supporting nitrogen fixation in Medicago

This repository contains scripts used for spatial transcriptomic analysis of rhizobium-legume interactions in Medicago. The scripts cover data preprocessing, clustering, visualization, trajectory inference, cell-type deconvolution, gene co-expression analysis, and regulatory network exploration.

## Scripts

* **Giotto.R**: R script for preprocessing spatial transcriptomic data using the Giotto framework. It generates the basic Giotto object used in downstream analyses.

* **Harmony.R**: R script for batch correction and data integration using Harmony. It removes technical variation between samples and generates integrated embeddings for clustering and visualization.

* **Cell_cluster.R**: R script for unsupervised clustering and cell-type identification. It performs clustering analysis and generates UMAP visualizations of identified cell populations.

* **Plot_UMAP.R**: R script for UMAP visualization of spatial transcriptomic datasets. It produces dimensional reduction figures colored by cluster, sample, or annotation.

* **Dotplot.R**: R script for visualization of marker gene expression across cell types. It generates dot plots for selected genes and cell populations.

* **AverageExpression.R**: R script for calculating and visualizing average gene expression across clusters or cell types. It generates heatmaps and expression summary plots.

* **RCTD.R**: R script for robust cell-type decomposition using RCTD. It integrates single-cell reference datasets with spatial transcriptomic data and generates cell-type proportion maps.

* **Monocle3.R**: R script for pseudotime trajectory analysis using Monocle3. It reconstructs developmental trajectories and generates pseudotime visualization figures.

* **RNA_Velocity.py**: Python script for RNA velocity analysis. It estimates transcriptional dynamics and generates velocity stream and trajectory plots.

* **Mfuzz.R**: R script for fuzzy c-means clustering of dynamic gene expression patterns. It identifies temporal expression modules and generates cluster trend visualizations.

* **WGCNA.R**: R script for weighted gene co-expression network analysis (WGCNA). It identifies co-expression modules and hub genes associated with developmental processes.

* **Sankey.R**: R script for Sankey diagram visualization. It illustrates relationships among cell types.

