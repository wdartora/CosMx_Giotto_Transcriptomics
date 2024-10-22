# Spatial Transcriptomics Pipeline Using Giotto <br>
This pipeline is designed for the analysis of spatial transcriptomics data using the Giotto package in R. It processes CosMx subcellular transcriptomics data, performs data visualization, and conducts clustering and dimension reduction analyses. The pipeline includes various steps such as data loading, feature extraction, normalization, clustering, and visualization.
 <br> <br>
## Setup <br>
Install required packages: <br>
Install and load the Giotto package along with any other dependencies like data.table, and ensure a working Python installation with required libraries. <br>
 <br>
Configure paths: <br>
## Set the working directory and provide paths for: <br>
 <br>
Data files (transcript coordinates, field of vision positions, etc.) <br>
Output directory for saving results. <br>
Python environment: <br>
Optionally, specify the path to a Python executable within a conda environment. <br>
 <br>

## Steps in the Pipeline <br>
1. Create Giotto CosMx Object <br>
Loads spatial transcriptomics data from the specified directory. <br>
Creates a Giotto object containing subcellular transcript data and field of views (FOVs). <br>
2. Data Exploration and Feature Extraction <br>
Loads and inspects subcellular detection information. <br>
Filters and separates features from negative probes. <br>
Visualizes the spatial distribution of features and probes. <br>
3. Create Giotto Objects for FOVs <br>
Processes each FOV, creating Giotto objects that include subcellular data, segmentation masks, and feature coordinates. <br>
4. Join FOVs <br>
Combines data from multiple FOVs into a single Giotto object, aligning them based on provided offsets. <br>
5. Visualize Cells and Features <br>
Generates various spatial and in situ plots of feature distributions across cells. <br>
Visualizes cell centroids and feature-specific expression. <br>
6. Aggregate Features and Normalize Data <br>
Aggregates subcellular features into cell-level data. <br>
Performs data filtering and normalization (log-normalization and Pearson residuals method). <br>
7. Dimension Reduction <br>
Detects highly variable genes and performs PCA and UMAP. <br>
Plots PCA and UMAP results to visualize feature distributions. <br>
8. Clustering <br>
Constructs a nearest-neighbor network and applies Leiden clustering. <br>
Visualizes clustering results spatially and across dimension-reduced space. <br>
9. Small Subset Visualization <br>
Extracts a small region of interest from the Giotto object based on spatial coordinates. <br>
Visualizes genes and clustering results within the subset. <br>
Custom Color Palettes <br>
The pipeline uses custom color palettes for plotting, defined at the beginning of the script (pal10, viv10, pal13). You can modify these to customize the appearance of the plots.
 <br>


## The pipeline produces a variety of outputs: <br>
 <br>
Plots for feature and cell distribution in spatial and dimension-reduced space. <br>
Aggregated and normalized expression matrices. <br>
Clustered cell populations, visualized through UMAP and spatial plots. <br>
Saving Results <br>
Plots and output files are automatically saved to the specified results_folder. You can adjust the save directory and other plot options using the instrs object, which controls the behavior of Giotto's plotting functions.
