# Clear all
rm(list = ls())
library(Giotto)


# 1. Setup ---------------------------------------------------------------------->>>>
# Custom color palettes from rcartocolor
# pal10 = rcartocolor::carto_pal(n = 10, name = 'Pastel')
pal10 = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#87C55F",
          "#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B3B3B3")
# viv10 = rcartocolor::carto_pal(n = 10, name = 'Vivid')
viv10 = c("#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0",
          "#24796C","#DAA51B","#2F8AC4","#764E9F","#A5AA99")

pal13 = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#87C55F","#FC8D62", 
            "#8DA0CB", "#E78AC3",
          "#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B3B3B3")
# set working directory
results_folder = 'D:/William/Project/Pipeline_Spatial/pipeline_01/results'

# Optional: Specify a path to a Python executable within a conda or miniconda
# environment. If set to NULL (default), the Python executable within the previously
# installed Giotto environment will be used.
my_python_path = "c:/users/wjd4002/appdata/local/programs/python/python312/python.exe"

## Set object behavior
# by directly saving plots, but not rendering them you will save a lot of time
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE,
                                  python_path = my_python_path)



# 1.1 CosMx Project loading function
## provide path to nanostring folder
data_path = 'D:/William/Project/Pipeline_Spatial/pipeline_01/data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/'

## create giotto cosmx object
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular', # only subcellular
                                   FOVs = c(2,3,4),
                                   instructions = instrs)

showGiottoFeatInfo(fov_join)
showGiottoSpatialInfo(fov_join)


# 2. Data exploration and loading  ------------------------------------------------->>>>
# 2.1 Subcellular detections (points info)
# tx_file.csv contains the subcellular detections information. 
# It contains information on each of the individual feature detections within the sample.
# - fov which FOV the detection happened in
# - cell_ID the ID of the cell the detection happened in
# - x_global_px the global spatial x location in pixels
# - y_global_px the global spatial y location in pixels
# - x_local_px the spatial x location in pixels within the FOV
# - y_local_px the spatial y location in pixels within the FOV
# - z the z plane the detection was called in (-1 to 16)
# - target the feature the probe is targeted against
# - CellComp Cellular compartment the detection happened in (0, Cytoplasm, Membrane, Nuclear)
## provide path to nanostring folder
data_path = 'D:/William/Project/Pipeline_Spatial/pipeline_01/data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/'

# load transcript coordinates
tx_coord_all = data.table::fread(paste0(data_path, 'Lung5_Rep1_tx_file.csv'))

colnames(tx_coord_all)
cat('\n')
# z planes
tx_coord_all[, table(z)]
cat('\n')
# Cell compartment
tx_coord_all[, table(CellComp)]


# 2.2 Split detections by features vs negative probes
all_IDs = tx_coord_all[, unique(target)]
# negative probe IDs
neg_IDs = all_IDs[grepl(pattern = 'NegPrb', all_IDs)]
cat('Negative Probe IDs\n')
neg_IDs
cat('\nFeature IDs\n')
feat_IDs = all_IDs[!all_IDs %in% neg_IDs]
length(feat_IDs)

# split detections
feat_coords_all = tx_coord_all[target %in% feat_IDs]
neg_coords_all = tx_coord_all[target %in% neg_IDs]

cat('\nFeatures: ', feat_coords_all[, .N], '\n',
    'NegProbes: ', neg_coords_all[, .N])

# 2.2.1 Preview negative probes (optional)
# Previewing the probe information can be done by converting to giottoPoints and then 
# using plot(). Here we show a preview of the negative probes.
# Note: if previewing the rna expression information, it is highly recommended to set a 
# subset of features using the feats param. The default is to plot all points, 
# which can be very slow for large data.
neg_points = createGiottoPoints(
  x = neg_coords_all[, .(target, x_global_px, y_global_px)]
)
#png(filename="D:/William/Project/Pipeline_Spatial/pipeline_01/results/neg_probs.png")
plot(neg_points, point_size = 0.2, feats = neg_IDs)
#dev.off()


# 2.3 FOV shifts
# fov_positions_file.csv contains information on the x and y shifts needed 
# in order to put the FOVs tiles together into a cohesive whole. 
# This information is needed during the image attachment and alignment process.

#  load field of vision (fov) positions
fov_offset_file = data.table::fread(paste0(data_path, 'Lung5_Rep1_fov_positions_file.csv'))
fov_offset_file

# 2.4 Choose field of view for analysis
# CosMx data is large and Giotto loads in the subcellular information by FOV. 
# This dataset includes 28 FOVs which can be difficult for most computers to handle at once.
# This tutorial will use FOVs ‘02’, ‘03’, and ‘04’ 
# which correspond to the 3 FOVs visible on the bottom right in the negative probe preview above.
gobjects_list = list()

id_set = c('02', '03', '04')



# 3. Create a Giotto Object for each FOV ------------------------------------------>>>
for(fov_i in 1:length(id_set)) {

  fov_id = id_set[fov_i]


  # 1. original composite image as png
  original_composite_image = paste0(data_path, 'CellComposite/CellComposite_F0', fov_id,'.jpg')

  # 2. input cell segmentation as mask file
  segmentation_mask = paste0(data_path, 'CellLabels/CellLabels_F0', fov_id, '.tif')

  # 3. input features coordinates + offset
  feat_coord = feat_coords_all[fov == as.numeric(fov_id)]
  neg_coord = neg_coords_all[fov == as.numeric(fov_id)]
  feat_coord = feat_coord[,.(x_local_px, y_local_px, z, target)]
  neg_coord = neg_coord[,.(x_local_px, y_local_px, z, target)]
  colnames(feat_coord) = c('x', 'y', 'z', 'gene_id')
  colnames(neg_coord) = c('x', 'y', 'z', 'gene_id')
  feat_coord = feat_coord[,.(x, y, gene_id)]
  neg_coord = neg_coord[,.(x, y, gene_id)]


  fovsubset = createGiottoObjectSubcellular(
    gpoints = list('rna' = feat_coord,
                   'neg_probe' = neg_coord),
    gpolygons = list('cell' = segmentation_mask),
    polygon_mask_list_params = list(
      mask_method = 'guess',
      flip_vertical = TRUE,
      flip_horizontal = FALSE,
      shift_horizontal_step = FALSE
    ),
    instructions = instrs
  )


  # cell centroids are now used to provide the spatial locations
  fovsubset = addSpatialCentroidLocations(fovsubset,
                                          poly_info = 'cell')

  # create and add Giotto images
  composite = createGiottoLargeImage(raster_object = original_composite_image,
                                     negative_y = FALSE,
                                     name = 'composite')

  fovsubset = addGiottoImage(gobject = fovsubset,
                             largeImages = list(composite))


  fovsubset = convertGiottoLargeImageToMG(giottoLargeImage = composite,
                                          #mg_name = 'composite',
                                          gobject = fovsubset,
                                          return_gobject = TRUE)

  gobjects_list[[fov_i]] = fovsubset

}


# 4. Join FOV Giotto Objects -------------------------------------------------->>>>
new_names = paste0("fov0", id_set)

id_match = match(as.numeric(id_set), fov_offset_file$fov)
x_shifts = fov_offset_file[id_match]$x_global_px
y_shifts = fov_offset_file[id_match]$y_global_px

# Create Giotto object that includes all selected FOVs
fov_join = joinGiottoObjects(gobject_list = gobjects_list,
                             gobject_names = new_names,
                             join_method = 'shift',
                             x_shift = x_shifts,
                             y_shift = y_shifts)




# 5. Visualize Cells and Genes of Interest ----------------------------------->>>>>
# When plotting subcellular data, Giotto uses the spatInSituPlot functions. 
# Spatial plots showing the feature points and polygons are plotted using spatInSituPlotPoints().

showGiottoImageNames(fov_join)

# Set up vector of image names
id_set = c('02', '03', '04')
new_names = paste0("fov0", id_set)
image_names = paste0(new_names, '-image')

spatInSituPlotPoints(fov_join,
                     show_image = TRUE,
                     image_name = image_names,
                     feats = list('rna' = c('MMP2', 'VEGFA', 'IGF1R',
                                            'MKI67', 'EPCAM', 'KRT8')),
                     feats_color_code = viv10,
                     spat_unit = 'cell',
                     point_size = 0.01,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.03,
                     save_param = list(base_height = 3,
                                       save_name = '1_inSituFeats'))




# 5.1 Visualize Cell Centroids
# The standard spatPlot2D() function can also be used, but this works off only the 
# aggregated information that is assembled based on the subcellular information. 
# Plotting information based on cell centroids can be done through this function.
spatPlot2D(gobject = fov_join,
           image_name = image_names,
           show_image = TRUE,
           point_shape = 'no_border',
           point_size = 0.01,
           point_alpha = 0.5,
           coord_fix_ratio = 1,
           save_param = list(base_height = 2,
                             save_name = '2_spatCentroids'))


# 6. Aggregate subcellular features
# Giotto supports working directly with the subcellular features in order to 
# generate cell by feature matrices. The data generated this way is then given 
# the spatial unit 'cell'. This workflow is recommended over loading the provided
# cell by feature (aggregated expression) matrix and then including the 
# subcellular information as secondary data.
# When both the raw subcellular information and the pre-made expression matrix 
# are loaded in at the same time, the subcellular data and all data generated 
# from it should be given the spatial unit 'cell' and the pre-generated 
# aggregated information should be given a different spatial unit such as 
# 'cell_agg' to differentiate between the two sources of information.

# In this step, we will be aggregating the feature points of 'rna' and '
# neg_probe' into the 'cell' spatial unit.
# Find the feature points overlapped by polygons. This overlap information is then
# returned to the relevant giottoPolygon object's overlaps slot.
fov_join = calculateOverlapRaster(fov_join, feat_info = 'rna')
fov_join = calculateOverlapRaster(fov_join, feat_info = 'neg_probe')

# Convert the overlap information into a cell by feature expression matrix which
# is then stored in the Giotto object's expression slot
fov_join = overlapToMatrix(fov_join, feat_info = 'rna')
fov_join = overlapToMatrix(fov_join, feat_info = 'neg_probe')

showGiottoExpression(fov_join)


# 6.1 Plot histograms of total counts per cell
filterDistributions(fov_join,
                    plot_type = 'hist',
                    detection = 'cells',
                    method = 'sum',
                    feat_type = 'rna',
                    nr_bins = 100,
                    save_param = list(base_height = 3,
                                      save_name = '3.1_totalexpr'))

filterDistributions(fov_join,
                    plot_type = 'hist',
                    detection = 'cells',
                    method = 'sum',
                    feat_type = 'neg_probe',
                    nr_bins = 25,
                    save_param = list(base_height = 3,
                                      save_name = '3.2_totalnegprbe'))


# 6.2 2D Density Plots
spatInSituPlotDensity(gobject = fov_join,
                      feats = c("MMP2", "VEGFA", "IGF1R",
                                'MKI67', 'EPCAM', 'KRT8'),
                      cow_n_col = 2,
                      save_param = list(base_height = 4,
                                        save_name = '4_inSituDens'))



# 6.3 Extract Data from Giotto Object
# combine cell data
morphometa = combineCellData(fov_join,
                             feat_type = 'rna')
morphometa

# combine feature data
featmeta = combineFeatureData(fov_join,
                              feat_type = c('rna'))
featmeta

# combine overlapping feature data
featoverlapmeta = combineFeatureOverlapData(fov_join,
                                            feat_type = c('rna'))
featoverlapmeta



# 7. Filtering and normalization --------------------------------------------->>>>
# After the expression matrix is generated from the subcellular information, 
# analysis proceeds through data filtering and normalization.
# For the normalization step, we will employ two types.
# standard normalization method: library size normalization and log normalization. 
# This method will produce both normalized and scaled values that are be returned 
# as the ‘normalized’ and ‘scaled’ expression matrices respectively. In this 
# tutorial, the normalized values will be used for generating expression 
# statistics and plotting expression values. The scaled values will be ignored. 
# We will also generate normalized values for the negative probes for 
# visualization purposes during which the library normalization step will be 
# skipped.
# pearson residuals: A normalization that uses the method described in 
# Lause/Kobak et al. 2021. This produces a set of values that are most 
# similar in utility to a scaled matrix and offer improvements to both 
# HVF detection and PCA generation. These values should not be used for 
# statistics, plotting of expression values, or differential expression analysis.

# filter (feat_type = 'rna' by default)
fov_join <- filterGiotto(gobject = fov_join,
                         feat_type = 'rna',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5)

# normalize
# standard method of normalization (log normalization based)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            norm_methods = 'standard',
                            verbose = TRUE)

fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'neg_probe',
                            norm_methods = 'standard',
                            library_size_norm = FALSE,
                            verbose = TRUE)

# new normalization method based on pearson correlations (Lause/Kobak et al. 2021)
# this normalized matrix is given the name 'pearson' using the update_slot param
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            scalefactor = 5000,
                            verbose = TRUE,
                            norm_methods = 'pearson_resid',
                            update_slot = 'pearson')

showGiottoExpression(fov_join)


# add statistics based on log normalized values for features rna and negative probes
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'rna')

fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'neg_probe')

# View cellular data (default is feat = 'rna')
showGiottoCellMetadata(fov_join)
# View feature data
showGiottoFeatMetadata(fov_join)


# 8. View Transcript Total Expression Distribution
# 8.1 Histogram of log normalized data
filterDistributions(fov_join,
                    detection = 'cells',
                    feat_type = 'rna',
                    expression_values = 'normalized',
                    method = 'sum',
                    nr_bins = 100,
                    save_param = list(base_height = 3,
                                      save_name = '5.1_rna_norm_total_hist'))


filterDistributions(fov_join,
                    detection = 'cell',
                    feat_type = 'neg_probe',
                    expression_values = 'normalized',
                    method = 'sum',
                    nr_bins = 20,
                    save_param = list(base_height = 3,
                                      save_name = '5.2_neg_norm_total_hist'))


# 8.2 Plot spatially as centroids
spatPlot2D(gobject = fov_join,
           cell_color = 'total_expr',
           color_as_factor = FALSE,
           show_image = TRUE,
           image_name = image_names,
           point_size = 0.9,
           point_alpha = 0.75,
           save_param = list(base_height = 2,
                             save_name = '5.3_color_centroids'))


# 8.3 Plot spatially as color-scaled polygons
spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     save_param = list(base_height = 2,
                                       save_name = '5.4_rna_color_polys'))

spatInSituPlotPoints(fov_join,
                     feat_type = 'neg_probe',
                     show_polygon = TRUE,
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     save_param = list(base_height = 2,
                                       save_name = '5.5_neg_color_polys'))



# 9. Dimension Reduction ------------------------------------------------------->>>>>>
# 9.1 Detect highly variable genes and generate PCA
# Detect highly variable genes using the pearson residuals method based on the ‘pearson’ expression matrix. These results will be returned as a new ‘hvf’ column in the ‘rna’ feature metadata.
# PCA generation will also be based on the ‘pearson’ matrix. Scaling and centering 
# of the PCA which is usually done by default will be skipped since the pearson 
# matrix is already scaled.
fov_join = calculateHVF(fov_join,
                        method = 'var_p_resid',
                        expression_values = 'pearson',
                        save_param = list(base_height = 5,
                                          save_name = '6.1_pearson_HVF'))

# print HVFs
gene_meta = fDataDT(fov_join)
gene_meta[hvf == 'yes', feat_ID]



fov_join = runPCA(fov_join,
                  scale_unit = FALSE,
                  center = FALSE,
                  expression_values = 'pearson')

# screeplot uses the generated PCA. No need to specify expr values
screePlot(fov_join, ncp = 20, save_param = list(save_name = '6.2_screeplot'))

plotPCA(fov_join,
        cell_color = 'nr_feats', # (from log norm statistics)
        color_as_factor = FALSE,
        point_size = 0.1,
        point_shape = 'no_border',
        save_param = list(save_name = '6.3_PCA'))


# 9.2 Run UMAP 
# Generate UMAP from PCA
fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    n_threads = 4)

plotUMAP(gobject = fov_join, save_param = list(save_name = '6.4_UMAP'))


# 9.3 Plot features on expression space
dimFeatPlot2D(gobject = fov_join,
              feat_type = 'rna',
              feats = c('MKI67', 'CD8A', 'CD4',
                        'COL1A1', 'MS4A1', 'MZB1'),
              expression_values = 'normalized',
              point_shape = 'no_border',
              point_size = 0.01,
              cow_n_col = 3,
              save_param = list(base_height = 5,
                                save_name = '6.5_UMAP_feats'))


# 10. Cluster ------------------------------------------------------------->>>>>>
# Verifiquei que era 13 clusters, entao usei o de 13 cores: pal13
# 10.1 Visualize clustering
fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)


fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.07,
                            n_iterations = 1000)


# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         cell_color_code = pal13,
         show_NN_network = TRUE,
         point_size = 2,
         save_param = list(save_name = '7.1_UMAP_leiden'))


# 10.2 Visualize clustering on expression and spatial space
# visualize UMAP and spatial results
spatDimPlot2D(gobject = fov_join,
              show_image = TRUE,
              image_name = image_names,
              cell_color = 'leiden_clus',
              cell_color_code = pal13,
              spat_point_size = 1,
              save_param = list(save_name = '7.2_spatdim_leiden'))


# 10.3 Map clustering spatially
spatInSituPlotPoints(fov_join,
                     feats = list('rna' = c('MMP2', 'VEGFA', 'IGF1R',
                                            'MKI67', 'EPCAM', 'MZB1')),
                     point_size = 0.15,
                     feats_color_code = viv10,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.01,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal13,
                     save_param = list(base_height = 5,
                                       save_name = '7.3_spatinsitu_leiden'))



# 11. Small Subset Visualization ----------------------------------------------->>>>>
# Acessando as coordenadas
coords <- fov_join@spatial_locs[[1]]$raw@coordinates
summary(coords)

#subset a Giotto object based on spatial locations
smallfov <- subsetGiottoLocs(fov_join,
                             x_max = 25070,   # Usando o valor máximo de x
                             x_min = 8669,    # Usando o valor mínimo de x
                             y_max = 158855,  # Usando o valor máximo de y
                             y_min = 155222)  # Usando o valor mínimo de y

#extract all genes observed in new object
smallfeats <- fDataDT(smallfov)[, feat_ID]

#plot all genes
spatInSituPlotPoints(smallfov,
                     feats = list(smallfeats),
                     point_size = 0.15,
                     polygon_line_size = 0.1,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     show_image = TRUE,
                     largeImage_name = 'fov002-composite',
                     show_legend = FALSE,
                     save_param = list(save_name = '8.1_smallfov_points'))

# plot only the polygon outlines
spatInSituPlotPoints(smallfov,
                     polygon_line_size = 0.1,
                     polygon_alpha = 0,
                     polygon_color = 'white',
                     show_polygon = TRUE,
                     show_image = TRUE,
                     largeImage_name = 'fov002-composite',
                     show_legend = FALSE,
                     save_param = list(save_name = '8.2_smallfov_poly'))

# plot polygons colorlabeled with leiden clusters
spatInSituPlotPoints(smallfov,
                     polygon_line_size = 0.1,
                     show_polygon = TRUE,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     show_image = TRUE,
                     largeImage_name = 'fov002-composite',
                     show_legend = FALSE,
                     save_param = list(save_name = '8.3_smallfov_leiden'))
