######################################################################
#
#          DOWNLOAD AND CLEAN ENVIRONMENTAL DATA
#                    
#
# Based on 
# Physalia course on SDMs (Oct 2021)
# by Marcia Barbosa
#
# Author: Daniel Romero Mujalli
# University of Greifswald
#
#######################################################################

#----------------------------------------------------------------------
# STEP 0:                 

# WARNING!
# some variables defined when preparing the biodiversity data are used
# in this script

#-----------------------------------------------------

# make sure that the working directory is appropriately set
getwd()

# create an "outputs" folder to store data
# Note that the new folder is created one level up
if(!file.exists("../outputs"))
    dir.create("../outputs")

# END OF STEP 0
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# STEP 1:           DOWNLOAD ENVIRONMENTAL VARIABLES

# We will use the sdmpredictors package to access different online
# datasets

pred_datasets <- sdmpredictors::list_datasets(terrestrial = TRUE
                                             ,marine      = TRUE
                                             )
# inspect the data available for download using sdmpredictors package
pred_datasets
names(pred_datasets)
# one can check the provided urls for more information
pred_datasets[ , 1:4]

# Remember to properly cite the actual data source
# not just the package
pred_datasets[, c("dataset_code", "citation")]

# check the sources of available datasets
pred_layers <- sdmpredictors::list_layers(datasets = pred_datasets)
unique(pred_layers$dataset_code)

# example of terrestrial variables dataset 
#(predictors as given by WorldClim)
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$name)

# example of marine variables dataset
# as given by MARSPEC
unique(pred_layers[pred_layers$dataset_code == "MARSPEC", ]$name)

# Let us choose one dataset (e.g., WorldClim) and one particular
# set of variables (e.g., altitude and the bioclimatic ones,
# which are in rows 1 to 20)
layers_choice <- unique(pred_layers[pred_layers$dataset_code=="WorldClim"
                                   , c("name", "layer_code")]
                       )
layers_choice
layers_choice <- layers_choice[1:20, ]
layers_choice

# create directory to store / fetch predictors data (map layers)
#options(sdmpredictors_datadir = "../outputs/sdmpredictors")
pred_dir <- "../outputs/sdmpredictors/"
if(!file.exists(pred_dir))
    dir.create(pred_dir)

# load / download selected layers
# We set rasterstack = FALSE, since TRUE gives error when there
# are layers with different extent
layers <- sdmpredictors::load_layers(layers_choice$layer_code
                                    ,rasterstack = FALSE
                                    ,datadir = pred_dir
                                    )
layers # a list of raster maps
# check how many elements in layers. This should match our choice above
length(layers)
# names of predictor variables 
names(layers)

# Before plotting some layers, we convert layers to SpatRaster class
# from terra package, which is much faster to process
layers <- lapply(layers, terra::rast)

terra::plot(layers[[1]], main = names(layers)[1])
terra::plot(layers[[5]], main = names(layers)[5])

# find out if layers have different extents or resolutions
# If there are different extents (Which typically does not happen
# with WorldClim, but may happen with other datasets, one will need
# to crop all layers to the minimum common extent before proceeding)
# For example, if the first layer has the smallest extent:
# layers <- lapply(layers, crop, extent(layers[[1]]))
# Resolutions:
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$cellsize_lonlat)

# Extents:
# we apply the ext function of terra package to the converted layers
sapply(layers, terra::ext)

# Once all layers have the same extent and resolution, one can combine them
# into a single multilayer raster object
# using the rast function of terra package
layers <- terra::rast(layers)
layers
terra::plot(layers)

# END OF STEP 1
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 2:             DELIMIT THE MODELLING REGION

# Let us convert species occurrences table to a spatial object.
# Similar to importing a delimited text file into a GIS, one
# needs to specify which columns of the table contain the spatial
# coordinates (point geometry geom) and what is the projection /
# coordinate reference system
# We use the vect function of terra package
# Refer to the documentation to understand the input values that go
# into each argument
occurrence_points <- terra::vect(occurrences
                                ,geom = c("decimalLongitude" # the x
                                         ,"decimalLatitude") # the y
                                ,crs = "+proj=longlat"
                                )
# plot the occurrence points
terra::plot(occurrence_points, col = "darkblue")
# add countries map
terra::plot(countries, border = "tan" # color of the borders
           ,lwd = 2, add = TRUE)

# In my case, I will exclude records of Central America
# thus,
toremove <- which(occurrences$decimalLatitude > 12)
dim(occurrences)
occurrences <- occurrences[-toremove, ]
dim(occurrences)

# save this cleaned file to disk
# write the cleaned data to *.csv file
write.csv(occurrences
         ,paste0(occurrence_dir, "/occurrences_", myspecies, "_cleaned.csv")
         ,row.names = FALSE
         )

occurrence_points <- terra::vect(occurrences
                                ,geom = c("decimalLongitude" # the x
                                         ,"decimalLatitude") # the y
                                ,crs = "+proj=longlat"
                                )

# Select the modelling region using the countries where this species
# has occurrence points (which means that the species was surveyed
# in those countries)

# from countries, select those in occurrence_points
countries_with_points <- countries[occurrence_points, ]

terra::plot(countries_with_points, border = "tan"
           ,col = scales::alpha(colour = "grey", alpha = 0.5)
           )
terra::plot(occurrence_points, col = "darkorange", cex = 0.3, add = TRUE)

# On the map, one should be able to assess whether the species has 
# been undersampled in some countries. It is important to compare
# the plot with other sources, e.g., iucnredlist.org, national atlases,
# data papers ...

# Using the default species, from the course: Spain and Portugal 
# seem to be the only evenly surveyed countries in this GBIF dataset
# also, only the mainland should be included in the modelling 
# region, as on islands species may be limited by factors 
# other than climate variables (e.g. dispersal)
# so, assuming we can't complement the dataset with occurrences 
# from other sources, let's select mainland Spain for modelling:

# Create a unique polygon identifier for "countries_with_points"
# and add it as labels to the map, to see which polygon(s) we
# want to select:
countries_with_points$my_id <- 1:length(countries_with_points)
terra::text(countries_with_points
           ,labels = countries_with_points$my_id
           ,col = "black"
           ,font = 2
           ,halo = TRUE
           )

# From the map, I want to exclude contries / polygons
# 3, 5, 14
#(which I will use to test the model when performing spatial projection)

# Select only the desired polygon(s) for the modelling region
# (from example data, we select polygons 1 and 2)
#selected_polygons <- c(1,2) # course example
polygons_to_remove <- c(3,5,14)
#model_region <- subset(countries_with_points
#                      ,countries_with_points$my_id %in% selected_polygons
#                      )
model_region <- subset(countries_with_points
                      ,!countries_with_points$my_id %in% polygons_to_remove
                      )

# check if desired polygons outline in green
terra::plot(model_region, border = "darkslategray4", lwd = 2, add = TRUE)

# Select the points in the model_region
occurrence_points <- occurrence_points[model_region, ]
terra::plot(occurrence_points, cex = 0.3, col = "darkslategray4", add = TRUE)

# if you can't select evenly surveyed countries 
# (e.g. if you're working with marine species), 
# you can delimit the modelling region as a buffer 
# of a given distance -- e.g. 1 geographic degree, 
# or 100 km, or the mean distance among points:

# takes time if there are many points; 
# may give "x$.self$finalize()" errors, 
# but check if object is successfully 
# created nonetheless:
#mean_dist <- mean(terra::distance(occurrence_points))
#mean_dist

##############################################################
# E.g., using a buffer of 100km
# Since the width argument of the terra::buffer function use meters
mean_dist <- 100 * 1000

# create predefined buffer using the buffer and aggregate functions
# of terra package
predefined_buffer <- terra::aggregate(terra::buffer(occurrence_points
                                                   ,width = mean_dist
                                                   )
                                     )
# take a look at the buffer (custom polygon)
terra::plot(countries[model_region,]
           ,col = scales::alpha(colour = "grey", alpha = 0.5)
           ,border = "grey", lwd = 2
           )
terra::plot(model_region, border = "darkgrey", lwd = 2, add = TRUE)
terra::plot(predefined_buffer, lwd = 2, border = "lightskyblue4", add = TRUE)
terra::plot(occurrence_points, col = "lightskyblue4", add = TRUE)


# if the buffer is to be combined with previously selected countries
# or polygons, the modelling region should be:
# (using the intersect function of terra package)
mod_region <- terra::intersect(predefined_buffer, model_region)
terra::plot(mod_region, border = "lightblue", lwd = 3, add = TRUE)

# if there are no previously selected polygons and one wants to
# use only the buffer, the modelling region should be:
#mod_region <- predefined_buffer

# IMPORTANT!!
# If one used a limited window of coordinates to download the occurrence
# data, one needs to intersect or crop with that too!
#mod_region <- terra::crop(model_region, terra::ext(my_window))
#terra::plot(mod_region, border = "green", lwd = 3, add = TRUE)

# aggregate modelling region into a single polygon:
mod_region <- terra::aggregate(mod_region)
terra::plot(mod_region,col = scales::alpha(colour = "grey", alpha = 0.5))

# Now import and cut (mask and trim) the variable maps to the
# extent of the modelling region
# load again the layers (that were downloaded previously) into a single
# raster object using rast function of terra package
layers <- terra::rast(list.files("../outputs/sdmpredictors"
                                ,pattern = "\\.tif$"
                                ,full.names = TRUE)
                     )
terra::plot(layers)

# crop(trim) layers to model region,
# then mask the values using the model region as template
model_region <- mod_region

# end of using buffer
################################################################

layers_cut <- terra::mask(terra::crop(layers, model_region)
                         ,mask = model_region
                         )
terra::plot(layers_cut)

# names of the environmental layers (predictors)
names(layers_cut)

# Given the case,
# let us simplify the names by removing the tail "_lonlat"
#names(layers_cut) <- sub(pattern = "_lonlat", replacement = ""
#                        ,x = names(layers_cut)
#                        )

# plot first predictor
terra::plot(layers_cut[[1]])
# add countries borders
terra::plot(countries, border = "darkgrey", add = TRUE)
# highlight selected modelling region
terra::plot(model_region, add = TRUE)
# add occurrence points
terra::plot(occurrence_points, pch = 20, cex = 0.1, add = TRUE)
# Make sure that everything overlaps correctly

# END OF STEP 2
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 3:          SET THE APPROPRIATE SPATIAL RESOLUTION

# closely inspect the species data points vs the size of the variables'
# pixels
# xlim and ylim range depends on the coordinate values!
# For the default species:
#terra::plot(layers_cut[[1]], xlim = c(-5, -3), ylim = c(38, 40))
#terra::points(occurrence_points, cex = 0.5)
# in case of latin america map:
terra::plot(layers_cut[[1]], xlim = c(-78, -76), ylim = c(-2, 0))
terra::points(occurrence_points, cex = 0.5)

# plot within different x/y limits if necessary to see if your 
# presence point resolution matches pixel resolution 
# (i.e., if you don't have many evenly spaced points with 
# pixels in between -- see lecture "04_predictor_variables.pdf")
# notice that, in the example data, pixels have approximately 
# the same spatial resolution as the presence points, but 
# they don't match exactly (i.e., points are not at the centroids 
# of these pixels); ideally, you should find the original grid 
# over which (most of) these presences were sampled, and extract 
# the raster values to that grid

# if necessary (which is not the case for the example data), 
# you can aggregate the layers, to e.g. a 5-times coarser resolution 
# (choose the 'fact' value that best matches your presence data 
# resolution to your variables' resolution):
# Resolution before
terra::res(layers_cut) # returns x and y resolution
layers_aggr <- terra::aggregate(layers_cut
                               ,fact = 5   # aggregation factor
                               ,fun = mean # function used to aggregate values
                               )
# Resolution after aggregation
terra::res(layers_aggr)

terra::plot(layers_aggr[[1]]
           ,xlim = range(occurrence_points$decimalLongitude)
           ,ylim = range(occurrence_points$decimalLatitude)
           )
terra::points(occurrence_points, pch = ".")

# run the command below only if you did need to aggregate the layers:
#layers_cut <- layers_aggr

# save the cut layers to a folder on disk for later use:

outdir_layers_cut <- paste0("../outputs/sdmpredictors/layers_cut_"
                           ,myspecies
                           )
if(!file.exists(outdir_layers_cut))
    dir.create(outdir_layers_cut)
terra::writeRaster(layers_cut
                  ,filename = paste0(outdir_layers_cut, "/layers_cut.tif")
                  ,overwrite = TRUE
                  )

# now we create a dataframe of the species occurrence data gridded to
# the resolution of the raster variables
# i.e., one row per pixel with the values of the predictor variables
# and the presence/absence of species records
head(occurrences)
#?gridRecords from fuzzySim package
presence_data_coord <- occurrences[, c("decimalLongitude", "decimalLatitude")]
gridded_data <- fuzzySim::gridRecords(rst = layers_cut # predictors rasobject
                                     ,pres.coords = presence_data_coord
                                     )
# inpect the data
head(gridded_data)
dim(gridded_data)
# rows should be the same as that for the predictor layers and occurrence 
# data after removing nas
sum(!is.na(terra::values(layers_cut[[1]])))

# check column names of dataframe that includes occurrence data and
# the predictor variables, with corresponding coordinate values (grid)
names(gridded_data)

# plot the gridded records
terra::plot(layers_cut[[1]])
# plot the absences (pixels without presence records)
terra::points(gridded_data[gridded_data[, "presence"] == 0, c ("x", "y")]
             ,col = "red", cex = 0.1
             )
# plot the presences (pixels with presence records)
terra::points(gridded_data[gridded_data[, "presence"] == 1, c ("x", "y")]
             ,col = "blue", cex = 0.1
             )

# plot a narrower coordinate range to take a closer look
# see the coordinate range on the plot axes, and pick some limits
# within which to look closer
# course default species:
#terra::plot(layers_cut[[1]], xlim = c(-5, -3), ylim = c(38, 40))
#terra::points(occurrence_points, cex = 0.5)

# new world monkeys:
terra::plot(layers_cut[[1]], xlim = c(-75, -65), ylim = c(-10, 0))
terra::points(occurrence_points, cex = 0.5)

# absence points
terra::points(gridded_data[gridded_data[, "presence"] == 0
                          , c("x", "y")
                          ]
             , col = "red"
             , pch = 1, cex = 0.5
             )
# presence points
terra::points(gridded_data[gridded_data[, "presence"] == 1
                          , c("x", "y")
                          ]
             , col = "blue"
             , pch = 20, cex = 0.5
             )

# It is possible to see that there is one data point per pixel
#----------

# save the modelling dataframe to a *.csv file on disk for future use
write.csv(gridded_data
         ,paste0("../outputs/dat_gridded_", myspecies, ".csv")
         ,row.names = FALSE
         )

#----------------------------------------------------------------------
# END OF DOWNLOAD AND CLEAN ENVIRONMENTAL DATA