######################################################################
#
#                           MODEL PROJECTION
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

# Required packages
pkgs <- c("rgbif", "sdmpredictors", "fuzzySim", "terra"
        , "modEvA", "dismo", "maxnet", "gam", "randomForest", "gbm"
        , "corrplot", "ecospat", "blockCV", "raster", "plotmo"
         )

# check if pkgs are available
my_packages <- rownames(installed.packages())

for (pkg in pkgs)
{
    if(!pkg %in% my_packages)
        print(pkg)
}

# install any missing package

# make sure that the working directory is appropriately set
getwd()

# create an "outputs" folder to store data
# Note that the new folder is created one level up
if(!file.exists("../outputs"))
    dir.create("../outputs")

# WARNING! We are using informaiton from previous steps / scripts!

# define a colour palette and colour breaks (intensity), so that all maps of model 
# predictions show a comparable gradient
# HCL Hue Chroma Luminance
mypalette <- hcl.colors( n = 10               # the number of colors
                       , palette = "viridis" # desired palette (check ?hcl.colors)
                       )
colorbreaks <- seq(from = 0, to = 1, by = 0.1)

#### IMPORT DATA FROM PREVIOUS STEPS

# Select the species
myspecies <- "Callitrichidae" # family of marmosets and tamarins
#myspecies <- "Pleurodeles waltl" # from the course

# set directory where models are stored
outdir_models <- paste0("../outputs/models_",myspecies)
# read models from data using readRDS base function to read info on r objects
# uncomment as necessary
#bioclim_mod <- readRDS(paste0(outdir_models, "/bioclim_mod.rds"))
#domain_mod <- readRDS(paste0(outdir_models, "/domain_mod.rds"))
#maxent_mod <- readRDS(paste0(outdir_models,"/maxent_mod.rds"))
maxnet_mod <- readRDS(paste0(outdir_models,"/maxnet_mod.rds"))
#glm_mod <- readRDS(paste0(outdir_models,"/glm_mod.rds"))
#gam_mod <- readRDS(paste0(outdir_models,"/gam_mod.rds"))
#rf_mod <- readRDS(paste0(outdir_models, "/rf_mod.rds"))
#gbm_mod <- readRDS(paste0(outdir_models, "/gbm_mod.rds"))

dat <- read.csv(paste0("../outputs/dat+preds_",myspecies,".csv"))
head(dat)

# END OF STEP 0
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# STEP 1: RESPONSE CURVES

# before projecting the models, it is recommended to check the response 
# curves to make sure that they make ecological sense
# Keep in mind that individual response curves may be affected by 
# interactions with other variables. Our visualization capability is
# limited, however useful to have an overview

# We use the response function from dismo package for 'dismo' models:
# bioclim, domain, maxent; plot for maxnet model; plotmo 
# function from plotmo package for the other models
# uncomment as necessary
# Single variable response curves (x: the model object)
#dismo::response(x = bioclim_mod) # bioclim
#dismo::response(x = domain_mod)  # domain
#dismo::response(x = maxent_mod)  # maxent

plot(x = maxnet_mod, type = "cloglog")

# the plotmo function plots the model's response when varying one or 
# two predictors while holding the other predictors constant
# check the documentation ?plotmo::plotmo for more information
plotmo::plotmo(object = glm_mod # the model object
              ,trace = 1        # summary trace? 1:TRUE
              ,all1 = TRUE      # plot all predictors?
              ,all2 = TRUE      # plot all pairs of predictors?
              )
plotmo::plotmo(object = gam_mod # the model object
              ,trace = 1        # summary trace? 1:TRUE
              ,all1 = TRUE      # plot all predictors?
              ,all2 = TRUE      # plot all pairs of predictors?
              )
# follow the same logic to produce the response curves from the
# other models

# END OF STEP 1:
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 2: SPATIAL PROJECTION

### FIRST, IMPORT OCCURRENCE RECORDS, PREDICTORS MAPS, AND BACKGROUND
### COUNTRY BORDERS
# This will help us creating the maps

# import the cleaned occurrence records
occurrences <- read.csv(paste0("../outputs/species_occurrences/"
                              ,"occurrences_"
                              ,myspecies
                              ,"_cleaned.csv")
                       )

# convert occurrence records dataframe to spatial vector
occurrence_points <- terra::vect(x = occurrences    # dataframe / file..
                                ,geom = c("decimalLongitude"
                                         ,"decimalLatitude"
                                         )
                                ,crs = "+proj=longlat"
                                )
# visual check
terra::plot(occurrence_points)

# import the countries map
countries <- terra::vect(x = "../data/countries/world_countries.shp")
terra::plot(x = countries
           ,xlim = range(occurrences$decimalLongitude)
           ,ylim = range(occurrences$decimalLatitude)
           )
terra::points(x = occurrence_points, col = "blue")

# Import the complete (global coverage) predictor / variable layers
# $: indicate that the pattern is at end of the string
# the character "\\" ensures to match extension files ending with
# the pattern ".tif" (in this case), in case unwanted files with
# similar ending exist (e.g., ".atif" <- excluded)
layers <- terra::rast(list.files(path = "../outputs/sdmpredictors/"
                                ,pattern = "\\.tif$" # find pattern
                                ,full.names = TRUE # names with path
                                )
                     )
# check that varible names make sense!
names(layers)
# adjust if necessary
names(layers) <- sub(pattern = "_lonlat", replacement = "", x = names(layers))

# IMPORTANT!!
# We may need to select the same environmental variables as those 
# used for the model training! Remember we selected a subset of 
# no correlated variables (see previous steps)
# Otherwise, projected predictions will not work
# We use paste0, instead of paste function, such that strings are
# concatenated without spaces
selected_variables <- read.csv(paste0("../outputs/variables_sel_"
                                    ,myspecies
                                    ,".csv"
                                    )
                              )
selected_variables <- as.vector(selected_variables[,1]) # as vector
# select only the layers corresponding to the selected variables,
# use for modelling training
layers <- layers[[selected_variables]]

# visual check
terra::plot(layers[[1:4]]) # plot the first 4 predictor layers

### SECOND, DEFINE THE PROJECTION EXTENT OF THE REGION OF INTEREST

# select the projection extent: this is, the region where to project
# model predictions
# Here one can include those regions with records that were excluded
# when training the model
# This is the one I used for the example with the new world monkeys
proj_extent <- terra::ext(x = occurrence_points) # x: SpatRaster
# alternatively, one can select a specific spatial window of interest
#proj_extent <- terra::ext(c(-15, 10, 30, 45)) # course example

# crop layers (x) to the extent of interest given by (y)
layers_proj <- terra::crop(x = layers, y = proj_extent)
# visual check
terra::plot(x = layers_proj, col = mypalette)

### THIRD, PREDICT WITH EACH MODEL THE EXPECTED DISTRIBUTION ON THE
### REGION OF INTEREST
# Note that, for some, we perform the prediction and convert the 
# outcome to SpatRaster object
# bioclim
bioclim_proj <- dismo::predict(object = layers_proj # a raster object
                              ,model = bioclim_mod # model object
                              )
# domain
domain_proj <- dismo::predict(object = layers_proj
                             ,model = domain_mod
                             )
# maxent
maxent_proj <- terra::rast(dismo::predict(raster::stack(layers_proj)
                                         ,maxent_mod
                                         )
                          )
# maxnet
maxnet_proj <- terra::rast(raster::predict(raster::stack(layers_proj)
                                          ,model = maxnet_mod
                                          ,type = "cloglog"
                                          )
                          )
# glm
glm_proj <- raster::predict(object = layers_proj
                           ,model = glm_mod
                           ,type = "response"
                           )
# gam
gam_proj <- raster::predict(object = layers_proj
                           ,model = gam_mod
                           ,type = "response"
                           )
# rf
rf_proj <- raster::predict(object = layers_proj
                          ,model = rf_mod
                          ,type = "prob")["X1"]
# gbm
gbm_proj <- reaster::predict(object = layers_proj
                            ,model = gbm_mod
                            ,type = "response"
                            )

# Convert probability to favourability (only presence-absence models)
# glm
my_prevalence <- modEvA::prevalence(mod = glm_mod)
glm_proj_fav <- fuzzySim::Fav(pred = glm_proj
                             ,sample.preval = my_prevalence
                             )
# gam
my_prevalence <- modEvA::prevalence(mod = gam_mod)
gam_proj_fav <- fuzzySim::Fav(pred = gam_proj
                             ,sample.preval = my_prevalence
                             )
# Follow the same logic to get the favourability of the other models 

# FOURTH, MAP THE EXTRAPOLATED PROJECTIONS
# One can add country borders (optional)

# bioclim
terra::plot(bioclim_proj, col = mypalette, breaks = colorbreaks
           ,main = "Bioclim"
           )
# maxent
terra::plot(maxent_proj, col = mypalette, breaks = colorbreaks
           ,main = "Maxent"
           )
# maxnet
terra::plot(maxnet_proj, col = mypalette, breaks = colorbreaks
           ,main = "Maxnet"
           )
# glm
terra::plot(glm_proj_fav, col = mypalette, breaks = colorbreaks
           ,main = "GLM_FAV"
           )
# gam
terra::plot(gam_proj_fav, col = mypalette, breaks = colorbreaks
           ,main = "GAM_FAV"
           )
# follow the same logic to plot the projections for the other
# modelling methods

# add contries (optional)
terra::plot(countries, xlim = range(occurrences$decimalLongitude)
            ,ylim = range(occurrences$decimalLatitude), add = TRUE)

# END OF STEP 2:
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 3: ANALYSIS OF ENVIRONMENTAL DISSIMILARITY

# we use the function mess from package dismo (Multivariate 
# Environmental Similarity Surface)
# Compute multivariate environmental similarity surfaces (MESS), as
# described by Elith et al., 2010

# import the layers used for training the model
path_layers_cut <- paste0("../outputs/sdmpredictors/layers_cut_"
                         ,myspecies
                         )
layers_cut <- terra::rast(list.files(path = path_layers_cut
                                    ,pattern = "\\.tif$"
                                    ,full.names = TRUE
                                    )
                         )
# visual check
terra::plot(x = layers_cut[[1:4]], col = mypalette)

# compute dissimilarity of projection layers relative to modelled
# layers
# x: raster object (projected region of interest)
# y: matrix or dataframe containing the reference values (such that
# values ranges can be extracted from those used for the model 
# training)
mess_proj <- dismo::mess(x = raster::stack(layers_proj)
                        ,v = na.omit(terra::values(layers_cut))
                        )
mess_proj

# visual inspection
# negative values indicate environmental dissimilarity from the 
# reference region
terra::plot(x = mess_proj, col = heat.colors(10), main = "MESS")

# Let us have a better look at those regions with negative mess 
# values, since this means that predictions on these areas are
# less reliable
mess_proj_unreliable <- mess_proj
mess_proj_unreliable[mess_proj_unreliable > 0] <- NA
# visual inspection
terra::plot(mess_proj_unreliable, col = heat.colors(10), main = "MESS")
terra::plot(terra::crop(countries, mess_proj_unreliable), add = TRUE)

# END OF STEP 3:
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 4: MODEL ENSEMBLES

# CURRENTLY, I AM NOT COVERING THIS TOPIC
# YOU MAY VISIT THE COURSE MATERIALS FOR MOR INFORMATION HERE

# END OF STEP 4:
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 5: TIME PROJECTION

# There are many possiblities to download data about future scenarios.
# What is important is to find future data information for the predictor
# variables used for modelling training. This is so, because the model
# can provide an answer based on future observations of the variables
# already screened (learned) by the model. E.g., if we train the model
# to predict species distribution based on temperature, we need to get
# data on future temperature scenarios to evaluate model projections

# Here we use the 'sdmpredictors' package which was previously used to
# extract information on current climate. The same package provides
# the possibility to extract data on future climate as well.

# get data
future_layers <- sdmpredictors::list_layers_future()
# inspect data
head(future_layers)
names(future_layers)
# check available datasets for future climate data
unique(future_layers$dataset_code)

# extract worldclim future data
# we use the subset base function
future_worldclim <- subset(x = future_layers  # object to be subsetted
                          ,dataset_code == "WorldClim" # logical expr
                          )
# extract bio-oracle future data
future_biooracle <- subset(x = future_layers  # object to be subsetted
                          ,dataset_code == "Bio-ORACLE" # logical expr
                          )

# check what layers/scenarios of future climate are available
# we use the function unique of base R to extract unique values
# for the requested column names
# worldclim:
unique(future_worldclim[,c("model", "scenario", "year", "version")])

# biooracle:
unique(future_biooracle[,c("model", "scenario", "year", "version")])
# whenever there are more than one version, it is recommended to select
# the most recent one

# select the most recent version in biooracle dataset
future_biooracle <- subset(future_biooracle
                          ,version == max(future_biooracle$version)
                          )
# check that the selection was correctly done
unique(future_biooracle[, c("model", "scenario", "year", "version")])

### SELECT CLIMATE MODEL(S), EMISSION SCENARIO(S) AND YEAR(S) TO PROJECT
### THE MODEL

# check available scenarios
unique(future_worldclim[, c("model", "scenario", "year")])

unique(subset(future_worldclim, model == "CCSM4" 
              & scenario == "rcp26"
              & year == 2050
              )
      )

# variable selection
variables_2050_rcp26 <- subset(x = future_worldclim
                              ,subset = model    == "CCSM4"
                                      & scenario == "rcp26"
                                      & year     == 2050
                              )

# layers of selected variables
data_dir <- "../outputs/sdmpredictors/future"
if(!file.exists(data_dir))
    dir.create(data_dir)

layers_2050_rcp26 <- sdmpredictors::load_layers(layercodes = variables_2050_rcp26
                                               ,rasterstack = FALSE
                                               ,datadir = data_dir
                                               )
# combine into rasterstack
layers_2050_rcp26 <- raster::stack(layers_2050_rcp26)

# Before using these layers, just as it was done for spatial projection,
# it is important to make sure that variable names are the same, for
# the future layer stack and for the models

names(layers) # names used when creating the model
names(layers_2050_rcp26) # names in layers of future scenarios

# It is possible to see that they differ in the ending strings

# adjust the names in future layers, so that they match the corresponding
# names of current climate data layers
# we use the gsub base function to find and replace given pattern
names(layers_2050_rcp26) <- gsub(pattern = "_cc26_2050"
                                ,replacement = "_lonlat"
                                ,x = names(layers_2050_rcp26)
                                )
# check names
names(layers_2050_rcp26) # should match now

# check in what variables layers of current climate differ to those
# of future climate
# we use the setdiff function of base R to identify the value of the
# element(s) that differ between the two vectors of names
# we ask what is missing in y that exists in x?
setdiff(x = names(layers_cut), y = names(layers_2050_rcp26))
# It is possible to see that altitude data "WC_alt_lonlat" is missing
# for the future climate. However, we can assume that altitude values
# will not significantly change for the projected time window, and use
# the exact same values used for current climate

# add altitude values to the stack of layers of future climate
layers_2050_rcp26 <- raster::stack(raster::raster(layers[["WC_alt_lonlat"]])
                                  , layers_2050_rcp26
                                  )
# check that all variables used for the modelling are present in the
# stack of layers of future climate
setdiff(x = names(layers_cut), y = names(layers_2050_rcp26))
# visual inspection
terra::plot(layers_2050_rcp26, col = mypalette)

### CROP LAYERS TO REGIONS OF INTEREST AND PROJECT MODEL IN TIME

# crop layers. We use the crop function of raster package
# x: raster object
# y: extent object

# we use the first layer of the region of interest to get the
# extent values that can be used to crop the layers of future
# scenarios.
# Only selected one layer is enough, and computationally more
# efficient than taking the full object
raster_obj_with_extent <- raster::raster(layers_proj[[1]])
layers_2050_rcp26 <- raster::crop(x = layers_2050_rcp26
                                 ,y = raster_obj_with_extent
                                 )

# Now we can proceed to project in time using the same procedure used
# for the spatial projection task, above
# Here only one example is provided. For time projections using the
# other models, please refer to the spatial projection performed above
gam_2050_rcp26 <- dismo::predict(object = layers_2050_rcp26
                                ,model = gam_mod
                                ,type = "response"
                                )
# visualization
# let us plot present vs future projections
par(mfrow = c(2,1))
raster::plot(raster::raster(gam_proj)
            ,col = mypalette
            ,main = "GAM present"
            )
raster::plot(gam_2050_rcp26
            ,col = mypalette
            ,main = "GAM 2050 (RCP 2.6)"
            )
# END OF STEP 5:
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# END OF MODEL PROJECTION