######################################################################
#
#                           MODELLING METHODS
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
# STEP 0:                 INSTALLING PACKAGES

# Required packages
pkgs <- c("rgbif", "sdmpredictors", "fuzzySim", "terra"
        , "modEvA", "dismo", "maxnet", "gam", "randomForest", "gbm"
        , "corrplot", "raster", "sp", "Rcpp"
         )

# check if pkgs are available
my_packages <- rownames(installed.packages())

for (pkg in pkgs)
{
    if(!pkg %in% my_packages)
        print(paste("missing_package: ",pkg))
}

# proceed to install any missing package
# LINUX USERS
# the dismo package requires
# sudo apt-get install libgdal-dev
# and the "rJava" package (see the bottom of this file)
#----------------------------------------------------------------------

# make sure that the working directory is appropriately set
getwd()

# create an "outputs" folder to store data
# Note that the new folder is created one level up
if(!file.exists("../outputs"))
    dir.create("../outputs")

# END OF STEP 0
#----------------------------------------------------------------------
# STEP 1:       IMPORT DATA

# select the species
# e.g.,
#myspecies <- "Pleurodeles waltl" # the species of interest (course)
myspecies <- "Callitrichidae" # family of marmosets and tamarins

# import the data prepared from the previous scripts
# 1) dataframe gridded data of occurrences and predictors
filepath <- paste0("../outputs/dat_gridded_", myspecies, ".csv")
dat <- read.csv(filepath)
head(dat)
# 2) environmental layers cut and merged into a single layer object
filepath <- paste0("../outputs/sdmpredictors/layers_cut_",myspecies,"/layers_cut.tif")
layers_cut <- terra::rast(filepath)
# visual check
terra::plot(layers_cut)

# define a colour palette and colour breaks (intensity), so that all maps of model 
# predictions show a comparable gradient
# HCL stands for Hue Chroma Luminance
mypalette <- hcl.colors( n = 10               # the number of colors
                       , palette = "viridis" # desired palette (check ?hcl.colors)
                       )
colorbreaks <- seq(from = 0, to = 1, by = 0.1)

# For the case of modelling presence data only, we make a version of layers_cut with
# pixel values only at the locations of the presence points
names(dat) # check the name of the coordinates information
# create mask
presence_centroids <- terra::vect(dat[dat$presence == 1, ], geom = c("x","y")
                                 ,crs = "+proj=longlat")
terra::plot(presence_centroids) # visual check

# create layers_cut for presence only data using the mask presence_centroids
# Here the values of the rasterobject x will become NAs except for those overlapping
# with the mask
layers_cut_presonly <- terra::mask(x = layers_cut, mask = presence_centroids)
terra::plot(layers_cut_presonly, col = mypalette)
terra::plot(layers_cut_presonly[[1]])

# Some of the modelling methods required specific format of Spatial/Raster data
# here we convert terra format maps (Spatvector and SpatRaster) to SpatialPoints or
# RasterStack, as required by the modelling methods
presence_centroids_spatial <- as(presence_centroids, "Spatial")
# we use the stack function of the package raster to convert to stack raster object
# A RasterStack is a collection of RasterLayer objects with the same spatial extent 
# and resolution. A RasterStack can be created from RasterLayer objects, or from 
# raster files, or both. It can also be created from a SpatialPixelsDataFrame or a 
# SpatialGridDataFrame object.
layers_cut_stack <- raster::stack(layers_cut)
terra::plot(layers_cut_stack, col = mypalette)
layers_cut_presonly_stack <- raster::stack(layers_cut_presonly)
terra::plot(layers_cut_presonly_stack, col = mypalette)

# END OF STEP 1
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 2:            SELECTION OF VARIABLES / PREDICTORS

# check column names in dataframe
names(dat)

# let us simplify the names by removing the "_lonlat" tail
#names(dat) <- sub(pattern = "_lonlat", replacement = "", x = names(dat))

# we use pattern matching (function grep ) to select the predictor variables
selection <- grep(pattern = "WC_alt|WC_bio", x = names(dat)) # | -> "or" logic
variables <- names(dat)[selection]
print(variables) # check if variables predictors were properly selected

# In the case, we want to exclude some variables, one could
# toremove <- names(dat)[grep(pattern = "bio3|bio14|bio15", x = names(dat))]
# we use the setdiff function that excludes from x any match with y
# selection <- setdiff(x = vars, y = toremove)

# Keep in mind that many methods are affected by having too many variables
# and / or high correlations among them
# Let us check if there are highly correlated values (see the arbitrary 
# threshold) and exclude them from the analysis. This is, if there are
# two highly correlated values, we keep only one of them for the modelling
# We use the AIC as the criterion to decide which variable to exclude when
# a high correlation is found
# AIC: akaike information criterion, an estimator of prediction error
variables_sel <- fuzzySim::corSelect(data = dat
                          ,sp.cols = "presence"  # species / response variable colname
                          ,var.cols = variables  # predictors variables
                          ,cor.thresh = 0.8      # arbitrary corr threshold
                          ,select = "AIC"        # criterion for excluding correlated vars
                          )
variables_sel   # list with results from the correlation analysis
variables_sel <- variables_sel$selected.vars # keep only the selected variable names
variables_sel

# write selected variables to *.csv file
write.csv(variables_sel, paste0("../outputs/variables_sel_",myspecies,".csv"),
          row.names = FALSE
          )

# to read the selected variables (we apply transpose function and read as vector):
variables_sel <- t(read.csv(paste0("../outputs/variables_sel_",myspecies,".csv")))
variables_sel <- as.vector(variables_sel)

# END OF STEP 2
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 1:              PRESENCE ONLY MODELLING: BIOCLIM

# When calling the functions of modelling methods, consider the following
# x: predictors
# p: occurrence data

# we use the bioclim function of the dismo package
# this is a model based on presence only (po) data and its predictions are not that
# "powerful" as those from other modelling methods
# x: raster object; y: spatialPoints object or two column matrix
bioclim_mod_po <- dismo::bioclim(x = layers_cut_presonly_stack[[variables_sel]]
                                , p = presence_centroids_spatial
                                )

# We can test and confirm that we get the same model performance whether or not 
# we provide pixels without presence records, showing that Bioclim is a true
# presence only method (it uses only the presence pixels)
bioclim_mod <- dismo::bioclim(x = layers_cut_stack[[variables_sel]]
                             ,p = presence_centroids_spatial
                             )
# we use all.equal() native function to test whether two given r objects are
# equal
all.equal(target = bioclim_mod_po, current = bioclim_mod)

# Plot model response: it includes the first two variables, but the response
# includes all variables in the model
par(mar = c(4,4,2,1))
terra::plot(bioclim_mod)

# compute and map Bioclim predictions accross the study region
# and overlay the original presences
bioclim_pred <- dismo::predict(layers_cut, bioclim_mod)

terra::plot(bioclim_pred, col = mypalette, breaks = colorbreaks, main = "Bioclim")
terra::plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# END OF MODEL 1: BIOCLIM
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 2:              PRESENCE ONLY MODELLING: DOMAIN

# presence only po
domain_mod_po <- dismo::domain(x = layers_cut_presonly_stack[[variables_sel]]
                              ,p = presence_centroids_spatial
                              )
# pixels without presence
domain_mod <- dismo::domain(x = layers_cut_stack[[variables_sel]]
                           ,p = presence_centroids_spatial
                           )
all.equal(target = domain_mod_po, current = domain_mod)
# domain also returns the same model independently of using presence only 
# pixels of the predictor layers or not (i.e., it is a presence only model)

# END OF MODEL 2: DOMAIN
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 3:          PRESENCE BACKGROUND MODELLING: MAXENT

# WARNING! IT CAN TAKE LONG TIME TO RUN!

# here we use the maxent function in dismo package
# dimo maxent function uses the maxent.jar file
# read documentation on how to download and where to save the file
# On Linux you will need to setup java and install the rJava package
# Instructions from the following link: 
# https://datawookie.dev/blog/2018/02/installing-rjava-on-ubuntu/
# are found at the end of this file

# The following call tries using presence only data and therefore results in an
# error. Maxent needs background data without presence points in order to
# train the model
#maxent_mod_po <- dismo::maxent(x = layers_cut_presonly_stack[[variables_sel]]
#                              ,p = presence_centroids_spatial
#                              )

# Maxent also takes the dataframe instead of the raster maps + presence points:
maxent_mod <- dismo::maxent(x = dat[, variables_sel]    # dataframe instead of rasterslack
                           ,p = dat$presence            # occurrence data column in df
                           )           

# let us examine the maxent model object
str(maxent_mod) # there are presence and absence slots
head(maxent_mod@presence)
head(maxent_mod@absence) # it creates absences from the pixels wihtout presence points
nrow(maxent_mod@presence)
nrow(maxent_mod@absence)

# compute and map maxent predictions
maxent_pred <- dismo::predict(object = maxent_mod       # a fitted model
                             ,x = layers_cut_stack      # raster object / dataframe
                             )
# plot predictions
terra::plot(maxent_pred, col = mypalette, breaks = colorbreaks, main = "MaxEnt")
terra::plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# END OF MODEL 3: MAXENT
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 4:          PRESENCE BACKGROUND MODELLING: MAXNET

# WARNING! IT CAN TAKE LONG TIME TO RUN!

# the maxnet function in maxnet package uses glmnet rather than the java version of
# maxent. Thus, it runs independently from java, in contrast to the method above
# Maxnet works on a matrix or dataframe! rather than raster maps + presence points
# As demonstrated above, the provided data needs to contain presence and no presence
# points. Otherwise it will not work
# f: a formula to determine the features to be used (below we use the deafault)
myformula <- maxnet::maxnet.formula(p = dat$presence,data = dat[,variables_sel])

maxnet_mod <- maxnet::maxnet(p = dat$presence # a vector: 1 (presence), 0 (background)
                            ,data = dat[,variables_sel] # the predictors
                            ,f = myformula
                            )

# compute and map maxnet predictions
# we use the predict function of raster package
maxnet_pred <- raster::predict(object = layers_cut_stack # rasterobj where to predict on
                              ,model = maxnet_mod   # model object with predict attribute
                              ,type = "cloglog"     # complementary log-log type of output
                              )

terra::plot(maxnet_pred, col = mypalette, breaks = colorbreaks, main = "MaxNet")
terra::plot(presence_centroids, pch = ".", add = TRUE)

# END OF MODEL 4: MAXNET
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 5:            PRESENCE ABSENCE MODELLING: GLM

# we use glm function in package stats

# let us first define the formula required by the glm function
# collapse argument is the value that binds the string inputs
glm_formula <- as.formula(paste("presence~", paste(variables_sel, collapse = " + ")))
glm_formula

# we fit the model using a binomial family because our binary response variable
glm_mod <- glm(formula = glm_formula, family = binomial, data = dat)
summary(glm_mod)

# we use the predict function in raster package so that we can pass on a raster object
# where to perform the predictions. Type "response" transforms the linear prediction
# with the appropriate link function (e.g., logistic regression), 
# yielding results as "probability" values
glm_pred <- raster::predict(object = layers_cut, model = glm_mod, type = "response")

terra::plot(glm_pred, col = mypalette, breaks = colorbreaks, main = "GLM")
terra::plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# The predictions of GLM and other presence-absence models are of
# presence "probability", which considers the effects of the predictor variables
# and the baseline prevalence (proportion of presences) of the species in the
# model training data. Probabilities are generally low for rare species, and
# generally higher for widespread species

# We can convert probability to prevalence-independent favourability:
# we use the Fav function in fuzzySim package
# to the argument sample.preval we give the prevalence as given in our glm model
# The prevalence is computed using the prevalence function in ModEvA package
glm_mod_prevalence <- modEvA::prevalence(model = glm_mod)

glm_fav <- fuzzySim::Fav(pred = glm_pred, sample.preval = glm_mod_prevalence)

# plot glm predictions converted into favourability
terra::plot(glm_fav, col = mypalette, breaks = colorbreaks, main = "GLM_FAVOURABILITY")
terra::plot(presence_centroids, pch = ".", add = TRUE)

# END OF MODEL 5: GLM
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 6:           PRESENCE ABSENCE MODELLING: GAM

# using gam function in gam package
# Generalized Additive Models GAM
# The procedure is similar to that for GLM

# formula
# The "s" character is used to specify that each variable needs to undergo
# smoothing splines (function s in package gam)
gam_formula <- as.formula(paste("presence~"
                               ,paste0("gam::s(",variables_sel,")"
                                      ,collapse = " + "
                                      )
                                )
                         )
gam_formula

# fit model
gam_mod <- gam::gam(formula = gam_formula, family = binomial, data = dat)
summary(gam_mod)

# predictions
gam_pred <- raster::predict(object = layers_cut, model = gam_mod, type = "response")

# plot predictions
terra::plot(gam_pred, col = mypalette, breaks = colorbreaks, main = "GAM")
terra::plot(presence_centroids, col = "red", pch = ".", add = TRUE)

# convert probability to prevalence-independent favourability:
# we follow the same procedure as for GLM above
gam_mod_pevalence <- modEvA::prevalence(model = gam_mod)
gam_fav <- fuzzySim::Fav(pred = gam_pred, sample.preval = gam_mod_pevalence)

terra::plot(gam_fav, col = mypalette, breaks = colorbreaks, main = "GAM_FAVOURABILITY")
terra::plot(presence_centroids, col = "red", pch = ".", add = TRUE)

# END OF MODEL 6: GAM
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 7:        PRESENCE ABSENCE MODELLING: RANDOM FOREST

# WARNING! IT CAN TAKE LONG TIME TO RUN!

# we use the randomForest function in randomForest package

# first, the response variable (presence column in dat) needs to be passed
# to the function as factor, such that the function performs classification
# and not regression
dat$presence_factor <- as.factor(dat$presence)

rf_formula <- as.formula(paste("presence_factor~",paste(variables_sel
                                                       ,collapse = " + ")))
rf_formula

# check documentation ?randomForest::randomForest to understand parameter values
# and prepare a better model
rf_mod <- randomForest::randomForest(formula = rf_formula
                                    ,data = dat
                                    ,na.action = na.exclude
                                    ,ntree = 1000   # nr of trees to grow
                                    ,keep.forest = TRUE
                                    )
rf_mod

rf_pred <- raster::predict(object = layers_cut, model = rf_mod, type = "prob")
# plotting results in two layers, one with the prediction for each class (0 and 1)
terra::plot(rf_pred, col = mypalette, breaks = colorbreaks)
# we want the prediction of presence only which is named X1 in the figure
rf_pred <- rf_pred[["X1"]]

terra::plot(rf_pred, col = mypalette, breaks = colorbreaks, main = "RF")
terra::plot(presence_centroids, col = "red", pch = ".", add = TRUE)
# Note that Random Forest tend to overfit to observed presence records

# convert probability to prevalence-independent favourability
rf_fav <- fuzzySim::Fav(pred = rf_pred, sample.preval = modEvA::prevalence(model = rf_mod)) 
terra::plot(rf_fav, col = mypalette, breaks = colorbreaks, main = "RF_FAV")
terra::plot(presence_centroids, col = "red", pch = "-", add = TRUE)
# note that rf overfits to strips of presence missing from the main data source

# END OF MODEL 7: RANDOM FOREST
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# MODEL 8:  PRESENCE ABSENCE MODELLING: BOOSTED REGRESSION TREES BRT

# we use the function gbm in pacakge gbm
# we also recycle the formula used for the glm approach

gbm_formula <- glm_formula

# here we use the defualt values
# there are important parameters to explore in depth, like
# shinkage: learning rate; and interaction.depth = tree complexity
gbm_mod <- gbm::gbm(formula = gbm_formula, data = dat)

gbm_pred <- raster::predict(object = layers_cut, model = gbm_mod, type  ="response")
terra::plot(gbm_pred, col = mypalette, breaks = colorbreaks, main = "GBM")
terra::plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# convert probability to prevalence-independent favourability
gbm_fav <- fuzzySim::Fav(pred = gbm_pred
                        ,sample.preval = modEvA::prevalence(model = gbm_mod))
terra::plot(gbm_fav, col = mypalette, breaks = colorbreaks, main = "GBM_FAV")
terra::plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# END OF MODEL 8: BRT
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 3:             VISUALLY COMPARING MODELS
par(mfrow = c(4, 2))
terra::plot(bioclim_pred, col = mypalette, main = "Bioclim")
terra::plot(maxent_pred, col = mypalette, main = "MaxEnt")
terra::plot(maxnet_pred, col = mypalette, main = "Maxnet")
terra::plot(glm_fav, col = mypalette, main = "GLM_Fav")
terra::plot(gam_fav, col = mypalette, main = "GAM_fav")
terra::plot(rf_fav, col = mypalette, main = "RF_Fav")
terra::plot(gbm_fav, col = mypalette, main = "GBM_Fav")

# END OF STEP 3
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 4:    SAVE PREDICTIONS IN FORM OF RASTER MAPS ON DISK

# save prediction raster maps to a folder on disk:
outdir_preds <- paste0("../outputs/predictions_", myspecies)
if (!file.exists(outdir_preds)) 
    dir.create(outdir_preds)

terra::writeRaster(bioclim_pred, filename = paste0(outdir_preds, "/bioclim_pred.tif"))
terra::writeRaster(domain_pred, filename = paste0(outdir_preds, "/domain_pred.tif"))
raster::writeRaster(maxent_pred, filename = paste0(outdir_preds, "/maxent_pred.tif"))
raster::writeRaster(maxnet_pred, filename = paste0(outdir_preds, "/maxnet_pred.tif"))
terra::writeRaster(glm_pred, filename = paste0(outdir_preds, "/glm_pred.tif"))
terra::writeRaster(gam_pred, filename = paste0(outdir_preds, "/gam_pred.tif"))
terra::writeRaster(rf_pred, filename = paste0(outdir_preds, "/rf_pred.tif"))
terra::writeRaster(gbm_pred, filename = paste0(outdir_preds, "/gbm_pred.tif"))
terra::writeRaster(glm_fav, filename = paste0(outdir_preds, "/glm_fav.tif"))
terra::writeRaster(gam_fav, filename = paste0(outdir_preds, "/gam_fav.tif"))
terra::writeRaster(rf_fav, filename = paste0(outdir_preds, "/rf_fav.tif"))
terra::writeRaster(gbm_fav, filename = paste0(outdir_preds, "/gbm_fav.tif"))

# END OF STEP 4
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 5:              SAVE MODEL OBJECTS ON DISK

outdir_models <- paste0("../outputs/models_", myspecies)
if (!file.exists(outdir_models)) 
    dir.create(outdir_models)

# we use the native function saveRDS to write objects to file
saveRDS(bioclim_mod, paste0(outdir_models, "/bioclim_mod.rds"))
saveRDS(domain_mod, paste0(outdir_models, "/domain_mod.rds"))
saveRDS(maxent_mod, paste0(outdir_models, "/maxent_mod.rds"))
saveRDS(maxnet_mod, paste0(outdir_models, "/maxnet_mod.rds"))
saveRDS(glm_mod, paste0(outdir_models, "/glm_mod.rds"))
saveRDS(gam_mod, paste0(outdir_models, "/gam_mod.rds"))
saveRDS(rf_mod, paste0(outdir_models, "/rf_mod.rds"))
saveRDS(gbm_mod, paste0(outdir_models, "/gbm_mod.rds"))

# END OF STEP 5
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 6:           APPLY PREDICTIONS TO THE DATA TABLE

# important! variables must have exactly the same name as in the models!!
# Note that this time we perform predictions based on the dataframe
# instead of raster maps
# It is more convinient to then store the outcome into dataframe for
# further analyses
# Note that below we use default values to build each model!

dat$bioclim_pred <- dismo::predict(object = bioclim_mod, x = dat)
dat$domain_pred  <- dismo::predict(object = domain_mod, x = dat)
dat$maxent_pred  <- dismo::predict(object = maxent_mod, x = dat)
dat$maxnet_pred  <- as.vector(raster::predict(object = maxnet_mod
                                             ,newdata = dat
                                             ,type = "cloglog"
                                             )
                             )
dat$glm_pred     <- as.vector(raster::predict(object = glm_mod
                                             ,newdata = dat
                                             ,type = "response"
                                             )
                             )
dat$gam_pred     <- as.vector(raster::predict(object = gam_mod
                                             ,newdata = dat
                                             ,type = "response"
                                             )
                             )
# for randomForest, consider that predictions include probabilities for each class
# and we want to take predictions corresponding to presence (i.e., "1")
# The result has two columns "0" for absence; "1", presence
dat$rf_pred     <- raster::predict(object = rf_mod, newdata = dat, type = "prob")[,"1"]
dat$gbm_pred     <- raster::predict(object = gbm_mod, newdata = dat, type = "response")

# convert probabilities by presence-absence models to favourability
# we use the predictions calculated from the dataframe, which values are also stored 
# in the dataframe
dat$glm_fav <- fuzzySim::Fav(pred = dat$glm_pred
                            ,sample.preval = modEvA::prevalence(model = glm_mod)
                            )
dat$gam_fav <- fuzzySim::Fav(pred = dat$gam_pred
                            ,sample.preval = modEvA::prevalence(model = gam_mod)
                            )
dat$rf_fav <- fuzzySim::Fav(pred = dat$rf_pred
                            ,sample.preval = modEvA::prevalence(model = rf_mod)
                            )
dat$gbm_fav <- fuzzySim::Fav(pred = dat$gbm_pred
                            ,sample.preval = modEvA::prevalence(model = gbm_mod)
                            )

# save data to disk
write.csv(dat, paste0("../outputs/dat+preds_",myspecies,".csv"),row.names = FALSE)

# compare models predictions visually

# column names with prediction values
pred_cols <- names(dat)[grep(pattern = "_pred|_fav", x = names(dat))]

# map predictions from the datatable

# first we need to convert dat to vector map of points
dat_points <- terra::vect(dat, geom = c("x", "y"))

# map all prediction columns in one window
terra::plot(dat_points, pred_cols, col = mypalette, breaks = colorbreaks, cex = 0.1
            ,axes = FALSE, legend = FALSE, nr = 4, nc = 2, mar = c(0,0,1,0))

# one can also map predictions one by one
terra::plot(dat_points, "bioclim_pred", col = mypalette, breaks = colorbreaks, cex = 0.5
            ,main = "Bioclim")
# and so on...

# END OF STEP 6
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 7:         COMPARE PREDICTIONS OF DIFFERENT MODELS

# we check for correlations in the predicted values by each model
pred_corrs <- cor(dat[,pred_cols])
pred_corrs
min(pred_corrs, na.rm = TRUE)

# visual correlation matrix
# using package corrplot
par(mfrow = c(1,1))
corrplot::corrplot(pred_corrs, method = "ellipse", type = "upper"
                  ,addCoef.col = "wheat3"
                  ,addCoefasPercent = TRUE, tl.col = "black"
                  )

# END OF STEP 7
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# END OF MODELLING METHODS
#----------------------------------------------------------------------
#               LINUX INSTALLING RJAVA

# -- BASH
#Update the repository listings.
#sudo apt update -y

#Install the Java Runtime Environment (JRE) and Java Development Kit (JDK).
#sudo apt install -y openjdk-8-jdk openjdk-8-jre

#Update where R expects to find various Java files.
#sudo R CMD javareconf

#If you get an error about jni.h not being found, then try this:
#sudo R CMD javareconf JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/

# -- R
#Install the package.
#install.packages("rJava")

#If you have a RStudio session open, then exit and restart it. 
#This is important (a running RStudio session will not pick up these changes!).
