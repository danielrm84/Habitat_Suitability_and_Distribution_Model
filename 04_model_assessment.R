######################################################################
#
#                           MODEL ASSESSMENT
#                    
#
# Based on 
# Physalia course on SDMs (Oct 2021)
# by Marcia Barbosa
#
# Author: Daniel Romero Mujalli
# University of Greifswald
#
# Important! openCV was not working for me on Linux. There was a
#            problem with my R version (too old!). Therefore, the 
#            code below for block-cross-validation might need some 
#            small adjustments to make it work.
#
#
#######################################################################

#----------------------------------------------------------------------
# STEP 0:                 

# Required packages
pkgs <- c("rgbif", "sdmpredictors", "fuzzySim", "terra"
        , "modEvA", "dismo", "maxnet", "gam", "randomForest", "gbm"
        , "corrplot", "ecospat", "blockCV", "raster"
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

# WARNING! We are using information from previous steps / scripts!

# define a colour palette and colour breaks (intensity), so that all maps of model 
# predictions show a comparable gradient
# HCL Hue Chroma Luminance
mypalette <- hcl.colors( n = 10               # the number of colors
                       , palette = "viridis" # desired palette (check ?hcl.colors)
                       )
colorbreaks <- seq(from = 0, to = 1, by = 0.1)

# END OF STEP 0
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# STEP 1: IMPORT DATA

# select the species
#myspecies <- "Pleurodeles waltl" # the species of interest
myspecies <- "Callitrichidae" # family of marmosets and tamarins

dat <- read.csv(paste0("../outputs/dat+preds_", myspecies, ".csv"))
head(dat)

# load the maps (tif files)
pred_files <- list.files(paste0("../outputs/predictions_", myspecies),full.names = TRUE)

# convert to raster maps
pred_maps <- terra::rast(pred_files)
terra::plot(pred_maps, col = mypalette)

# let us inform better the layer names
# we remove the filepath and keep the base name of the file
pred_names <- basename(tools::file_path_sans_ext(pred_files))
names(pred_maps) <- pred_names
terra::plot(pred_maps, col = mypalette)

# END OF STEP 1
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 2: EVALUATE MODELS ON THE TRAINING DATA
names(dat)

## AREA UNDER THE ROCK CURVE AUC
# we use the AUC function in package modEvA

# set plot margins and 3x2 plots on same window
par(mfrow = c(4,2), mar = c(3,2,2,1))

# with native R function
# Evaluates an R expression in an environment constructed from data
# data: data to use for constructing an environment
# expr: expression to evaluate possibly modifying (a copy of) the original data.
with(data = dat, expr = modEvA::AUC(obs = presence, pred = bioclim_pred, main = "Bioclim"))
with(dat, modEvA::AUC(obs = presence, pred = domain_pred, main = "Domain"))
with(dat, modEvA::AUC(obs = presence, pred = maxent_pred, main = "Maxent"))
with(dat, modEvA::AUC(obs = presence, pred = maxnet_pred, main = "Maxnet"))
with(dat, modEvA::AUC(obs = presence, pred = glm_pred, main = "GLM"))
with(dat, modEvA::AUC(obs = presence, pred = gam_pred, main = "GAM"))
with(dat, modEvA::AUC(obs = presence, pred = rf_pred, main = "RF"))
with(dat, modEvA::AUC(obs = presence, pred = gbm_pred, main = "GBM"))

## THRESHOLD-BASED CLASSIFICATION METRICS
# we use the threshMeasures function in package modEvA

# select the preferred metrics
classif_metrics <- c("CCR", "Sensitivity", "Specificity", "Precision", "Recall"
                    ,"TSS", "kappa"
                    )

# one can choose a threshold that optimizes classification performance for a
# particular metric for each method
# using the function optiThresh from modEvA
# threshold for each model:
my_thresh <- with(dat, modEvA::optiThresh(obs = dat$presence
                                 ,pred = maxnet_pred#glm_pred
                                 ,interval = 0.001
                                 ,measures = "TSS"
                                 ,optimize = "each"
                                 )
                )
my_thresh <- my_thresh$optimals.each$threshold
my_thresh

# alternatively, one can use species prevalence as the threshold for all models
# threshold for all models
my_thresh <- "preval"

par(mfrow = c(4, 2), mar = c(5, 2, 2, 1))

with(dat, modEvA::threshMeasures(obs = presence
                        ,pred = bioclim_pred
                        #,standardize = FALSE (uncomment in case data isnt standardized)
                        ,thresh = my_thresh
                        ,measures = classif_metrics
                        ,ylim = c(0, 1)
                        ,main = "Bioclim"
                        )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = domain_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "Domain"
                                )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = maxent_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "Maxent"
                                )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = maxnet_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "Maxnet"
                                )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = glm_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "GLM"
                                )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = gam_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "GAM"
                                )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = rf_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "RF"
                                )
    )
with(dat, modEvA::threshMeasures(obs = presence
                                ,pred = gbm_pred
                                #,standardize = FALSE
                                ,thresh = my_thresh
                                ,measures = classif_metrics
                                ,ylim = c(0, 1)
                                ,main = "GBM"
                                )
    )

# MILLER CALIBRATION
# calculate and plot Miller calibration line
# we use the MillerCalib function in pacakge modEvA

par(mfrow = c(4, 2), mar = c(3, 2, 2, 1))

with(dat, modEvA::MillerCalib(obs = presence, pred = bioclim_pred, main = "Bioclim"))
with(dat, modEvA::MillerCalib(obs = presence, pred = domain_pred, main = "Domain"))
with(dat, modEvA::MillerCalib(obs = presence, pred = maxent_pred, main = "Maxent"))
with(dat, modEvA::MillerCalib(obs = presence, pred = maxnet_pred, main = "Maxnet"))
with(dat, modEvA::MillerCalib(obs = presence, pred = glm_pred, main = "GLM")) 
# mind that, by definition, for GLM MillerCalib is always perfect with the 
# training data; it is useful for evaluating calibration on external 
# (out-of-sample) data
with(dat, modEvA::MillerCalib(obs = presence, pred = gam_pred, main = "GAM"))
with(dat, modEvA::MillerCalib(obs = presence, pred = rf_pred, main = "RF"))
with(dat, modEvA::MillerCalib(obs = presence, pred = gbm_pred, main = "GBM"))

## COMPUTE BOYCE INDEX

# Continuous Boyce Index (CBI) is an additional assessment tool. 
# The Boyce index requires presence data only and measures by how much model 
# predictions differ from random distribution of observed presence across 
# the prediction gradient.
# The continuous values of the Boyce index vary between −1 and +1. 
# It is thus the most appropriate metric in the case of presence-only
# models.
# Positive values indicate a model where predictions are consistent with 
# the distribution of actual presence data, values close to zero mean that 
# the model is no different from a random model and negative values indicate 
# counter predictions (e.g. predicting no occurrence in areas where actual 
# presence is recorded)

# we use the ecospat.boyce function in ecospat package
# The index is calculated on raster variables + presence points
presence_coords <- dat[dat$presence == 1, c("x", "y")] #x,y coordinate info
# Bioclim
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$bioclim_pred)
                      ,obs = presence_coords)
title("Bioclim")
# Domain
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$domain_pred)
                      ,obs = presence_coords)
title("Domain")
# Maxent
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$maxent_pred)
                      ,obs = presence_coords)
title("Maxent")
# Maxnet
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$maxnet_pred)
                      ,obs = presence_coords)
title("Maxnet")
# GLM
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$glm_pred)
                      ,obs = presence_coords)
title("GLM")
# GAM
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$gam_pred)
                      ,obs = presence_coords)
title("GAM")
# RF
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$rf_pred)
                      ,obs = presence_coords)
title("RF")
# GBM
ecospat::ecospat.boyce(fit = raster::raster(pred_maps$gbm_pred)
                      ,obs = presence_coords)
title("GBM")

# Alternatively, the same index can be computed using the dataframe
# with the observed (presence) and predicted (model) values
# Bioclim
ecospat::ecospat.boyce(obs = dat[dat$presence == 1, "bioclim_pred"]
                      ,fit = dat[ , "bioclim_pred"]
                      )
title("Bioclim")

# and so on.. I will not write the code for the other model predictions

# So far, the model assessment has been performed on the same data on
# which the models were trained
# Another possibilitiy of model assessment is the use of a test-
# independent dataset (cross-validation)
# Block cross-validation is regarded as the most appropriate method
# when using SDMs
#----------------------------------------------------------------------
# END OF STEP 2:
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 3: SPATIAL-BLOCK CROSS-VALIDATION

#### FIRST, CREATE THE BLOCK PARTITIONS

# we use functions spatialBlock and spatialAutoRange from the blockCV 
# package

# import the variable layers over which to calculate spatial autocorrelation:
layers_cut_folder <- paste0("../outputs/sdmpredictors/layers_cut_", myspecies)
layers_cut <- terra::rast(list.files(layers_cut_folder
                                     ,pattern = "\\.tif$"
                                     , full.names = TRUE)
                          )
terra::plot(layers_cut)

# convert 'dat' to a spatial vector map with its coordinate reference system:
# given that the 'dat' coordinates have same CRS as 'layers_cut', we can
# pass the latter as a reference
names(dat)
dat_points <- terra::vect(dat, geom = c("x", "y"), crs = terra::crs(layers_cut))  
par(mfrow = c(1, 1))  # 1 plot per window
terra::plot(dat_points, cex = 0.1)

# calculate the range of spatial autocorrelation in the modelling variables, 
# in the study area using the spatialAutoRange function of blockCV package
# Also, import the vector of selected variables saved with the previous 
# script:
vars_sel <- read.csv(paste0("../outputs/vars_sel_", myspecies , ".csv")
                    ,stringsAsFactors = FALSE)[ , 1]
vars_sel
sarange <- blockCV::spatialAutoRange(subset(raster::stack(layers_cut), vars_sel))
sarange$range

# One can use sarange$range as 'theRange' if it isnt too large for the species
# of interest.
# However, blocks based on the autocorrelation range may be too large for the 
# models to be able to capture the species-environment relationship adequately
# An alternative could be to get spatial blocks of 150 km instead:
# One can set a particular seed of random numbers, so the next command always 
# provides the same set of blocks:
set.seed(42)  # 42: the ultimate secret of the Universe
blocks <- blockCV::spatialBlock(speciesData = as(dat_points, "Spatial")
                               ,rasterLayer = raster::raster(layers_cut[[1]])
                               , theRange = 150000, k = 5, selection = "random"
                               )  
# The function also admits using argument species = "presence"; 
# see ?spatialBlock for more info

blocks$folds
blocks$foldID

# add the fold ID to the data frame:
dat$foldID <- blocks$foldID
head(dat)

dat_points$foldID <- blocks$foldID
par(mfrow = c(1, 1))
terra::plot(dat_points, "foldID")  # each fold has a different colour

#### SECOND, COMPUTE MODELS AND GET PREDICTIONS LEAVING OUT ONE BLOCK
#### FOLD EACH TIME

folds <- sort(unique(dat$foldID))
folds

names(dat)

# Here I will give the example for only one model (i.e., glm)
# Additional models can be considered within the loop, following the
# same logic

for (f in folds)
{
    # inform the progress
    message("modelling outside fold ", f, "...\n")

    # select current training data
    dat_train <- subset(dat, foldID != f) # leave f out

    # remember that some models, e.g., bioclim, use presence-only data
    # uncomment accordingly
    # dat_train_presonly <- subset(dat_train, presence == 1)

    # define formula for the glm model
    glm_formula <- as.formula(paste("presence ~"
                             ,paste(vars_sel, collapse = " + ")
                                   )
                            )
    # fit the model
    glm_mod_fold <- glm(formula = glm_formula
                       ,family = binomial
                       ,data = dat_train
                       )
    # write predictions based on the current testdata fold f
    # into the dataframe
    dat[,paste0("glm_fold", f, "_pred")] <- raster::predict(object = glm_mod_fold
                                                           ,newdata = dat
                                                           ,type = "response"
                                                           ) 

    # call the garbage collector to clean up memory before next loop iteration
    gc()

}# end of loop, over folds

# see the new predictions added to the the data frame:
head(dat)

#### THIRD, EVALUATE EACH MODEL ON ITS VALIDATION FOLD

# Again, I only use the example on the glm model

# select model predictions
fold_cols <- grep(pattern = "_fold", x = names(dat))
names(dat)[fold_cols] # a check

# models in use
models <- c("glm") # include here other models, if necessary
# performance metrics
measures <- c("AUC", "TSS", "MCS")

# create and empty table to receive the cross-validation results
crossval <- as.data.frame(matrix(nrow = length(folds)
                                ,ncol = length(measures) * length(models)
                                )
                         )

# we define the column names of the empty dataframe taking advantage
# of the outer function (see ?outer)
colnames(crossval) <- c(outer(X = models, Y = measures, FUN = paste, sep  = "_"))
# check the dataframe where data will be stored (should be empty: i.e.,NAs)
crossval

# set plot area
par(mfrow = c(5, 2), mar = c(2, 2, 2, 1))
# loop over models and folds, compute the metrics and write results into
# crossval dataframe
for (m in models)  for (f in folds) {
  fold_name <- paste0("fold", f)
  fold_col <- names(dat)[grep(paste0(m, "_fold", f), names(dat))]
  fold_dat <- subset(dat, foldID == f)
  # get AUC (also, ROC curve)
  roc_col <- paste(m, "ROC curve", sep = "_")
  crossval[f, roc_col] <- modEvA::AUC(obs = fold_dat[ , "presence"]
                                     ,pred = fold_dat[ , fold_col]
                                     ,simplif = TRUE
                                     ,plot = TRUE
                                     ,main = paste(m, "AUC")
                                     )
  # get TSS
  tss_col <- paste(m, "TSS", sep = "_")
  crossval[f, tss_col] <- modEvA::threshMeasures(obs = fold_dat[ , "presence"]
                                                ,pred = fold_dat[ , fold_col]
                                                ,thresh = my_thresh
                                                ,measures = "TSS"
                                                ,simplif = TRUE
                                                ,standardize = FALSE
                                                ,main = paste(m, "TSS")
                                                )
  # get MCS
  mcs_col <- paste(m, "MCS", sep = "_")
  crossval[f, mcs_col] <- modEvA::MillerCalib(obs = fold_dat[ , "presence"]
                                             ,pred = fold_dat[ , fold_col]
                                             ,main = paste(m, "Miller line")
                                             )$slope
}

crossval
# order by names
crossval <- crossval[ , order(names(crossval))]
crossval

# boxplots of cross-validation metrics:
# remember Miller calibration slope (MCS) should ideally be close to 1 
# (it does not mean the bigger the value = better)
par(mfrow = c(1, 1), mar = c(7, 3, 2, 1))
boxplot(crossval, col = rep(1:length(models), each = length(measures)), las = 2)
abline(h = 1, col = "darkgrey", lty = 2)  

# we can also plot just the mean cross-validation performance for each model
# (but note that the mean is quite limited information):
crossval_means <- sapply(crossval, mean, na.rm = TRUE)
barplot(crossval_means, ylim = c(0, max(crossval, na.rm = TRUE))
       ,col = rep(1:length(models)
       ,each = length(measures))
       ,las = 2
       )
abline(h = 1, col = "darkgrey", lty = 2)

# save the cross-validation results to a file on disk:
write.csv(crossval, paste0("../outputs/crossval_", myspecies, ".csv")
         ,row.names = FALSE
         )

# END OF STEP 3:
#----------------------------------------------------------------------
# END OF MODEL ASSESSMENT