######################################################################
#
#          DOWNLOAD AND CLEAN DATA FROM BIODIVERSITY REPOSITORY
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
# deprecated
#pkgs <- c("rgbif", "scrubr", "sdmpredictors", "fuzzySim", "terra")
# new list
pkgs <- c("CoordinateCleaner", "rgbif", "sdmpredictors", "fuzzySim", "terra")

# check if pkgs are available
my_packages <- rownames(installed.packages())

for (pkg in pkgs)
{
    if(!pkg %in% my_packages)
        print(paste("missing_package: ",pkg))
}

# proceed to install any missing package
# Personally, I like to use the mirror
# https://www.stats.bris.ac.uk/R/ 
# (visit https://cran.r-project.org/mirrors.html for more alternatives)

# LINUX USERS:
# CoordinateCleaner relies on rgeos package, which requieres
# sudo apt-get install libgeos-dev
#
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
# STEP 1:            DOWNLOAD SPECIES OCCURRENCE DATA

# Select the species of interest
# Course default species = the Iberian ribbed newt: Pleurodeles waltl
# https://www.gbif.org/species/9620
# For this example, I decided to work on a family of new world
# monkeys
#myspecies <- "Pleurodeles waltl"
myspecies <- "Callitrichidae" # family of marmosets and tamarins

# get the data from GBIF repository using the function
# occ_data from the package rgbif
gbif_data <- rgbif::occ_data(scientificName = myspecies
                            ,hasCoordinate  = TRUE
                            ,limit = 50000 # nr of records to return
                            )

# create spatvector from existing file of countries political borders
# using the vect function of terra package
countries <- terra::vect("../data/countries/world_countries.shp")

# plot the countries on a global map
terra::plot(countries)

# save to png file
# Remember to check the R documentation to learn more on how to create
# image/pdf/csv ... files
png("../outputs/map_countries.png")
    terra::plot(countries)
dev.off()

# Based on the plot axes values, select the xmin, xmax, ymin, ymax
# coordinates of the desired extent
# We use the ext function of the terra package
# For the default species we aim at the Iberian Peninsula:
#my_window <- terra::ext(c(-20, 10, 20, 50))
# In my case, my_window tries to capture part of Latin America
my_window <- terra::ext(c(-100, -30, -40, 20))


# check if the polygon of the extent is on the region of interest
# note that we use the function as.polygon from terra package
terra::plot(countries)
terra::plot(terra::as.polygons(my_window) # extent as polygon
          , border = "red"
          , lwd = 2
          , add = TRUE             # add to current plot
          )

# plot country borders for the region of interest (i.e., 
# within "my_window" only)
terra::plot(countries
           ,xlim = my_window[1:2] # first two elements (xmin, xmax)
           ,ylim = my_window[3:4] # third and fourth element (ymin,ymax)
           )

# Alternatively, one can download biodiversity data within the window
# (extent) of interest
# The records search can be constrained to a Region of interest by
# providing the corresponding range of coordinates (min, max), for both
# latitude and longitude, based on the 
# World Geodetic System 1984 (WGS 84) 
# We take advantage of the collapse argument of the paste0 native 
# function such that the longitude range is reported in "xmin, xmax" 
# format, as required by occ_data function of the rgbif package
my_longitude <- paste0(my_window[1:2]
                      , collapse = ", " # character value to separate 
                                        # result
                      )
my_latitude <- paste0(my_window[3:4]
                     ,collapse = ", "
                     )

# Request biodiversity data within the region of interest:
gbif_data <- rgbif::occ_data(scientificName = myspecies
                            ,hasCoordinate = TRUE
                            ,limit = 50000
                            ,decimalLongitude = my_longitude
                            ,decimalLatitude = my_latitude
                            )
# take a look at the data
gbif_data
# Note: if "records found" is larger than "records returned",
# increase the "limit" argument (or decrease the window size of the
# extent)
# examining data structure
names(gbif_data)
names(gbif_data$meta)
names(gbif_data$data)
# Remember to properly cite the biodiversity data!
rgbif::gbif_citation(gbif_data)

# store occurrence records into a variable
occurrences <- gbif_data$data

# map the occurrence records to see how they look:
# plot the territorial divisions of countries
terra::plot(countries
           ,xlim = range(occurrences$decimalLongitude)
           ,ylim = range(occurrences$decimalLatitude)           
           )
# add the occurrence records to the plot
points(occurrences[, c("decimalLongitude", "decimalLatitude")]
      , pch = 20
      , col = "red"
      )
# compare the distribution with the range map reported for the
# species, for example, at https://www.iucnredlist.org 
# Make sure that the distribution makes sense

# Export data as *.csv file
# (remove lists variable types from the data before exporting)
# check the class of variables in the data:
unique(sapply(occurrences, class))
toremove <- which(sapply(occurrences, class) == "list")
names(occurrences)[toremove]
# remove those that are list
# you may consider making a backup of the original data (variable)
dim(occurrences)
occurrences <- occurrences[,-toremove]
dim(occurrences)
# create a directory and write data as *.csv
occurrence_dir <- "../outputs/species_occurrences"
if(!file.exists(occurrence_dir))
    dir.create(occurrence_dir)
write.csv(occurrences
         ,paste0(occurrence_dir, "/occurrences_", myspecies, "_raw.csv")
         ,row.names = FALSE
         )
# end of writing *.csv file

# From now on, there is no need to download the data again from GBIF
occurrences <- read.csv(paste0(occurrence_dir
                              ,"/occurrences_"
                              ,myspecies
                              ,"_raw.csv")
                        )
head(occurrences)

# END OF STEP 1
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# STEP 2:            CLEANING SPECIES OCCURRENCE DATA

# Data may contain many errors. Therefore, it is important to conduct
# careful mapping, inspection and cleaning

# first, we remove records of absence or zero abundance (if any)
names(occurrences)
# we look at the column "occurrenceStatus"

# check for different indications of "absent" which could be in
# different languages
sort(unique(occurrences$occurrenceStatus))
dim(occurrences)

absence_tags <- c("absent", "Absent", "ABSENT"
                 ,"ausente", "Ausente", "AUSENTE"
                 )
toremove <- which(occurrences$occurrenceStatus %in% absence_tags)
if(length(toremove) > 0)
    occurrences <- occurrences[-toremove, ]
dim(occurrences)

# Further cleaning: problematic coordinate values
# using the scrubr package (deprecated)
# can be all written in one line
#occurrences <- scrubr::coord_unlikely(occurrences)
#occurrences <- scrubr::coord_impossible(occurrences)
#occurrences <- scrubr::coord_imprecise(occurrences)
#occurrences <- scrubr::coord_incomplete(occurrences)

# using clean_coordinates method from CoordinateCleaner package
# check reliability of coordinate values according to the
# following test
# more info under ?CoordinateCleaner::clean_coordinates()
# Alternatively, you can perform manual cleaning of gbif data
# following the nice guide by John Waller 2021-02-17 GBIF
# https://data-blog.gbif.org/post/gbif-filtering-guide/

tests_to_perform <- c("capitals", "centroids", "equal", "gbif"
                     ,"institutions", "outliers","seas", "zeros")
flags <- CoordinateCleaner::clean_coordinates(x = occurrences
                                    ,lon = "decimalLongitude"
                                    ,lat = "decimalLatitude"
                                    ,tests = tests_to_perform
                                    )
# summary of data flagged as problematic records
summary(flags)

# It is recommended to explore in more detail the flagged data
# before deciding to remove those records
# Here we will skip that procedure and proceed to remove the 
# flagged data
dim(occurrences)
occurrences <- occurrences[flags$.summary,]
dim(occurrences)

# add cleaned data to the map (blue color)
# Excluded points remain in red
points(occurrences[, c("decimalLongitude", "decimalLatitude")]
      , pch = 20, col = "blue"
      )

# Now, remove presence data with a reported coordinate uncertainty
# (location error, spatial resolution) larger than 10x10 km2: Here we
# use distance (in meters) from centroid to corner (half the diagonal)
# of a square with 10km side (10000 meters)
max_uncertainty <- (10000 * sqrt(2)) / 2

# deprecated
#occurrences <- scrubr::coord_uncertain(x = occurrences
#                      , coorduncertainityLimit = max_uncertainty
#                                      )

# we remove records with uncertainty greater than the threshold given 
# by max_uncertainty
# We also keep missing values. One can consider removing them as well
# though they might be good records...
dim(occurrences)
flags <- is.na(occurrences$coordinateUncertaintyInMeters) | 
         occurrences$coordinateUncertaintyInMeters < max_uncertainty

# In case, missing uncertainty values are to be removed:
#flags <- !is.na(occurrences$coordinateUncertaintyInMeters) & 
#         occurrences$coordinateUncertaintyInMeters < max_uncertainty

occurrences <- occurrences[flags, ] 
dim(occurrences)        

# filter out records with coordinate precision > 0.01
# we keep missing data!
max_precision <- 0.01
flags <- is.na(occurrences$coordinatePrecision) | 
         occurrences$coordinatePrecision < max_precision
occurrences <- occurrences[flags,]
dim(occurrences)

# Note that this will keep missing data on coordinate precision and 
# uncertainty. Careful mapping and visual inspection are necessary

# add the less uncertain occurrence records to the map
points(occurrences[, c("decimalLongitude", "decimalLatitude")]
      ,pch = 20, col = "turquoise"
      )

# write the cleaned data to *.csv file
write.csv(occurrences
         ,paste0(occurrence_dir, "/occurrences_"
                ,myspecies, "_cleaned.csv"
                )
         ,row.names = FALSE
         )

# list the data on disk so far
list.files(occurrence_dir)

# from now on, there is no need to download and clean the data again

occurrences <- read.csv(paste0(occurrence_dir, "/occurrences_"
                              , myspecies, "_cleaned.csv")
                       )
head(occurrences)
dim(occurrences)

# END OF STEP 2
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# END OF DOWNLOAD AND CLEAN BIODIVERSITY DATA
