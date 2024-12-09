# Get Burlington land cover proportions 
# Because Julia packages are not loading
# And I am too lazy to troubleshoot

library(raster)
library(terra)
library(landscapemetrics)

# load raster
land <- raster('BurlingtonLandCover2016.grd')

# get proportions of each land cover class
amounts <- table(as.vector(land))

# remove water and combine some of the classes -------
amounts <- amounts[-1] # remove water

# Combine low-density urban types
amounts[1] <- sum(amounts[1], amounts[2])
amounts <- amounts[-2]
names(amounts)[1] <- "DevelopedLow"

# Combine high-density urban types
amounts[2] <- sum(amounts[2], amounts[3])
amounts <- amounts[-3]
names(amounts)[2] <- "DevelopedHigh"

# rename barren habitat
names(amounts)[3] <- "Other"

# Combine deciduous forest
amounts[4] <- sum(amounts[4], amounts[6])
amounts <- amounts[-6]
names(amounts)[4] <- "DeciduousForest"

# rename coniferous forest
names(amounts)[5] <- "Conifers"

# combine shrubby habitats
amounts[6] <- sum(amounts[6], amounts[7])
amounts <- amounts[-7]
names(amounts)[6] <- "Shrub"

# rename pasture and crops
names(amounts)[7] <- "Pasture"
names(amounts)[8] <- "Crops"

# Combine wetlands
amounts[9] <- sum(amounts[9], amounts[10])
amounts <- amounts[-10]
names(amounts)[9] <- "Wetlands"

# Get proportions
props <- amounts/sum(amounts)

# Reclassify land cover -------------
land <- reclassify(land, rcl = matrix(data = c(11,22,24,43,71,95,
                                             NA,21,23,41,52,90),
                                    byrow = F))

# Calculate aggregation index ---------
lsm_l_ai(land)
