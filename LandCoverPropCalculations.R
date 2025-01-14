# Get Burlington land cover proportions 
# Because Julia packages are not loading
# And I am too lazy to troubleshoot

library(raster)
library(terra)
library(landscapemetrics)
library(tidyverse)
library(patchwork)

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
land <- reclassify(land, rcl = matrix(data = c(11,21,22,23,24,41,
                                               43,52,71,90,95,
                                             NA,21,21,23,23,41,41,
                                             52,52,90,90),
                                    byrow = F, ncol = 2))

# Calculate aggregation index ---------
lsm_l_ai(land)

# Plot initial raster ------------
land_spdf <- as(land, "SpatialPixelsDataFrame")
land_df <- as.data.frame(land_spdf)
colnames(land_df) <- c("value", "x", "y")

land_df <- land_df %>%
  mutate(value = factor(value))

levels(land_df$value) <- c("DevLow", "DevHigh", "Other", 
                           "Deciduous", "Conifers", "Shrub",
                           "Pasture", "Crops", "Wetlands")

burly_smallcell <- ggplot(data = land_df, 
                          aes(x = x, y = y, fill = value))+
  geom_tile()+
  scale_fill_viridis_d(name = "Cover Type", end = 0.9)+
  theme(axis.text = element_blank(), axis.title = element_blank())

# Increase cell size ----------------
land_agg <- aggregate(land, fact=(500/30), fun = modal)

agg_spdf <- as(land_agg, "SpatialPixelsDataFrame")
agg_df <- as.data.frame(agg_spdf)
colnames(agg_df) <- c("value", "x", "y")

agg_df <- agg_df %>%
  mutate(value = factor(value))

levels(agg_df$value) <- c("DevLow", "DevHigh", "Other", 
                           "Deciduous", "Conifers", "Shrub",
                           "Pasture", "Crops", "Wetlands")

# Plot rasters -------------------

burly_bigcell <- ggplot(data = agg_df, 
                          aes(x = x, y = y, fill = value))+
  geom_tile()+
  scale_fill_viridis_d(name = "Cover Type", end = 0.9)+
  theme(axis.text = element_blank(), axis.title = element_blank())

(burly_smallcell + burly_bigcell)+
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = 'a')

ggsave(filename = "Burly_and_pixels.jpeg", width = 10, height = 5,
       units = "in")

