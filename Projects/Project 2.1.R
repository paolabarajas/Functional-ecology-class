# 1. Load packages -------------------------------------------------------------

# install.packages('terra')
# install.packages('sf')
# install.packages("ggplot2")

library(dplyr)
library(RColorBrewer)
library(ggExtra)
library(terra)
library(sf)
library(ggplot2)

# 2. Set working directory and load data ---------------------------------------

setwd("G:/My Drive/_Teaching/2026_Class_functional ecology/MSc/Projects/P2 data")

traits <- read.csv("trait_plots_Tenerife.csv")
plots <- read.csv("plots_Tenerife.csv")
plots_shp <- read_sf("plot_VT.gpkg")

# Retrieve climatic data from
# https://www.chelsa-climate.org/datasets

raster_path <-  # raster path 
  'G:/My Drive/_Teaching/2026_Class_functional ecology/MSc/Projects/P2 data/CHELSACanaryClim_tasmin_01_1979-2013_V.1.0.tif'
raster_data <- rast(raster_path)  # load raster data


# 3. Data checks ---------------------------------------

plot(plots_shp['plot'])
plot(raster_data)

ggplot(trait, aes(Leaf_area)) + 
  geom_density() + 
  # scale_x_log10() + 
  theme_minimal()

ggplot(trait, aes(Seed_mass)) + 
  geom_density() + 
  # scale_x_log10() + 
  theme_minimal()
