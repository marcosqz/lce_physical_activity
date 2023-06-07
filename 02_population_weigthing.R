# WARNING: Data Unavailability in Repository
# Please note that the following code assumes the existence of a data file in the repository,
# but the data file is currently not available. Make sure to acquire the necessary data
# and update the file path accordingly before running this code.

# This code take the population data from eurostat in a format of grid cell of
# 1x1 km but in shapefile format and calculated the population in each cell of
# the climate date and store it in the same raster format than the climate data

library(sf)
library(raster)
library(ncdf4)
library(data.table)
library(dplyr)

path_publicdata <- "X:/"
year <- 2018 # Select the year of the eurostat population data (2006, 2011, 2018)

# Random climate raster as a model to extract main characteristic of the rasters
raster_temp <- brick(
  paste0(path_publicdata, 
         "adaptation/climate/reanalysis/era5_land/rh/2001", 
         "/era5-land_global_hourly_rh_20010101.nc"), varname = "rh")
raster_temp <- raster_temp[[15]] # raster at 14:00

if (year == 2018) {
  shp_popu <- st_read("indata/other/geostat/2018/JRC_POPULATION_2018.shp")
} else if (year == 2011) {
  shp_popu <- st_read(
    "indata/other/geostat/2011/GEOSTATReferenceGrid/Grid_ETRS89_LAEA_1K-ref_GEOSTAT_POP_2011_V2_0_1.shp")
} else if (year == 2006) {
  shp_popu <- st_read("indata/other/geostat/2006/Grid_ETRS89_LAEA_1K_ref_GEOSTAT_2006.shp")
}

shp_popu <- st_transform(shp_popu, 4326)

### Now we will create a raster object (identical to the one of the climate)
### with the population values in each cell

# European extension
extension <- extent(-32, 45, 27, 72)

# Crop the model raster (snap = "out" to ensure we take the left extreme of the
# longitudes in the extension domain)
raster_template <- crop(raster_temp, extension, snap = "out")
rm(raster_temp)

# Change climate variable for a vector from 1 to total_length that we will use
# as ids
values(raster_template) <- 1:length(values(raster_template))

# Get centroids of the grids
points_popu <- st_centroid(shp_popu) # data.frame with te id of the population cell and the value of the population
rm(shp_popu)

# Read (if needed the csv with the population data)

if (year == 2018){
  
  names(points_popu)[names(points_popu) == "TOT_P_2018"] <- "POPULATION"

} else if(year == 2011){
  
  pop_grid <- fread("indata/geostat/2011/GEOSTAT_grid_POP_1K_2011_V2_0_1.csv")
  pop_grid <- subset(pop_grid, select = c("TOT_P", "GRD_ID", "CNTR_CODE"))
  pop_grid <- aggregate(pop_grid$TOT_P, by = list(GRD_ID = pop_grid$GRD_ID), FUN = sum)
  
  names(pop_grid)[names(pop_grid) == "x"] <- "POPULATION"
  
  points_popu <- left_join(points_popu, pop_grid)

} else if (year == 2006) {
  
  names(points_popu)[names(points_popu) == "GRD_INSPIR"] <- "GRD_ID"
  
  pop_grid <- fread("indata/geostat/2006/GEOSTAT_grid_EU_POP_2006_1K_V1_1_1.csv")
  pop_grid <- subset(pop_grid, select = c("POP_TOT", "GRD_ID"))
  pop_grid <- aggregate(pop_grid$POP_TOT, by = list(GRD_ID = pop_grid$GRD_ID), FUN = sum)
  
  names(pop_grid)[names(pop_grid) == "x"] <- "POPULATION"
  
  points_popu <- left_join(points_popu, pop_grid)
  
}

# Tidy population data in the points data
points_popu <- na.omit(points_popu) # exclude grids with NA
points_popu <- points_popu[points_popu$POPULATION != 0,] # exclude pts with 0 population

# For each population point extract to which raster cell of the climate raster corresponds 
# (raster cell identified by the id 1:total_length generated before)
extract_id <- data.frame(raster::extract(raster_template, points_popu)); names(extract_id) <- "id_cell"
# assign to the object with the population points the id of the grid cell to which corresponds
points_popu$id_cell <- extract_id$id_cell
rm(extract_id)

# Sum the population in each point by the id of the cell in the climate data just assigned
data_pop <- sapply(split(points_popu$POPULATION, points_popu$id_cell), sum)
data_pop <- data.frame(
  id = names(data_pop),
  pop = data_pop
)
# > head(data_pop)
# id pop
# 6800 6800   5
# 6801 6801  21
# 6814 6814 112
# 6818 6818  89
# 6819 6819  19
# 6842 6842  19

# Create a new raster with the same features of the climate raster that will store
# the population in each grid cell

raster_pop <- raster_template; rm(raster_template)
values_pop <- rep(NA, length(raster_pop))
values_pop[as.numeric(data_pop$id)] <- data_pop$pop
values(raster_pop) <- values_pop
names(raster_pop) <- "population"

rm(values_pop)

# plot(!is.na(raster_pop))

# Save the resulting population raster

if (year == 2018) {
  writeRaster(raster_pop, filename = paste0(getwd(), "/outdata/file/raster_population_2018.nc"), overwrite = TRUE)
} else if(year == 2011) {
  writeRaster(raster_pop, filename = paste0(getwd(), "/outdata/file/raster_population_2011.nc"), overwrite = TRUE)
} else if(year == 2006) {
  writeRaster(raster_pop, filename = paste0(getwd(), "/outdata/file/raster_population_2006.nc"), overwrite = TRUE)
}