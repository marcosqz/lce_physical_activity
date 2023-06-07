### LOAD LIBRARIES

library(raster)
library(rgdal)
library(ncdf4)
library(ggplot2)
library(giscoR)
library(sf)
library(exactextractr)
library(data.table)
library(future.apply)
plan(multicore)

### SET WORKING DIRECTORY AND PATHS
# Define paths depending if the code is executed locally or in the cluster
# Variables: path_publicdata, path_adaptation

execution <- "cluster" # "local" / "cluster"
category_risk <- 5

if (execution == "local") {
  
  path_publicdata <- "X:/"
  path_adaptation <- "Z:/"
  
} else if (execution == "cluster") {
  
  path_publicdata <- "/PROJECTES/PUBLICDATA/"
  path_adaptation <- "/PROJECTES/ADAPTATION/"
  
  setwd(paste0(path_adaptation, 
               "proj/mquijal/p20230502_MQ_lce_physical_activity_from_temp"))
  
}

rm(execution)

### LOAD FUNCTIONS

source("codes/functions/function_daily_phisical_activity_indicator.R")

### CROP DATA TO THE EUROPEAN DOMAIN
# Identify cells (latitude and longitude) which are of interest for this study
# We assume that all temperature and relative humidity files have exactly the
# same structure; i.e. same resolution, number of cells, latitudes, longitudes,
# missing values, ...
# Output: ilat, ilon

# Define the European domain: 27N-72N x 32W-45E
european_extent <- extent(-32, 45, 27, 72)

# Open any raster file 
# (we are interested in latitudes and longitudes, no values)
netcdf_object <- nc_open(
  paste0(path_publicdata,
         "adaptation/climate/reanalysis/era5_land/rh/2018/",
         "era5-land_global_hourly_t2m-d2m_20180101.nc"))

lon <- ncvar_get(netcdf_object, varid = "longitude")
lat <- ncvar_get(netcdf_object, varid = "latitude")
lat <- rev(lat) # put latitudes from negative to positive values

ilon <- lon >= european_extent[1] & lon <= european_extent[2]
ilat <- lat >= european_extent[3] & lat <= european_extent[4]

rm(lat, lon, european_extent, netcdf_object)

### DEFINE CELLS WITHIN THE EUROPEAN BOUNDARIES
# Identify cells that are inside the European boundaries (e.g its shapefile) and
# the correspondence between cells and regions
# Output: out_matrix, regions_coverage_matrix, nreg

# Read the latest European shapefile at NUTS2 level at the highest resolution
nuts_polygons <- gisco_get_nuts(resolution = "01", year = "2021",
                                epsg = "3035", nuts_level = "2")
# Remove some regions in which we are no interested:
# Overseas France - FRY (out of the European extent)
# Turkey - TR (no population data)
# Svalbard - NO0B (out of the European extent)
nuts_polygons <- nuts_polygons[!(grepl("FRY", nuts_polygons$NUTS_ID) |
                                   nuts_polygons$CNTR_CODE == "TR" | 
                                   nuts_polygons$NUTS_ID == "NO0B"),]

# Transform polygons into the sf format (simple features)
# Unite regions to get the European boundaries
# Transform shapefiles to the resolution of the rasters
region_sf <- st_geometry(nuts_polygons)
domain_sf <- st_union(region_sf)

region_sf <- st_transform(region_sf, "+proj=longlat +datum=WGS84")
domain_sf <- st_transform(domain_sf, "+proj=longlat +datum=WGS84")

nreg <- length(region_sf)

rm(nuts_polygons)

# Open a raster object with the proper extension (our population data is 
# already in the proper cropped extension from the previous code) 
raster_object <- raster("outdata/file/raster_population_2018.nc")

# Identify all cells that are partially covered by the European boundaries
raster_coverage <- coverage_fraction(raster_object, domain_sf)
raster_coverage <- raster_coverage[[1]] > 0

# Define a function to transform raster values into vectors in a specific
# order
Raster_To_Vector <- function(raster) {
  
  # Transform raster into vectors with the same order that the vector of
  # the netcdf defined in this study (reverted latitudes)
  # (i.e order(Raster_To_Vector(raster)) == order(c(netcdf)))
  
  matrix <- as.matrix(raster)
  matrix <- matrix[nrow(matrix):1, ]
  vector <- c(t(matrix))
  
  return(vector)
  
}

# Boolean matrix with TRUE for region out the European boundaries
domain_vector <- Raster_To_Vector(raster_coverage)
out_matrix <- matrix(!domain_vector,
                        nrow = sum(ilon == TRUE),
                        byrow = FALSE)
# image(out_matrix) to check that it the creation of the out_matrix works as
# it was supposed

rm(raster_coverage, domain_vector)

# Generate a matrix where columns are regions, rows are cells and values
# define the coverage of each region in each cell
regions_coverage_raster <- lapply(1:length(region_sf), function(ireg) {
  coverage_fraction(raster_object, st_geometry(region_sf[[ireg]]))[[1]]
})

regions_coverage_matrix <- sapply(regions_coverage_raster, Raster_To_Vector)

# We assign each grid cell to the region that covers the most.
# Identify cells that are not completely filled by a regions (so, shared
# by region and sea or region 1 and region2) to avoid lose people in the sea.
# Limitation is the bias caused in the shared borders. Benefit if we consider
# all population
regions_coverage_matrix <- 
  t(sapply(1:nrow(regions_coverage_matrix), function(i) {
    x <- regions_coverage_matrix[i,]
    # the cell has some part of its area within a region x!=0 (if it is 1 all the cell already corresponds fully to a region)
    if(any(x!=0 & x!=1)){
      new_vector <- rep(0, nreg)
      # then we assign all the cell to the region to the region with more coverage
      new_vector[which.max(x)] <- 1 
      return(new_vector)
    } else {
      return(x)
    }
    
  }))

save(regions_coverage_matrix, file = "outdata/file/regions_coverage_matrix.RData")

# regions_coverage_matrix
# 0        1 
# 95403573    74029 
rm(regions_coverage_raster, raster_object, Raster_To_Vector,
   region_sf, domain_sf)

# LOAD POPULATION DATA (FIXED POPULATION DATA 2018)
# Output: population

netcdf_pop <- nc_open("outdata/file/raster_population_2018.nc")
population <- ncvar_get(netcdf_pop, varid = "population")
population <- population[, ncol(population):1]
# We are interested only in values inside the European domain
population[out_matrix] <- NA

rm(netcdf_pop)

# COMPLETE MISSING CLIMATE DATA IN GRID CELLS WITH POPULATION VALUES
# Find the closest point with climate data for the grid cells with population
# data and no climate data
# Here we keep assuming that all the climate raster have exactly the same
# structure and cells with missing values
# Output: position_missing, position_available

# Open any temperature raster file 
netcdf_temp <- nc_open(
  paste0(path_publicdata,
         "adaptation/climate/reanalysis/era5_land/rh/2018/",
         "era5-land_global_hourly_t2m-d2m_20180518.nc"))

temp <- ncvar_get(netcdf_temp, varid = "t2m")
temp <- temp[ilon, (ncol(temp):1)[ilat], 1]
longitude <- ncvar_get(netcdf_temp, varid = "longitude")[ilon]
latitude <- rev(ncvar_get(netcdf_temp, varid = "latitude")[ilat])

nc_close(netcdf_temp)

# We are interested only in values inside the European domain
temp[out_matrix] <- NA

# Identify cells with population data and no climate data
netcdf_missing <- is.na(temp) & !is.na(population)

# Create dataframes with the missing and the available grid cell
df_missing <- data.frame(
  row = longitude[which(netcdf_missing == TRUE, arr.ind = TRUE)[,1]],
  col = latitude[which(netcdf_missing == TRUE, arr.ind = TRUE)[,2]]
)
df_missing$id <- 1:nrow(df_missing)
df_missing$position <- which(netcdf_missing == TRUE, arr.ind = FALSE)

df_available <- data.frame(
  row = longitude[which(!is.na(temp), arr.ind = TRUE)[,1]],
  col = latitude[which(!is.na(temp), arr.ind = TRUE)[,2]]
)
df_available$position <- which(!is.na(temp), arr.ind = FALSE)

rm(temp, longitude, latitude)

# Use function which calculates Euclidian distances between the two datasets
# Source: https://stackoverflow.com/questions/55752064/finding-closest-coordinates-between-two-large-data-sets
mydist <- function(a, b, df1, x, y){

  dt <- data.table(sqrt((df1[[x]]-a)^2 + (df1[[y]]-b)^2))

  return(data.table(Closest.V1  = which.min(dt$V1),
                    Distance    = dt[which.min(dt$V1)]))
}

# Run the function and save closest grids in a dataframe
df_closest <- setDT(df_missing)[, j = mydist(row, col, setDT(df_available),
                                             "row", "col"),
                                by = list(id, row, col)]

# Save vectors with the position of the missing and assigned-available grids
position_missing <- df_missing$position
position_available <- df_available[df_closest$Closest.V1,]$position

rm(df_closest, df_missing, df_available, netcdf_temp, netcdf_missing, mydist)

### LOOP EACH SPORT RISK CATEGORY

# for(category_risk in 4:1) {
  
  ### COEFFICIENT DEFINING THE RISK LEVEL
  # Take the coefficient of the 3 quadratic functions defining the 4 levels of
  # risk in each sport category
  # Output: risk_coef
  load("outdata/list_coefficients_thresholds_risk_classification.RData")
  
  risk_coef <- list_coefficients[[category_risk]]
  risk_coef <- split(t(risk_coef), seq(ncol(risk_coef))) # List format
  rm(list_coefficients)
  
  ### LOOP EACH YEAR IN THE SPORT RISK CATEGORY
  
  # for(year in 1990:2022){
  for(year in 2022:2022){
    
    print(paste0("Category risk: ", category_risk, " - ", 
                 "Year: ", year))
    
    # Define the dates in the corresponding year
    dates <- seq(as.Date(paste0(year, "-01-01")),
                 as.Date(paste0(year, "-12-31")), by = "days")
    filenames_dates <- gsub('-', '', as.character(dates))
    
    gc()
    
    options(future.globals.maxSize = 5000 * 1024 ^ 2)
    # Parallelize the function which calculates the person-risk for all Europe 
    # and each region in the correspoding day
    risk <- future_lapply(seq(length(filenames_dates)), function(iday) {
      
      gc()
      
      daily_phisical_activity_indicator(
        iday = iday, year = year,
        ilon = ilon, ilat = ilat,
        risk_coef = risk_coef,
        path_publicdata = path_publicdata, 
        filenames_dates = filenames_dates,
        population = population,
        out_matrix = out_matrix,
        nreg = nreg,
        regions_coverage_matrix = regions_coverage_matrix,
        position_missing = position_missing,
        position_available = position_available)
      
      
    }, future.seed = TRUE)
    
    save(risk, file = paste0("outdata/file/parallel_risk_sport_category",
                             category_risk,"_", year, ".RData"))
    
  }
# }