daily_phisical_activity_indicator <- function(
  iday, 
  year,
  ilon,
  ilat,
  risk_coef,
  path_publicdata, 
  filenames_dates,
  population,
  out_matrix,
  nreg,
  regions_coverage_matrix,
  position_missing,
  position_available) {
  
  # Open temperature and relative humidity netcdfs for the corresponding day
  
  netcdf_temp <- nc_open(
    paste0(path_publicdata,
           "adaptation/climate/reanalysis/era5_land/rh/", year,
           "/era5-land_global_hourly_t2m-d2m_", filenames_dates[iday], ".nc"))
  
  netcdf_rh <- nc_open(
    paste0(path_publicdata,
           "adaptation/climate/reanalysis/era5_land/rh/", year,
           "/era5-land_global_hourly_rh_", filenames_dates[iday], ".nc"))
  
  # Format temperature and relative humidity
  
  temp <- ncvar_get(netcdf_temp, varid = "t2m",
                    start = c(1, 1 ,1), count = c(-1, -1, -1))
  temp <- temp[ilon, (ncol(temp):1)[ilat], ] - 273.15
  temp[temp < 20] <- 20 # Temperatures under 20?C will be always risk 0
  temp[temp > 47] <- 47 # Temperatures above 50?C Will be always risk 3

  rh <- ncvar_get(netcdf_rh, varid = "rh",
                  start = c(1, 1, 1), count = c(-1, -1, -1))
  rh <- rh[ilon, (ncol(rh):1)[ilat], ] * 100
  
  # Put NAs in the cells we are not interested in
  
  for(i in 1:dim(temp)[3]) {
    temp[,,i][out_matrix] <- NA
    rh[,,i][out_matrix] <- NA
  }
  rm(i)
  
  nc_close(netcdf_temp)
  nc_close(netcdf_rh)
  
  # Risk categories are defined by quadratic functions which coefficients
  # are stored in the object "risk_coef"
  temp2 <- temp*temp; temp3 <- temp2*temp; temp4 <- temp3*temp
  
  threshold <- lapply(risk_coef, function(x)
    x[1] + x[2] * temp + x[3] * temp2 + x[4] * temp3 + x[5] * temp4)
  rm(temp2, temp3, temp4)
  
  # Create a netcdg with the risk at each grid cell for the 24 hours in that day
  # Values are 0 (low risk), 1 (moderate risk), 2 (high risk), 3 (extreme risk)
  netcdf_risk <-
    (rh >= threshold[[1]]) +
    (rh >= threshold[[2]]) +
    (rh >= threshold[[3]])
  
  rm(threshold, temp, rh)
  
  # Use closest risk for those cells with population but not risk result
  for(i in 1:dim(netcdf_risk)[3]) {
    netcdf_risk[,,i][position_missing] <-
      netcdf_risk[,,i][position_available]
  }
  rm(i)
  
  # Multiply netcdf of risk with netcdf of population to get the risk weighted
  # by population for all Europe in each level of risk
  
  # total, no considering regions
  risk_weighted <- sapply(0:3, function(risk_level) {
    apply(netcdf_risk, 3, function(x) # the third dimension of netcdf_risk is 24, 1 per hour of the day
      sum((x == risk_level) * population, na.rm = TRUE))
  })
  
  # Risk x person in each region (risk greater than 0) in each region
  regional_exposure <- sapply(1:24, function(ihour) {
    
    vector_risk <- c(netcdf_risk[,,ihour] > 0) # > 0, ie some risk
    
    # Repeat risk vector to the same dimensions as the matrix regions
    # coverage
    matrix_risk <- matrix(rep(vector_risk, nreg),
                          ncol = nreg)
    matrix_population <- matrix(rep(c(population), nreg),
                                ncol = nreg)

    if(any(vector_risk == TRUE, na.rm = TRUE)) {
      
      map <- colSums(matrix_population * matrix_risk * regions_coverage_matrix, 
                     na.rm = TRUE)
      
    } else {
      
      map <- rep(0, nreg)
      
    }
    
    return(map)
    
  })
  
  # Return the person-risk in Europe and all the regions
  
  return(
    list("europe_risk" = risk_weighted,
         "regions_exposure" = regional_exposure))
  
}