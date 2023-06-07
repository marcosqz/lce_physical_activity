library(ggplot2)
library(ggpubr)
library(dplyr)
library(ncdf4)

### LOAD DATA ####

rm(list = ls())

category_pa <- 1
type_plot <- "jpg" # "pdf", "jpg"
geo_class <- "european_regions" # "european_regions", "countries", "nuts"

# Read and process the shapefile with the classification of regions
# into climatic zones

nuts_polygons <- giscoR::gisco_get_nuts(resolution = "01", year = "2021",
                                        epsg = "3035", nuts_level = "2")
nuts_polygons <- nuts_polygons[!(grepl("FRY", nuts_polygons$NUTS_ID) |
                                   nuts_polygons$CNTR_CODE == "TR" | 
                                   nuts_polygons$NUTS_ID == "NO0B"),]
nreg <- nrow(nuts_polygons)
nuts_id <- nuts_polygons$NUTS_ID
countries_id <- substr(nuts_polygons$NUTS_ID, 1, 2)

if(geo_class == "european_regions") {
  data_countries <- readxl::read_excel(
    "indata/other/[LCDE 2024] Country names and groupings.xlsx", skip = 1)
  
  zones <- data_countries[match(countries_id, data_countries$`Eurostat code`),]$`European sub-region (UN geoscheme)`
  name_zones <- unique(zones)[order(unique(zones))]
  name_zones <- name_zones[name_zones %in% c("Eastern Europe",  
                                             "Northern Europe", 
                                             "Southern Europe",
                                             "Western Europe")]
} else if(geo_class == "countries") {
  zones <- countries_id
  name_zones <- unique(zones)[order(unique(zones))]
} else if(geo_class == "nuts") {
  zones <- nuts_id
  name_zones <- unique(zones)[order(unique(zones))]
}

rm(nuts_polygons)

# Read the raster with the population data
netcdf_pop <- nc_open("outdata/file/raster_population_2018.nc")
population <- ncvar_get(netcdf_pop, varid = "population")
population <- population[, ncol(population):1]

rm(netcdf_pop)

matrix_population <- matrix(rep(c(population), nreg), ncol = nreg)

load("outdata/file/regions_coverage_matrix.RData")
regional_population <- 
  colSums(matrix_population * regions_coverage_matrix, na.rm = TRUE)
rm(matrix_population, population, regions_coverage_matrix, nreg, nuts_id)

year_period <- list(
  1990:2000,
  2001:2011,
  2012:2022,
  2022)

# for(category_pa in 1:5) {
  
  ### PLOT1: INTRA-DAY. CIRCULAR. REGIONAL. YEARLY - WARMING STRIPES ####
  
  years <- year_period[[1]][1]:tail(year_period[[3]], 1)
  num_colors <- length(years)  # Number of colors in the palette
  original_palette <- c("#FFF3F3",
                        "#FEDFD1",
                        "#F5CAAE",
                        "#D7B38F",
                        "#B7A88A",
                        "#9CA592",
                        "#82A19C",
                        "#6598A4",
                        "#4C8AA7",
                        "#3B78A5",
                        "#2F649E",
                        "#285193",
                        "#233C86",
                        "#1F2676",
                        "#1A0C65")
  
  palette <- colorRampPalette(original_palette)(num_colors)
  
  # Loop for each type of region-classification
  plot <- lapply(1:length(name_zones), function(j) {
    
    izone <- name_zones[j]
    
    pop_zone <- sum(regional_population[zones %in% izone])
    
    data <- lapply(years, function(iyear) {
      load(paste0(
        "outdata/file/parallel_risk_sport_category",
        category_pa, "_", iyear,".RData"))
      risk_daily <- sapply(risk, function(x) colSums(x$regions_exposure[zones %in% izone,,drop = FALSE]))
      y <- rowSums(risk_daily)
      data.frame(x = 0:24, y = c(y, y[1]) / pop_zone, year = iyear)
    })
    
    data <- do.call(rbind, data)
    
    # y_breaks <- round(seq(0, ceiling(max(data$y)), length = 5), 1)
    # Round number to the closest multiply of 5
    if(max(data$y) > 0) {
      x1 <- max(data$y) / 4
      if(x1 < 10) {
        x2 <- x1 / (10 ^ floor(log10(x1)))
      } else {
        x2 <- x1 / 10
      }
      proportion <- x1/x2
      distance_breaks <- ceiling(x2*2) / 2 * proportion
    } else {
      distance_breaks <- 0
    }
    
    y_breaks <- seq(0, 4 * distance_breaks, by = distance_breaks)
    
    text_data <- data.frame(x = 0,
                            y = y_breaks,
                            label = y_breaks,
                            year = "all")
    
    
    data_mean <- data %>%
      group_by(x) %>%
      summarise(y = mean(y))
    data_mean$year <- "all"
    
    p <- ggplot(data, aes(x = x, y = y, group = factor(year))) +
      geom_line(aes(col = factor(year)), linewidth = 1.5) +
      geom_line(data = data_mean, linewidth = 1, col = "red") +
      geom_point(data = data_mean, aes(x = x, y = y), size = 3, col = "red") +
      theme_minimal() + 
      theme(panel.grid.major = element_line(linewidth = 1.00, linetype = 1)) +
      scale_x_continuous("", breaks = 0:23) +
      scale_y_continuous("",
                         breaks = y_breaks,
                         limits = c(0, max(y_breaks))) +
      coord_polar() +
      scale_color_manual(values = palette) + 
      theme(
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0.1, -0.5, -0.5, -0.5), "cm")) +
      labs(color = "") +
      geom_text(data = text_data, aes(x = x, y = y, label = label)) +
      ggtitle(izone)
    
    return(p)
    
  })
  
  if(type_plot == "jpg") {
    jpeg(paste0("outdata/plot/", geo_class, "/pa", category_pa, "_intraday_regional_yearly.jpg"),
         width = 6.5, height = 6
         , units = "in", res = 500, quality = 100
    )
  } else if(type_plot == "pdf") {
    pdf(paste0("outdata/plot/", geo_class, "/pa", category_pa, "_intraday_regional_yearly.pdf"),
        width = 6.5, height = 6)
  }
  
  ggarrange(
    plotlist = plot, ncol = 2, nrow = 2,
    common.legend = TRUE, legend = "right"
  )
  dev.off()
  
  rm(plot, num_colors, original_palette, palette, years)
  
  ### PLOT2: INTRA-DAY. CIRCULAR. REGIONAL. 10-YEAR PERIODS ####
  
  plot <- lapply(1:length(name_zones), function(j) {
    
    izone <- name_zones[j]
    
    pop_zone <- sum(regional_population[zones %in% izone])
    
    data <- lapply(1:3, function(iperiod) {
      
      years <- year_period[[iperiod]]
      nyears <- length(years)
      factor_name <- paste0(years[1], "-", years[length(years)])
      
      y <- rowSums(
        sapply(years, function(iyear) {
          load(paste0(
            "outdata/file/parallel_risk_sport_category", 
            category_pa, "_", iyear,".RData"))
          y <- rowSums(sapply(risk, function(x) 
            colSums(x$regions_exposure[zones %in% izone,,drop = FALSE]))) / (pop_zone * nyears)
          y
        })
      )
      
      data.frame(x = 0:24, y = c(y, y[1]), year = factor_name)
      
    })
    
    data <- do.call(rbind, data)
    
    # Round number to the closest multiply of 5
    if(max(data$y) > 0) {
      x1 <- max(data$y) / 4
      if(x1 < 10) {
        x2 <- x1 / (10 ^ floor(log10(x1)))
      } else {
        x2 <- x1 / 10
      }
      proportion <- x1/x2
      distance_breaks <- ceiling(x2*2) / 2 * proportion
    } else {
      distance_breaks <- 0
    }
    
    y_breaks <- seq(0, 4 * distance_breaks, by = distance_breaks)
    text_data <- data.frame(x = 0,
                            y = y_breaks,
                            label = y_breaks,
                            year = "all")
    
    p <- ggplot(data, aes(x = x, y = y, group = factor(year))) +
      geom_line(aes(col = factor(year)), linewidth = 1.5) +
      geom_point(aes(x = x, y = y, col = factor(year)), size = 3) +
      theme_minimal() + 
      theme(panel.grid.major = element_line(linewidth = 1.00, linetype = 1)) +
      scale_x_continuous("", breaks = 0:23) +
      scale_y_continuous("",
                         breaks = y_breaks,
                         limits = c(0, max(y_breaks))) +
      coord_polar() +
      scale_color_manual(values = c("black", "blue", "red")) + 
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_blank(),
            plot.margin = unit(c(0.1, -0.5, -0.5, -0.5), "cm")) +
      labs(color = "") +
      geom_text(data = text_data, aes(x = x, y = y, label = label)) +
      ggtitle(izone)
    
    return(p)
    
  })
  
  if(type_plot == "jpg") {
    jpeg(paste0("outdata/plot/", geo_class, "/pa", category_pa, "_intraday_regional_10years.jpg"),
         width = 6, height = 6.2
         , units = "in", res = 500, quality = 100
    )
  } else if(type_plot == "pdf") {
    pdf(paste0("outdata/plot/", geo_class, "/pa", category_pa, "_intraday_regional_10years.pdf"),
        width = 6, height = 6.2)
  }
  
  ggarrange(
    plotlist = plot, ncol = 2, nrow = 2,
    common.legend = TRUE, legend = "bottom"
  )
  dev.off()
  
  rm(plot)
  
  ### PLOT3: INTER-SEASONAL. CIRCULAR. REGIONAL. 10-YEAR PERIODS ####
  
  plot <- lapply(1:length(name_zones), function(j) {
    
    izone <- name_zones[j]
    
    pop_zone <- sum(regional_population[zones %in% izone])
    
    data <- lapply(1:3, function(iperiod) {
      
      years <- year_period[[iperiod]]
      nyears <- length(years)
      factor_name <- paste0(years[1], "-", years[length(years)])
      
      y <- rowSums(
        sapply(years, function(iyear) {
          load(paste0(
            "outdata/file/parallel_risk_sport_category", 
            category_pa, "_", iyear,".RData"))
          dates <- lubridate::month(seq(as.Date(paste0(iyear, "-01-01")), as.Date(paste0(iyear, "-12-31")), by = "days"))
          y <- colSums(sapply(risk, function(x) colSums(x$regions_exposure[zones %in% izone,,drop = FALSE])))  / (pop_zone * nyears)
          sapply(1:12, function(imonth) sum(y[dates == imonth]))
        })
      )
      
      data.frame(x = c(1:12, 13), y = c(y, y[1]), year = factor_name)
    })
    
    data <- do.call(rbind, data)
    
    # Round number to the closest multiply of 5
    if(max(data$y) > 0) {
      x1 <- max(data$y) / 4
      if(x1 < 10) {
        x2 <- x1 / (10 ^ floor(log10(x1)))
      } else {
        x2 <- x1 / 10
      }
      proportion <- x1/x2
      distance_breaks <- ceiling(x2*2) / 2 * proportion
    } else {
      distance_breaks <- 0
    }
    
    y_breaks <- seq(0, 4 * distance_breaks, by = distance_breaks)
    text_data <- data.frame(x = 1,
                            y = y_breaks,
                            label = y_breaks,
                            year = "all")
    
    p <- ggplot(data, aes(x = x, y = y, group = factor(year))) +
      geom_line(aes(col = factor(year)), linewidth = 1.5) +
      geom_point(aes(x = x, y = y, col = factor(year)), size = 3) +
      theme_minimal() + 
      theme(panel.grid.major = element_line(linewidth = 1.00, linetype = 1)) +
      scale_x_continuous("", breaks = 1:12, labels = month.abb[1:12]) +
      scale_y_continuous("",
                         breaks = y_breaks,
                         limits = c(0, max(y_breaks))) +
      coord_polar() +
      scale_color_manual(values = c("black", "blue", "red")) + 
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_blank(),
            plot.margin = unit(c(0.1, -0.5, -0.5, -0.5), "cm")) +
      labs(color = "") +
      geom_text(data = text_data, aes(x = x, y = y, label = label)) +
      ggtitle(izone)
    
    return(p)
    
  })
  
  if(type_plot == "jpg") {
    jpeg(paste0("outdata/plot/", geo_class, "/pa", category_pa, "_interseasonal_regional_10years.jpg"),
         width = 6, height = 6.2
         , units = "in", res = 500, quality = 100
    )
  } else if(type_plot == "pdf") {
    pdf(paste0("outdata/plot/", geo_class, "/pa", category_pa, "_interseasonal_regional_10years.pdf"),
        width = 6, height = 6.2
    )
  }
  
  ggarrange(
    plotlist = plot, ncol = 2, nrow = 2,
    common.legend = TRUE, legend = "bottom"
  )
  dev.off()
  
  rm(plot)
  
  ### PLOT4: BOTH INTRA-DAY AND INTER-SEASONAL. CIRCULAR. EUROPE. 10-YEAR PERIODS ####
  
  if(geo_class == "european_regions") {
    
    # Population for Europe as the sum of the population of all regions
    pop_zone <- sum(regional_population)
    
    # Dataset with the intra-day PA risky hours per person
    data_day <- lapply(1:3, function(iperiod) {
      
      years <- year_period[[iperiod]]
      nyears <- length(years)
      factor_name <- paste0(years[1], "-", years[length(years)])
      
      y <- rowSums(
        sapply(years, function(iyear) {
          load(paste0("outdata/file/parallel_risk_sport_category", 
                      category_pa, "_", iyear,".RData"))
          y <- rowSums(sapply(risk, function(x) colSums(x$regions_exposure)))  / (pop_zone * nyears)
          y
        })
      )
      
      data.frame(x = 0:24, y = c(y, y[1]), year = factor_name)
      
    })
    
    data_day <- do.call(rbind, data_day)
    
    # Dataset with the inter-seasonal PA risky hours per person.
    # We calculate the risky hours per person in each month
    data_season <- lapply(1:3, function(iperiod) {
      
      years <- year_period[[iperiod]]
      nyears <- length(years)
      factor_name <- paste0(years[1], "-", years[length(years)])
      
      y <- rowSums(
        sapply(years, function(iyear) {
          load(paste0(
            "outdata/file/parallel_risk_sport_category", category_pa, 
            "_", iyear,".RData"))
          dates <- lubridate::month(seq(as.Date(paste0(iyear, "-01-01")), as.Date(paste0(iyear, "-12-31")), by = "days"))
          y <- colSums(sapply(risk, function(x) colSums(x$regions_exposure)))  / (pop_zone * nyears)
          sapply(1:12, function(imonth) sum(y[dates == imonth]))
        })
      )
      
      data.frame(x = c(1:12, 13), y = c(y, y[1]), year = factor_name)
    })
    
    data_season <- do.call(rbind, data_season)
    
    # Parameters for the plots
    
    Generate_Breaks <- function(max_y) {
      x1 <- max_y / 4
      if(x1 < 10) {
        x2 <- x1 / (10 ^ floor(log10(x1)))
      } else {
        x2 <- x1 / 10
      }
      proportion <- x1/x2
      distance_breaks <- ceiling(x2*2) / 2 * proportion
      return(distance_breaks)
    }
    
    distance_breaks_day <- Generate_Breaks(max(data_day$y))
    distance_breaks_season <- Generate_Breaks(max(data_season$y))
    
    y_breaks_day <- seq(0, 4 * distance_breaks_day, by = distance_breaks_day)
    y_breaks_season <- seq(0, 4 * distance_breaks_season, by = distance_breaks_season)
    
    text_data_day <- data.frame(x = 0,
                                y = y_breaks_day,
                                label = y_breaks_day,
                                year = "all")
    text_data_season <- data.frame(x = 1,
                                   y = y_breaks_season,
                                   label = y_breaks_season,
                                   year = "all")
    
    # Plot for day and season
    
    p_day <- ggplot(data_day, aes(x = x, y = y, group = factor(year))) +
      geom_line(aes(col = factor(year)), size = 1.5) +
      geom_point(aes(x = x, y = y, col = factor(year)), size = 3) +
      theme_minimal() + 
      theme(panel.grid.major = element_line(size = 1.00, linetype = 1)) +
      scale_x_continuous("", breaks = 0:23) +
      scale_y_continuous("",
                         breaks = y_breaks_day,
                         limits = c(0, max(y_breaks_day))) +
      coord_polar() +
      scale_color_manual(values = c("black", "blue", "red")) + 
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_blank(),
            plot.margin = unit(c(0.1, -0.5, -0.5, -0.5), "cm")) +
      labs(color = "") +
      geom_text(data = text_data_day, aes(x = x, y = y, label = label)) +
      ggtitle("Intra-day")
    
    p_season <- ggplot(data_season, aes(x = x, y = y, group = factor(year))) +
      geom_line(aes(col = factor(year)), size = 1.5) +
      geom_point(aes(x = x, y = y, col = factor(year)), size = 3) +
      theme_minimal() + 
      theme(panel.grid.major = element_line(size = 1.00, linetype = 1)) +
      scale_x_continuous("", breaks = 1:12, labels = month.abb[1:12]) +
      scale_y_continuous("",
                         breaks = y_breaks_season,
                         limits = c(0, max(y_breaks_season))) +
      coord_polar() +
      scale_color_manual(values = c("black", "blue", "red")) + 
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_blank(),
            plot.margin = unit(c(0.1, -0.5, -0.5, -0.5), "cm")) +
      labs(color = "") +
      geom_text(data = text_data_season, aes(x = x, y = y, label = label)) +
      ggtitle("Inter-sesonal")
    
    # Save the plot
    
    if(type_plot == "jpg") {
      jpeg(paste0("outdata/plot/pa", category_pa, "_both_intraday_interseaonal_europe_10years.jpg"),
           width = 6*1.25, height = 3.1*1.25
           , units = "in", res = 500, quality = 100
      )
    } else if(type_plot == "pdf") {
      pdf(paste0("outdata/plot/pa", category_pa, "_both_intraday_interseaonal_europe_10years.pdf"),
          width = 6*1.25, height = 3.1*1.25)
    }
    ggarrange(p_day, p_season, common.legend = TRUE, legend = "bottom")
    dev.off()
    
    rm(pop_zone, 
       data_day, data_season,
       distance_breaks_day, distance_breaks_season,
       y_breaks_day, y_breaks_season, 
       text_data_day, text_data_season,
       p_day, p_season)
    
  }
  
  ### TABLE1. NUMBERS of relative change and absolute number of late-early hours ####
  
  plot <- lapply(1:length(name_zones), function(j) {
    
    izone <- name_zones[j]
    
    pop_zone <- sum(regional_population[zones %in% izone])
    
    data <- lapply(1:4, function(iperiod) {
      
      years <- year_period[[iperiod]]
      nyears <- length(years)
      factor_name <- paste0(years[1], "-", years[length(years)])
      
      y <- rowSums(
        sapply(years, function(iyear) {
          load(paste0(
            "outdata/file/parallel_risk_sport_category", 
            category_pa, "_", iyear,".RData"))
          y <- rowSums(sapply(risk, function(x) colSums(x$regions_exposure[zones %in% izone,,drop = FALSE]))) / (pop_zone * nyears)
          y
        })
      )
      
      data.frame(time = 0:23, risky_hour_person = y, year = factor_name)
      
    })
    
    data <- do.call(rbind, data)
    data$zone <- izone
    
    data1 <- data[data$year == paste0(year_period[[1]][1], "-", tail(year_period[[1]], 1)),]
    max_hour <- data1[order(data1$risky_hour_person, decreasing = TRUE),]$time[1:4]
    min_hour <- (0:23)[!(0:23 %in% max_hour)]
    
    rm(data1)
    
    data_aggregated1 <- data %>%
      filter(time %in% max_hour) %>%
      group_by(year) %>%
      summarize(hottest = sum(risky_hour_person))
    
    data_aggregated2 <- data %>%
      filter(time %in% min_hour) %>%
      group_by(year) %>%
      summarize(not_hottest = sum(risky_hour_person))
    
    data_aggregated <- data.frame(
      year = data_aggregated1$year,
      hottest = data_aggregated1$hottest,
      not_hottest = data_aggregated2$not_hottest,
      zone = izone
    )
    
    return(list(
      max_hour = c(izone, paste0(year_period[[1]][1], "-", tail(year_period[[1]], 1)), max_hour),
      data = data,
      data_aggregated = data_aggregated
    ))
    
  })
  
  table_max_hour <- do.call(rbind, lapply(plot, function(x) x$max_hour))
  table_data <- do.call(rbind, lapply(plot, function(x) x$data))
  table_data_hottest <- do.call(rbind, lapply(plot, function(x) x$data_aggregated))
  
  write.csv(table_max_hour, row.names = FALSE,
            file = paste0("outdata/plot/", geo_class, "/pa", category_pa, "_table_intraday_maxhours_regional_10years.csv"))
  write.csv(table_data, row.names = FALSE,
            file = paste0("outdata/plot/", geo_class, "/pa", category_pa, "_table_intraday_data_regional_10years.csv"))
  write.csv(table_data_hottest, row.names = FALSE,
            file = paste0("outdata/plot/", geo_class, "/pa", category_pa, "_table_intraday_data_hottest_regional_10years.csv"))
  
  # table interseasonal
  
  plot <- lapply(1:length(name_zones), function(j) {
    
    izone <- name_zones[j]
    
    pop_zone <- sum(regional_population[zones %in% izone])
    
    data <- lapply(1:3, function(iperiod) {
      
      years <- year_period[[iperiod]]
      nyears <- length(years)
      factor_name <- paste0(years[1], "-", years[length(years)])
      
      y <- rowSums(
        sapply(years, function(iyear) {
          load(paste0(
            "outdata/file/parallel_risk_sport_category", 
            category_pa, "_", iyear,".RData"))
          dates <- lubridate::month(seq(as.Date(paste0(iyear, "-01-01")), as.Date(paste0(iyear, "-12-31")), by = "days"))
          y <- colSums(sapply(risk, function(x) colSums(x$regions_exposure[zones %in% izone,,drop = FALSE]))) / (pop_zone * nyears)
          sapply(1:12, function(imonth) sum(y[dates == imonth]))
        })
      )
      
      data.frame(month = c(1:12), risky_hour_person = y, year = factor_name)
    })
    
    data <- do.call(rbind, data)
    data$zone <- izone
    
    return(data)
    
  })
  
  table_data <- do.call(rbind, plot)
  
  write.csv(table_data, row.names = FALSE,
            file = paste0("outdata/plot/", geo_class, "/pa", category_pa, "_table_interseasonal_data_regional_10years.csv"))
  
  rm(list = ls())
  
# }