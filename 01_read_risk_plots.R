# https://sma.org.au/sma-site-content/uploads/2021/02/SMA-Extreme-Heat-Policy-2021-Final.pdf
# https://www.sciencedirect.com/science/article/pii/S1440244018300215?via%3Dihub
# https://journals.physiology.org/doi/full/10.1152/japplphysiol.00191.2018
# https://www.omnicalculator.com/physics/wet-bulb

library(imager)

function_interpolation <- function(data) {
  
  data$x1 <- data$row ^ 1
  data$x2 <- data$row ^ 2
  data$x3 <- data$row ^ 3
  data$x4 <- data$row ^ 4
  
  model <- lm(formula = col~x1+x2+x3+x4, 
              data = data)
  
  return(coefficients(model))
  
}

nrow.fig <- 2; ncol.fig <- 3
pdf("outdata/plot/risk_functions_physical_activity.pdf", 
    width = 4 * ncol.fig, height = 3 * nrow.fig)

layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))

list_coefficients <- list()

# This loop goes through the 5 risk classifications and the three levels of risk
# It reads the images from the original document, extract the lines defining the
# levels of risks.
# From the extracted points we fit a function of degree 4 to approximate the
# thresholds.
# We plot and save the results to have a visual look if everything worked
# correctly

for(sport_classification in 1:5) {
  
  data <- lapply(1:3, function(irisk){
    
    image <- load.image(
      paste0(getwd(), "/indata/other/risk_images/risk", sport_classification, 
             "_edited_threshold", irisk, ".png"))
    
    # plot(image) # Image. Width: 648 pix Height: 506 pix Depth: 1 Colour channels: 3
    mat <- image[,,1,1]
    mat <- mat != 1
    # image(mat)
    
    data <- data.frame(which(mat == TRUE, arr.ind = TRUE))
    
    # change to the coordinates in the original plots (26, 44) and (100, 0)
    x <- seq(26, 44, length.out = nrow(mat))
    y <- seq(100, 0, length.out = ncol(mat))
    data$row <- x[data$row]
    data$col <- y[data$col]
    
    return(data)
    
  })
  
  list_coefficients[[sport_classification]] <- sapply(data, function_interpolation)
  
  plot(1, 1, type = "n", xlim = c(20, 50), ylim = c(0, 100),
       xlab = "Air Temperature (ºC)", ylab = "Relative humidity (%)",
       main = paste0("Sport Risk Classification ", sport_classification))
  
  x <- 20:50
  
  col <- c("#92d050",
           "#ed7d31",
           "#ff0000")
  
  for(i in 1:3) {
    points(data[[i]])
    lines(x, 
          cbind(x^0, x^1, x^2, x^3, x^4) %*% 
            list_coefficients[[sport_classification]][, i, drop = FALSE], 
          col = col[i], lwd = 2)
  }
  
  abline(v = c(26, 44))
  abline(h = c(0, 100))
  
  
}

dev.off()

save(list_coefficients, file = "outdata/list_coefficients_thresholds_risk_classification.RData")