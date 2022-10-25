################################################################################
#
#   Risk Extrapolation project
#   
#   Parameters
#
################################################################################

#------------------------
# Parameters
#------------------------

#----- Data management

# Paths
path_nuts <- "V:/VolumeQ/AGteam/Eurostat/Regional by NUTS classification (reg)"
path_urau <- "V:/VolumeQ/AGteam/Eurostat/Urban Audit (urb_cgc)"

# Years to be considered in the analysis
year <- 2011:2021

# Age groups
agebreaks <- c(0, 45, 65, 75, 85)
agelabs <- c(paste(agebreaks[-length(agebreaks)], agebreaks[-1], sep = "-"), 
  sprintf("%i+", agebreaks[length(agebreaks)]))

#----- First-stage analysis

# Number of cities in first-stage
nobs <- 60

# Exposure dimension
varfun <- "bs"
varper <- c(10,75,90)
vardegree <- 2

# Lag dimension
maxlag <- 21
lagfun <- "ns"
lagknots <- logknots(maxlag, 3)

#----- Second-stage model

# List of metapredictors in composite indices
metaprednames <- c(Population = "pop", 'Population above 65' = "prop_65p", 
  'Population density' = "popdens", 'Life expectancy' = "lifexp_00", 
  'Social isolation' = "isol", GPD = "gdp", 'Unemployment rate' = "unempl", 
  'Educational level' = "educ", 'Deprivation rate' = "depriv", 
  'Hospital bed rate' = "bedrates", Imperviousness = "imperv", 
  'Tree cover density' = "tree", 'Grassland' = "grass", 
  'Wetness and water' = "water", 'Small woody features' = "woody",
  Elevation = "elevation", 'Distance to coast' = "coast_dist", NDVI = "ndvi", 
  'PM2.5' = "pm25", NO2 = "no2", 'Average temperature' = "tmean", 
  'Temperature range' = "trange")

# Number of simulations for eCI
nsim <- 20

#----- Results

# Temperature percentile grid for ERFs
predper <- c(seq(0.1,1,0.1), 2:98, seq(99,99.9,0.1))

# Acceptable MMP range
mmprange <- c(25, 99)

# Percentiles to display on x-axis
axisper <- c(1, 25, 50, 75, 99)

# Denominator for death rates
byrate <- 10^5
