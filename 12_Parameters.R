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
  'Population density' = "popdens", 
  'Social isolation' = "isol", GPD = "gdp", 'Unemployment rate' = "unempl", 
  'Educational level' = "educ", 'Deprivation rate' = "depriv", 
  'Hospital bed rate' = "bedrates", Imperviousness = "imperv", 
  'Tree cover density' = "tree", 'Grassland' = "grass", 
  'Wetness and water' = "water", 'Small woody features' = "woody",
  Elevation = "elevation", 'Distance to coast' = "coast_dist", NDVI = "ndvi", 
  'PM2.5' = "pm25", NO2 = "no2", 'Average temperature' = "tmean", 
  'Temperature range' = "trange")

# Number of composite vulnerability indices
npc <- 6

# Number of simulations for eCI
nsim <- 1000

#----- Results

# Temperature percentile grid for ERFs
predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))

# Acceptable MMP range
mmprange <- c(25, 99)

# Percentiles to display on x-axis
axisper <- c(1, 25, 50, 75, 99)

# Denominator for death rates
byrate <- 10^5

# Ages to which to display results
agebreaks <- c(45, 65, 75, 85)
