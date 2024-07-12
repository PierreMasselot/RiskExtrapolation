################################################################################
#
#   Risk Extrapolation project
#   
#   Prepare time series data
#
################################################################################

#---------------------------
# Read time series data
#---------------------------

# Read mortality data
mortdata <- fread("data/mortality.csv.gz")

# Read temperature series and merge
tempdata <- fread("data/tmean.csv.gz")
tsdata <- merge(mortdata, tempdata, all.x = T)

# Order
setkey(tsdata, city_code, date)

# Create date-related variables
tsdata[, ":="(year = year(date), dow = weekdays(date))]

#---------------------------
# Read Metadata
#---------------------------

# Create all city-age combinations
agelabs <- grep("deaths_[[:digit:]]", names(tsdata), value = T) |> 
  gsub(pattern = "deaths_", replacement = "")
metadf <- expand.grid(agegroup = agelabs, 
  city_code = unique(tsdata$city_code))

# Load data from EUcityTRM and merge
metadata_spatial <- read.csv("data/metadata_spatial.csv.gz")
metadf <- merge(metadf, metadata_spatial) 

# Load age-specific demographic data
metadata_age <- read.csv("data/metadata_age.csv.gz")


#---------------------------
# Separate observed and predicted cities
#---------------------------

# Draw observed cities for the first stage
set.seed(1)
obs <- sample(unique(tsdata$city_code), nobs)

# Separate time series data
tspred <- tsdata[!city_code %in% obs]
tsdata <- tsdata[city_code %in% obs]

# Include info in metadata
metadf <- mutate(metadf, obs = city_code %in% obs)

#---------------------------
# Additional data for results
#---------------------------

# Read map of Italy
italymap <- st_read("data/italymap.shp")

#----- Common basis to represent curves

# Estimate an overall empirical distribution of temperature
tmeandist <- tempdata[, .(per = predper, tper = quantile(tmean, predper / 100)), 
  by = city_code]
ovper <- tmeandist[, .(tmean = mean(tper)), by = per]$tmean

# Create average basis
ovknots <- ovper[predper %in% varper]
ovbasis <- onebasis(ovper, fun = varfun, degree = vardegree, knots = ovknots)

# Axis locations for plots
ovaxis <- ovper[predper %in% axisper]

# Create basis and pred for acceptable MMTs 
mmtbasis <- ovbasis[between(predper, mmprange[1], mmprange[2]),]
mmtper <- ovper[between(predper, mmprange[1], mmprange[2])]

#----- Useful objects

# Number of spline coefficients
nc <- ncol(ovbasis)

# Number of cities
n <- nrow(metadata_spatial)

# Number of city-ages
na <- nrow(metadf)

