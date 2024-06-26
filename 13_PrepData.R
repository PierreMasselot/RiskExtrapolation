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
fulldata <- merge(mortdata, tempdata, by.x = c("CITY_CODE", "date"), 
  by.y = c("URAU_CODE", "date"))

# Order
setkey(fulldata, CITY_CODE, date)

# Create date-related variables
fulldata[, ":="(year = year(date), dow = weekdays(date))]

# Split
dlist <- split(fulldata[,!"CITY_CODE"], fulldata$CITY_CODE)
dlist <- lapply(dlist, as.data.frame)

#---------------------------
# Read Metadata
#---------------------------

# Load data from EUcityTRM
metadata <- read.csv("data/metadata.csv.gz")

# Draw observed cities for the first stage
set.seed(1)
obs <- sample.int(nrow(metadata), nobs)
metadata$obs <- seq_len(nrow(metadata)) %in% obs

# Create all city-age combinations
agelabs <- grep("all_[[:digit:]]", names(fulldata), value = T) |> 
  gsub(pattern = "all_", replacement = "")
stage2df <- expand.grid(agegroup = agelabs, city = metadata$CITY_CODE)

# Add name, and geographical information
stage2df <- merge(stage2df, 
  metadata[, c("CITY_CODE" ,"CITY_NAME", "geozone", "lon", "lat", "obs")], 
  by.x = "city", by.y = "CITY_CODE")

# Separate data.frame
stage2_obs <- subset(stage2df, obs)
stage2_pred <- subset(stage2df, !obs)

#---------------------------
# Additional data for results
#---------------------------

# Read map of Italy
italymap <- st_read("data/italymap.shp")

#----- Common basis to represent curves

# Estimate an overall empirical distribution of temperature
tmeandist <- t(sapply(dlist, 
  function(x) quantile(x$tmean, predper / 100, na.rm = T)))
ovper <- colMeans(tmeandist)

# Create average basis
ovknots <- ovper[sprintf("%i.0%%",varper)]
ovbasis <- onebasis(ovper, fun = varfun, degree = vardegree, knots = ovknots)

# Axis locations for plots
ovaxis <- ovper[predper %in% axisper]

# Create basis and pred for acceptable MMTs 
mmtbasis <- ovbasis[between(predper, mmprange[1], mmprange[2]),]
mmtper <- ovper[between(predper, mmprange[1], mmprange[2])]

#----- Useful objects

# Number of spline coefficients
nc <- length(varper) + 2

# Number of cities
n <- nrow(metadata)

# Number of city-ages
na <- nrow(stage2df)

# Length of time series
nts <- nrow(dlist[[1]])
