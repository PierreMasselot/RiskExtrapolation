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

#----- Read

# File to read
ita_file <- "italy_cities_2022-10-07.RData"

# Read
load(paste0("data/", ita_file))

#----- Aggregate age groups

dlist <- lapply(dlist, function(d){
  
  # Get indices and age groups
  mort_vars <- grep("all_[[:digit:]]", names(d), value = T)
  agemax <- suppressWarnings(as.numeric(substr(mort_vars, 7, 8)))
  agemax[length(agemax)] <- Inf
  
  # Cut as groups
  agegrps <- cut(agemax, c(agebreaks, Inf), agelabs, include.lowest = T,
    right = F)
  
  # Aggregate
  dagg <- tapply(as.list(d[mort_vars]), agegrps, function(x) Reduce("+", x))
  dagg <- do.call(data.frame, dagg)
  names(dagg) <- sprintf("all_%s", agelabs)
  
  # put in original df
  d[mort_vars] <- NULL
  cbind(d, dagg)
})

#---------------------------
# Read Metadata
#---------------------------

#----- Load data from EUcityTRM

# File 
metafile <- "metadata.csv"

# Download from Zenodo
if (!file.exists(paste0("data/", metafile))){
  download_zenodo("10.5281/zenodo.7672108", path = "data",
    files = list("metadata.csv"))
}

# Load and select Italian cities
metadata <- read.csv("data/metadata.csv")

# Select italian cities
metadata <- merge(cities, dplyr::select(metadata, !c(CNTR_CODE, URAU_NAME, lon, 
  lat, region, LABEL, cntr_name, nmiss, mcc_code, cityname, country, inmcc)), 
  by.x = "CITY_CODE", by.y = "URAU_CODE")

#----- Prepare data for the second-stage

# Draw observed cities for the first stage
set.seed(1)
obs <- sample.int(nrow(metadata), nobs)
metadata$obs <- seq_len(nrow(metadata)) %in% obs

# Create all city-age combinations
stage2df <- expand.grid(agegroup = agelabs, city = metadata$CITY_CODE)

# Add name, and geographical information
stage2df <- merge(stage2df, 
  metadata[, c("CITY_CODE" ,"CITY_NAME", "geozone","lon", "lat", "obs")], 
  by.x = "city", by.y = "CITY_CODE")

# Separate data.frame
stage2_obs <- subset(stage2df, obs)
stage2_pred <- subset(stage2df, !obs)

#---------------------------
# Additional data for results
#---------------------------

# Extract map of Italy
italymap <- gisco_get_countries(year = "2020", epsg = "4326", country = "IT")

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

#----- Dimensions

# Number of spline coefficients
nc <- length(varper) + 2

# Number of cities
n <- nrow(metadata)

# Number of city-ages
na <- nrow(stage2df)

# Length of time series
nts <- nrow(dlist[[1]])
