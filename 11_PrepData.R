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
ita_path <- "V:/VolumeQ/AGteam/ISTAT/cities/processed"
ita_file <- "italy_cities_2022-10-07.RData"

# Read
load(paste0(ita_path, "/", ita_file))

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

load("data/EUcityTRM_metadata.Rdata")

# Select italian cities
metadata <- merge(cities, dplyr::select(metadata, !c(CNTR_CODE, URAU_NAME, lon, 
  lat, region, LABEL, cntr_name, nmiss, mcc_code, cityname, country, inmcc)), 
  by.x = "CITY_CODE", by.y = "URAU_CODE")

# Draw observed cities for the first stage
set.seed(1)
obs <- sample.int(nrow(metadata), nobs)
metadata$obs <- seq_len(nrow(metadata)) %in% obs
