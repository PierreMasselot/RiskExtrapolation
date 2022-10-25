################################################################################
#
#   New extensions
#   0. Prepare second stage
#
################################################################################

source("01_Packages.R")
load("data/FirstStage.RData")

#-------------------------
# Prepare second-stage data
#-------------------------

# Create all city-age combinations
stage2df <- expand.grid(agegroup = agelabs, city = names(stage1res))

# Add name, and geographical information
stage2df <- merge(stage2df, 
  metadata[, c("CITY_CODE" ,"CITY_NAME", "geozone","lon", "lat", "obs")], 
  by.x = "city", by.y = "CITY_CODE")

#----- Extract first stage results

# Get all results in big list
allresults <- unlist(stage1res, recursive = F)

# Get coefs and vcov
coefs <- t(sapply(allresults, "[[", "coef"))
vcovs <- lapply(allresults, "[[", "vcov")

#----- Meta-predictors

# Extract meta-predictors
metapred <- data.matrix(metadata[,metaprednames])

# Expand to several age groups
metapred <- metapred[match(stage2df$city, metadata$CITY_CODE),]

#-------------------------
# Common parameters to all illustrated models
#-------------------------

#----- Formulas

# Baseline fixed effects formula (might add some regional background)
fixform <- coefs ~ 1

# Random effects
ranform <- ~ 1|city

#----- Other useful parameters

# Number of coefficients
nc <- ncol(coefs)

# Number of cities
n <- length(dlist)

# Number of age groups
nage <- length(agebreaks - 1)

# Length of time series
nts <- nrow(dlist[[1]])

#-------------------------
# Construct basis representing an average temperature distribution
#-------------------------

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
