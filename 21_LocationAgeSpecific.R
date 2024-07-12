################################################################################
#
#  RiskExtrapolation
#
#  Location and age-specific first-stage model
#
################################################################################

#-------------------------
# Run first-stage
#-------------------------

#----- Prepare parallelisation
ncores <- detectCores()
cl <- makeCluster(max(1, ncores - 2))
registerDoParallel(cl)

#----- Perform first stage on each city

# Loop over cities, splitting the data table
stage1res <- foreach(dat = isplit(tsdata, tsdata$city_code), .combine = rbind, 
  .packages = c("dlnm", "splines", "foreach", "mixmeta", "dplyr")) %dopar% 
{
  
  # Extract from iterator
  city <- dat$key[[1]]
  dat <- as.data.frame(dat$value)
  
  # Define crossbasis
  argvar <- list(fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100, na.rm = T))
  cb <- crossbasis(dat$tmean, lag = maxlag, argvar = argvar,
    arglag = list(fun = lagfun, knots = lagknots))
  
  # Loop on deaths by age group
  agevars <- grep("deaths", names(dat), value = T)
  agegroups <- gsub("deaths_", "", agevars)
  out <- foreach(y = dat[, agevars], agegr = agegroups, .combine = rbind) %do% {
    
    # Run model
    res <- glm(y ~ cb + dow + ns(date, df = 7 * length(unique(year))), 
      dat, family = quasipoisson)
    
    # Reduce coefficients to overall cumulative
    redall <- crossreduce(cb, res, cen = median(dat$tmean, na.rm = T))
    
    # Output
    coefs <- coef(redall); names(coefs) <- sprintf("coef%i", 1:nc)
    vcovs <- vechMat(vcov(redall))
    names(vcovs) <- sprintf("vcov%i%i", 
      row(vcov(redall))[lower.tri(row(vcov(redall)), diag = T)],
      col(vcov(redall))[lower.tri(col(vcov(redall)), diag = T)])
    data.frame(t(coefs), t(vcovs), 
      totdeath = sum(y, na.rm = T), conv = res$converged, agegroup = agegr)
  }
  
  # Add info about city and return
  mutate(out, city_code = city)
}

# Stop parallel
stopCluster(cl)

# Add to metadata
metadf <- merge(metadf, stage1res, all.x = T)
