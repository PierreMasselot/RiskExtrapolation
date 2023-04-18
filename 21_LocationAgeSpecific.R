################################################################################
#
#  RiskExtrapolation
#   First-stage analysis
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
stage1res <- foreach(i = seq_len(nrow(metadata)), 
  .packages = c("dlnm", "splines")) %dopar% 
{
  
  dat <- dlist[[i]]  
  
  # Define crossbasis
  argvar <- list(fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100, na.rm = T))
  cb <- crossbasis(dat$tmean, lag = maxlag, argvar = argvar,
    arglag = list(fun = lagfun, knots = lagknots))
  
  # Prepare output
  out <- list()
  
  # Loop on ages
  for (a in agelabs){
    
    # Define outcome
    y <- dat[[sprintf("all_%s", a)]]
    
    # Run model
    res <- glm(y ~ cb + dow + ns(date, df = 7 * length(unique(year))), 
      dat, family = quasipoisson)
    
    # Reduce coefficients
    # Reduction to overall cumulative
    redall <- crossreduce(cb, res, cen = median(dat$tmean, na.rm = T))
    
    # Output
    out[[a]] <- list(coef = coef(redall), vcov = vcov(redall), 
      totdeath = sum(y, na.rm = T), conv = res$converged)
  }
  
  out
}

# Stop parallel
stopCluster(cl)

# Rename
names(stage1res) <- metadata$CITY_CODE

#-------------------------
# Extract coefficients and vcovs
#-------------------------

# Get all results in big list
allresults <- unlist(stage1res, recursive = F)

# Get coefs and vcov
coefs <- t(sapply(allresults, "[[", "coef"))
vcovs <- lapply(allresults, "[[", "vcov")

#----- Separate into observed and prediction for the illustration

# Extract observed locations
coefs_obs <- subset(coefs, stage2df$obs)
vcovs_obs <- vcovs[stage2df$obs]

# And 'unobserved' locations for later comparison
coefs_pred <- subset(coefs, !stage2df$obs)
vcovs_pred <- vcovs[!stage2df$obs]