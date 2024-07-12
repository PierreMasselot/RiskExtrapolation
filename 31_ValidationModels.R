################################################################################
#
#  RiskExtrapolation
#
#  First-stage unobserved ERFs for comparison
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
validres <- foreach(dat = isplit(tspred, tspred$city_code), .combine = rbind, 
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

#----- Create new full dataset for risk prediction

# Merge full dataset
metafull <- select(metadf, city_code, agegroup, lon, lat, age, obs,
    all_of(unname(metaprednames))) |>
  merge(rbind(stage1res, validres))

#-------------------------
# Create predictive model including "unobserved" locations
#-------------------------

#----- New composite indices of vulnerability

# Compute PLSR (basic PLS computed as a linear model)
fullpls <- plsr(as.formula(plsform), data = metafull, scale = T)

# Extract scores for all locations and age
newcomp <- predict(fullpls, newdata = metafull, type = "scores")
colnames(newcomp) <- sprintf("Comp%i", 1:ncol(newcomp))
metafull <- cbind(metafull, newcomp)

#----- New predictive model

# Extract variance-covariance 
sfull <- select(metafull, starts_with("vcov")) |> data.matrix()

# Fit model
validmeta <- mixmeta(as.formula(fixform), random = ranform, 
  data = metafull, S = sfull, subset = conv) 

# Extract BLUPs
validblups <- blup(validmeta, vcov = T)
blupcoefs <- sapply(validblups, "[[", "blup") |> t()
blupvcov <- sapply(validblups, function(x) vechMat(x$vcov)) |> t()
colnames(blupvcov) <- grep("vcov", names(stage1res), value = T)
validblups <- cbind(subset(metafull, conv, c(city_code, agegroup)), blupcoefs,
  blupvcov)


#-------------------------
# Alternative predictions
#-------------------------

#----- Fixed effect part only

# Rename and store into data.frame
names(fixcoef) <- coefvars
colnames(fixvcov) <- grep("vcov", names(stage1res), value = T)
fixedres <- cbind(metadf[, c("city_code", "agegroup")], 
  fixcoef, fixvcov)

#----- Age-only model

# Refit model
ageform <- sprintf("cbind(%s) ~ ns(age, knots = 60)", 
  paste(coefvars, collapse = ", "))
ageonlymeta <- mixmeta(as.formula(ageform), random = ranform, 
  data = metadf, S = smat, subset = conv) 

# Predict and store
ageonlypred <- predict(ageonlymeta, metadf, vcov = T)
ageonlycoef <- sapply(ageonlypred, "[[", "fit") |> t()
ageonlyvcov <- sapply(ageonlypred, function(x) vechMat(x$vcov)) |> t()
names(ageonlycoef) <- coefvars
colnames(ageonlyvcov) <- grep("vcov", names(stage1res), value = T)
ageonlyres <- cbind(metadf[, c("city_code", "agegroup")], 
  ageonlycoef, ageonlyvcov)

#----- Components only model

# Refit model
compform <- sprintf("cbind(%s) ~ %s", 
  paste(coefvars, collapse = ", "), 
  paste(colnames(comps)[1:npc], collapse = " + "))
compmeta <- mixmeta(as.formula(compform), random = ranform, 
  data = metadf, S = smat, subset = conv) 

# Predict and store
comppred <- predict(compmeta, metadf, vcov = T)
compcoef <- sapply(comppred, "[[", "fit") |> t()
compvcov <- sapply(comppred, function(x) vechMat(x$vcov)) |> t()
names(compcoef) <- coefvars
colnames(compvcov) <- grep("vcov", names(stage1res), value = T)
compres <- cbind(metadf[, c("city_code", "agegroup")], compcoef, compvcov)

#----- NULL model

# Refit model
nullform <- sprintf("cbind(%s) ~ 1", 
  paste(coefvars, collapse = ", "))
nullmeta <- mixmeta(as.formula(nullform), random = ranform, 
  data = metadf, S = smat, subset = conv) 

# Predict and store
nullpred <- predict(nullmeta, metadf, vcov = T)
nullcoef <- sapply(nullpred, "[[", "fit") |> t()
nullvcov <- sapply(nullpred, function(x) vechMat(x$vcov)) |> t()
names(nullcoef) <- coefvars
colnames(nullvcov) <- grep("vcov", names(stage1res), value = T)
nullres <- cbind(metadf[, c("city_code", "agegroup")], nullcoef, nullvcov)