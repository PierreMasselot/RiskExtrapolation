################################################################################
#
#  RiskExtrapolation
#
#  Extrapolation to unobserved locations
#
################################################################################

#---------------------------
# Predict curves from meta-regression and kriging models
#---------------------------

#----- Fixed effect prediction (age and vulnerability indices)

# Predict coefficients
fixpred <- predict(stage2res, metadf, vcov = T)

# Extract coefs and vcovs
fixcoef <- sapply(fixpred, "[[", "fit") |> t()
fixvcov <- sapply(fixpred, function(x) vechMat(x$vcov)) |> t()

#----- Random effects from Kriging model

# Predict at every city
spatpred <- predict(spatmod, metadf)

# Extract coefficients and covariance matrix
rancoef <- dplyr::select(spatpred, matches("pred$"))
ranvcov <- dplyr::select(spatpred, contains("var") | contains("cov")) 
ranvcov <- ranvcov[, order(gsub("(.var)|(cov.)", "", names(ranvcov)))]

#----- Create curves

# Add fixed and random effect predictions
predcoefs <- fixcoef + rancoef
predvcov <- fixvcov + ranvcov

# Put together into data.frame
names(predcoefs) <- coefvars
names(predvcov) <- grep("vcov", names(stage1res), value = T)
cityageres <- cbind(metadf[, c("city_code", "agegroup", "obs", "geozone")], 
  predcoefs, predvcov)

# Reconstruct curves
predERF <- foreach(co = iter(predcoefs, by = "row"), 
  v = iter(predvcov, by = "row")) %do% 
{
  uncentred <- mmtbasis %*% unlist(co)
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = unlist(co), vcov = xpndMat(unlist(v)), 
    model.link = "log", at = ovper, cen = mmt)  
}


####################################
# Plots
####################################

#-------------------------------
# Show predicted curves
#-------------------------------

# Indices for highest age-group
set.seed(1)
curveind <- with(metadf, which(agegroup == levels(agegroup)[4] & !obs)) |>
  sample(10)

# Palette (depends on latitude)
latitudes <- metadf[curveind, "lat"]
latbreaks <- seq(floor(min(latitudes)), ceiling(max(latitudes)), by = 1)
latcut <- cut(latitudes, breaks = latbreaks)
predpal <- mako(nlevels(latcut), end = .8)

# Initialize plot, draw grid and custom x-axis
par(mar = c(5, 4, 4, 7) + .1)
plot(NA, bty = "l", xaxt = "n", 
  xlab = "Temperature percentile", ylab = "RR",
  xlim = range(ovper), ylim = c(.95, 2.5))
abline(v = ovaxis, h = axTicks(2), lty = 2, col = "lightgrey")
axis(1, at = ovaxis, labels = axisper)

# Add curves
for (i in seq_along(curveind)){
  lines(predERF[[curveind[i]]], ptype = "overall", 
    col = predpal[latcut[i]], lwd = 2)
}

# Overlay the RR=1 line
abline(h = 1)

# Add legend
image.plot(legend.only = T, zlim = range(latbreaks), breaks = latbreaks, 
  col = predpal, legend.lab = "Latitude")

# Save
dev.print(pdf, file = "figures/Fig6_examplePred.pdf")

