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

# Create prediction data.frame including all extracted components
newdata <- list(pcs = comps[,1:npc], age = stage2df$age)

# Predict coefficients
fixpred <- predict(stage2res, newdata, vcov = T)

#----- Random effects from Kriging model

# Predict at every city
spatpred <- predict(spatmod, stage2df)

# Extract prediction and covariance matrix
ranpred <- foreach(x = iter(spatpred, by = "row")) %do% {
  fit <- dplyr::select(x, contains("pred")) |> unlist()
  vcov <- dplyr::select(x, contains("var") | contains("cov")) |> unlist()
  vcov <- vcov[order(gsub("(.var)|(cov.)", "", names(vcov)))]
  list(fit = fit, vcov = xpndMat(vcov))
}
names(ranpred) <- stage2df$city

#----- Create curves

# Add fixed and random effect predictions
predcoefs <- Map(function(fix, ran) Map(`+`, fix, ran),
  fixpred, ranpred)

# Reconstruct curves
predERF <- lapply(predcoefs, function(x){
  uncentred <- mmtbasis %*% x$fit
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = x$fit, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

####################################
# Plots
####################################

#-------------------------------
# Show predicted curves
#-------------------------------

# Indices for highest age-group
set.seed(1)
curveind <- with(stage2df, which(agegroup == levels(agegroup)[4] & !obs)) |>
  sample(10)

# Palette (depends on latitude)
latitudes <- stage2df[curveind, "lat"]
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
dev.print(pdf, file = "figures/Fig5_examplePred.pdf")

