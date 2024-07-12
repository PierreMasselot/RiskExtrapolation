################################################################################
#
#  RiskExtrapolation
#
#  Predictive meta-regression model
#
################################################################################

#--------------------------
# Second stage model
#--------------------------

# Fixed effect formula: Non linear effect of age + composite indices of vulnerability
# (see Appendices 1 and 2 for selection)
# NB: Additional terms such as a spatial background can be added
fixform <- sprintf("cbind(%s) ~ ns(age, knots = 60) + %s", 
  paste(coefvars, collapse = ", "), 
  paste(colnames(comps)[1:npc], collapse = " + "))

# Random effect formula: random effect of city
ranform <- ~ 1|city_code

# Extract variance-covariance 
smat <- select(metadf, starts_with("vcov")) |> data.matrix()

# Fit model
stage2res <- mixmeta(as.formula(fixform), random = ranform, 
  data = metadf, S = smat, subset = conv) 

####################################
# Plots
####################################

#------------------
# Age curves
#------------------

#----- Predict curves for several ages

# Create prediction data.frame
agedf <- cbind(data.frame(age = agebreaks), matrix(0, length(agebreaks), npc, 
  dimnames = list(NULL, sprintf("Comp%i", 1:npc))))

# Predict coefficients by age
agepreds <- predict(stage2res, agedf, vcov = T)

# Determine common MMT by age
uncentred <- mmtbasis %*% sapply(agepreds, "[[", "fit")
mmts <- mmtper[apply(uncentred, 2, which.min)]
mmt <- median(mmts)

# Construct centred ERFs by age
ageERF <- lapply(agepreds, function(x){
  crosspred(ovbasis, coef = x$fit, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

#----- Plot age-group specific curves

# Color palette
agepal <- mako(length(ageERF), end = .8, direction = -1)

# Initialize plot, draw grid and cusotm x-axis
plot(NA, bty = "l", xaxt = "n", 
  xlab = "Temperature percentile", ylab = "RR",
  xlim = range(ovper), ylim = c(.8, 2))
abline(v = ovaxis, h = axTicks(2), lty = 2, col = "lightgrey")
axis(1, at = ovaxis, labels = axisper)

# Add age curves
for (i in seq_along(ageERF)){
  lines(ageERF[[i]], ptype = "overall", col = agepal[i], ci = "area", 
    lwd = 2, ci.arg = list(col = adjustcolor(agepal[i], .2)))
}

# Overlay the RR=1 line
abline(h = 1)

# Add legend
legpars <- list(legend = agebreaks, col = agepal, cex = .8,
  lty = 1, lwd = 2, title = "Age", horiz = T, xpd = T)
legdim <- do.call(legend, c(legpars, list(x = "center", plot = F)))
do.call(legend, c(legpars, 
  list(x = mean(par("usr")[1:2]) - legdim$rect$w / 2, 
    y = par("usr")[4])))

# Save
dev.print(pdf, file = "figures/Fig3_age.pdf", width = 7, height = 5)


#------------------
# Evolution of curve for a city
#------------------

# Parameters on what city/age group combination to display
citypred <- "IT002C"
agecurve <- agelabs[4]

#----- Prepare ERFs

# Index for city
ind <- with(metadf, city_code == citypred & agegroup == agecurve)

# Fit models with various numbers of components and predict
comp_pred <- foreach(k = 0:npc) %do% {
  
  # Deal with no component and extract right number
  if (k == 0){
    kform <- sprintf("cbind(%s) ~ ns(age, knots = 60)", 
      paste(coefvars, collapse = ", "))
  } else {
    kform <- sprintf("cbind(%s) ~ ns(age, knots = 60) + %s", 
      paste(coefvars, collapse = ", "), 
      paste(colnames(comps)[1:k], collapse = " + "))
  }
  
  # Fit model
  fit <- mixmeta(as.formula(kform), random = ranform, 
    data = metadf, S = smat, subset = conv) 
  
  # Predict for one city
  predict(fit, metadf[ind,], vcov = T)
}

# Create ERFs
comp_erf <- lapply(comp_pred, function(x){
  firstpred <- mmtbasis %*% x$fit
  mmt <- mmtper[which.min(firstpred)]
  crosspred(ovbasis, coef = x$fit, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

# ERF for first-stage
coefsel <- unlist(metadf[ind, grep("coef", names(metadf))])
vcovsel <- metadf[ind, grep("vcov", names(metadf))] |> xpndMat()
firstpred <- mmtbasis %*% coefsel
mmt <- mmtper[which.min(firstpred)]
firsterf <- crosspred(ovbasis, coef = coefsel, vcov = vcovsel, 
  model.link = "log", at = ovper, cen = mmt)

#----- Plot

# Color palette
comppal <- mako(length(comp_erf), end = .8)

# Initialize plot, draw grid and custom x-axis
plot(NA, bty = "l", xaxt = "n", 
  xlab = "Temperature percentile", ylab = "RR",
  xlim = range(ovper), ylim = c(.8, 2))
abline(v = ovaxis, h = axTicks(2), lty = 2, col = "lightgrey")
axis(1, at = ovaxis, labels = axisper)

# Add age curves
for (i in seq_along(comp_erf)){
  lines(comp_erf[[i]], ptype = "overall", col = comppal[i], ci = "n", lwd = 2)
}

# Add first stage curve
lines(firsterf, ptype = "overall", col = 2, ci = "n", lty = 2, lwd = 3)

# Overlay the RR=1 line
abline(h = 1)

# Add legend
legdim <- legend("top", legend = 0:npc, col = comppal, cex = .8,
  lty = 1, lwd = 2, title = "# Components", horiz = T, xpd = T,
  box.col = "white")
legend(x = mean(par("usr")[1:2]), y = par("usr")[4] - legdim$rect$h, 
  legend = "First-stage", col = 2, lwd = 2, lty = 2, cex = .8,
  box.col = "white", xjust = 0.5)

# Save
dev.print(pdf, file = "figures/Fig4_compsPred.pdf")
