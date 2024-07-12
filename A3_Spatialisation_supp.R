################################################################################
#
#  Appendix 3.
#  Spatialisation details
#
################################################################################

# Run after the 25_Spatialisation.R script

#------------------
# Display variogram
#------------------

# Description of parameters


# Plot
plot(vario, spatmod)

# Save
dev.print(pdf, "figures/FigS2_variogram.pdf")


# #------------------
# # Sensitivity analysis
# #------------------
# 
# # Grid of ranges
# rnggrid <- c(0.5, 1, 1.5, 2)
# 
# # Loop on grid to fit variogram
# rngres <- lapply(rnggrid, function(rng){
#   vp <- variopars; vp$range <- rng
#   vm <- do.call(vgm, variopars)
#   fit.lmc(vario, cokrig, varmod, fit.method = 6, correct.diagonal = 1.01)
# })
# 
# # Loop on results to get ERFs
# rngerf <- lapply(rngres, function(krigmod){
#   
#   # Predict
#   sp <- predict(krigmod, metadf)
#   rc<- dplyr::select(sp, matches("pred$"))
#   rv <- dplyr::select(sp, contains("var") | contains("cov")) 
#   rv <- rv[, order(gsub("(.var)|(cov.)", "", names(rv)))]
#   
#   # Full coef
#   pc <- fixcoef + rc
#   pv <- fixvcov + rv
#   
#   # ERF
#   uncentred <- mmtbasis %*% t(data.matrix(pc))
#   mmt <- mmtper[apply(uncentred, 2, which.min)]
#   cenb <- onebasis(mmt, fun = varfun, degree = vardegree, knots = ovknots, 
#     Bound = range(ovper)) 
#   mapply(function(co, ce) scale(ovbasis, ce) %*% co, as.data.frame(t(pc)),
#     as.data.frame(t(cenb)))
# })
# 
# # Get cold and heat
# coldrr <- sapply(rngerf, function(x){
#   x[predper == 1, metadf$obs]
# })



#------------------
# Inverse distance weighting
# Replaces lines 23-48 in 25_Spatialisation.R
#------------------

# Power parameter for IDW
idpgrid <- seq(.5, 2, by = .1)

# Loop on coefficients
spatmod <- NULL; bestidp <- rep(NA, nc)
for (i in 1:nc) {
  
  # Specify model
  idwlist <- lapply(idpgrid, function(idp){
    gstat(spatmod, id = sprintf("b%i", i), 
      formula = as.formula(sprintf("b%i ~ 1", i)), locations = ~ lon + lat,
      data = blupres, set = list(idp = idp))
  })
  
  # Estimate prediction error by CV
  cve <- sapply(idwlist, function(x) sum(gstat.cv(x, verbose = F)$residual^2))
  
  # Keep best one
  spatmod <- idwlist[[which.min(cve)]]
  bestidp[i] <- idpgrid[which.min(cve)]
}