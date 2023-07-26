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