################################################################################
#
#  Appendix 3.
#  Alternative ways to extrapolate random effects
#
################################################################################

# Need the cityblup and blupgeo objects created

#----- Inverse distance weighting

# Power parameter for IDW
idpgrid <- seq(.5, 2, by = .1)

# Loop on coefficients
spatmod <- NULL; bestidp <- rep(NA, nc)
for (i in 1:nc) {
  
  # Specify model
  idwlist <- lapply(idpgrid, function(idp){
    gstat(spatmod, id = sprintf("b%i", i), 
      formula = as.formula(sprintf("b%i ~ 1", i)), 
      data = blupgeo, set = list(idp = idp))
  })
  
  # Estimate prediction error by CV
  cve <- sapply(idwlist, function(x) sum(gstat.cv(x, verbose = F)$residual^2))
  
  # Keep best one
  spatmod <- idwlist[[which.min(cve)]]
  bestidp[i] <- idpgrid[which.min(cve)]
}