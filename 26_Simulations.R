################################################################################
#
#   New extensions
#   6. Simulations for uncertainty
#
################################################################################

#----------------------------
# Simulations from 2nd stage
#----------------------------

#----- Simulate coefficients for each city age

set.seed(12345)

# Simulate meta-coefficients from multivariate normal distribution
metacoefsim <- mvrnorm(nsim, coef(stage2res), vcov(stage2res))

# Design matrix from second stage
design_mat <- model.matrix(delete.response(terms(stage2res)), stage2df) %x% 
  diag(nc)

# Generate fixed part of city-age coefficients
fixsim <- design_mat %*% t(metacoefsim)

# Generate random part from
# nearPD forces vcov matrix to be psotivie definite
ransim <- lapply(ranpred, function(x) mvrnorm(nsim, x$fit, nearPD(x$vcov)$mat))
ransim <- t(do.call(cbind, ransim))

# Add them and rearrange as array
coefsim <- fixsim + ransim
coefsim <- aperm(array(c(coefsim), dim = c(nc, n * nage, nsim)), c(2, 1, 3))

#----- Apply attribution to all simulated

afsim <- foreach(sim = iapply(coefsim, 1), dat = dlist[cityageres$city]) %do% {
    
    # Estimate MMT
    bvar <- onebasis(dat$tmean, fun = varfun, degree = vardegree, 
      knots = quantile(dat$tmean, varper / 100))
    erf <- bvar %*% sim
    inrange <- between(dat$tmean, quantile(dat$tmean, mmprange[1] / 100), 
      quantile(dat$tmean, mmprange[2] / 100))
    immt <- dat$tmean[inrange][apply(erf[inrange,], 2, which.min)]
    
    # Centre basis
    cenmat <- onebasis(immt, fun = varfun, degree = vardegree, 
      knots = quantile(dat$tmean, varper / 100), 
      Boundary.knots = range(dat$tmean))
    
    # Compute daily AF and AN
    afdaysim <- mapply(function(cen, coef) (1 - exp(-scale(bvar, center = cen, 
        scale = F) %*% coef)),
      as.data.frame(t(cenmat)), as.data.frame(sim))
    
    # Indicator of heat days
    heatindsim <- outer(dat$tmean, immt, ">=")
    
    # Sum all and return
    cbind(total = colSums(afdaysim), 
      cold = colSums(afdaysim * (!heatindsim)), 
      heat = colSums(afdaysim * heatindsim))
  }

#----- Compute death rates for each simulation

# Compute associated AN
ansim <- Map("*", afsim, cityageres$death / nts)

# Compute excess death rates
excess_sim <- Map("/", ansim, cityageres$agepop / byrate)

# Compute standardized excess death rates
stdexcess_sim <- tapply(excess_sim, cityageres$city, function(x){
  Reduce("+", Map("*", x, isp)) / sum(isp)
})

#----- Compute eCI

# City-age AN
anCI <- sapply(ansim, apply, 2, quantile, c(.025, .975), simplify = "array")

# City-age excess death rates
excessCI <- sapply(excess_sim, apply, 2, quantile, c(.025, .975), 
  simplify = "array")

# City standardized excess rates
stdexcessCI <- sapply(stdexcess_sim, apply, 2, quantile, c(.025, .975), 
  simplify = "array")

#----- Compute region level eCI

# Aggregate AN by region
anregion <- tapply(ansim, cityageres[, c("agegroup", "geozone")],
  Reduce, f = "+")

# Compute excess death rates
popcountry <- aggregate(agepop ~ agegroup + geozone, data = cityageres, sum)
excess_region <- Map("/", anregion, popcountry$agepop / byrate)

# Compute standardized excess death rates
stdexcess_sim_region <- tapply(excess_region, popcountry$geozone, 
  function(x) Reduce("+", Map("*", x, isp)) / sum(isp)
)

# Compute eCI
stdexcessCI_region <- sapply(stdexcess_sim_region, apply, 2, quantile, 
  c(.025, .975), simplify = "array")
dimnames(stdexcessCI_region)[[1]] <- c("low", "high")

#----------------------------
# Simulations from predicted coef as comparison
#----------------------------

afsim_old <- foreach(ifix = newcoefs, iran = ranpred, 
  dat = dlist[cityageres$city]) %do% {
  
  # Simulate coefficients
  fixsim <- mvrnorm(nsim, ifix$fit, ifix$vcov)
  ransim <- mvrnorm(nsim, iran$fit, nearPD(iran$vcov)$mat)
  sim <- fixsim + ransim

  # Estimate MMT
  bvar <- onebasis(dat$tmean, fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100))
  erf <- bvar %*% t(sim)
  inrange <- between(dat$tmean, quantile(dat$tmean, mmprange[1] / 100), 
    quantile(dat$tmean, mmprange[2] / 100))
  immt <- dat$tmean[inrange][apply(erf[inrange,], 2, which.min)]
  
  # Centre basis
  cenmat <- onebasis(immt, fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100), 
    Boundary.knots = range(dat$tmean))
  
  # Compute daily AF and AN
  afdaysim <- mapply(function(cen, coef) (1 - exp(-scale(bvar, center = cen, 
    scale = F) %*% coef)),
    as.data.frame(t(cenmat)), as.data.frame(t(sim)))
  
  # Indicator of heat days
  heatindsim <- outer(dat$tmean, immt, ">=")
  
  # Sum all and return
  cbind(total = colSums(afdaysim), 
    cold = colSums(afdaysim * (!heatindsim)), 
    heat = colSums(afdaysim * heatindsim))
}

#----- Compute region level death rates

# Compute associated AN
ansim_old <- Map("*", afsim_old, cityageres$death / nts)

# Aggregate by region
anregion_old <- tapply(ansim_old, cityageres[, c("agegroup", "geozone")],
  Reduce, f = "+", simplify = F)

# Compute excess death rates
popcountry <- aggregate(agepop ~ agegroup + geozone, data = cityageres, sum)
excess_region_old <- Map("/", anregion_old, popcountry$agepop / byrate)

# Compute standardized excess death rates
stdexcess_region_sim_old <- tapply(excess_region_old, popcountry$geozone, 
  function(x) Reduce("+", Map("*", x, isp)) / sum(isp)
)

# Compute eCI
stdexcessCI_region_old <- sapply(stdexcess_region_sim_old, apply, 2, quantile, 
  c(.025, .975), simplify = "array")
dimnames(stdexcessCI_region_old)[[1]] <- c("low", "high")

#----------------------------
# Plots
#----------------------------

#----- Gradient interval for the distribution of std excess at regional level

# Prepare data.frame for simulated excesses
excess_dist <- data.frame(do.call(rbind, stdexcess_sim_region), 
  geozone = rep(names(stdexcess_sim_region), each = nsim), type = "new")
excess_dist_old <- data.frame(do.call(rbind, stdexcess_region_sim_old), 
  geozone = rep(names(stdexcess_region_sim_old), each = nsim), type = "old")
simdf <- rbind(excess_dist, excess_dist_old)

# Prepare estimated and CI data.frame
estdf <- melt(stdexcess_region, id = "geozone", 
  measure = names(stdexcess_region)[-1], value.name = "est")
cidfnew <- dcast(
  melt(stdexcessCI_region, varnames = c("ci", "variable", "geozone")),
  variable + geozone ~ ci)
cidfold <- dcast(
  melt(stdexcessCI_region_old, varnames = c("ci", "variable", "geozone")),
  variable + geozone ~ ci)
alldf <- cbind(rbind(cidfnew, cidfold), type = rep(c("new", "old"), each = 15))
alldf <- merge(alldf, estdf)

# Plot
ggplot(simdf) + theme_classic() + 
  stat_gradientinterval(aes(x = geozone, y = total, fill = type), 
    position = "dodge", point_colour = NA, .width = .95) + 
  scale_fill_scico_d(palette = "cork", name = "Sampling", direction = -1,
    labels = c(new = "Meta-regression", old = "Predictions")) + 
  geom_point(aes(x = geozone, y = est, group = type), 
    data = subset(alldf, variable == "total"), size = 4,
    position = position_dodge2(width = 1)) + 
  labs(x = "Region", y = "Standardized excess mortaltiy rates (x 100,000)")

# Save
ggsave("figures/Fig6a_IntervalsComparison.png")
