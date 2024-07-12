################################################################################
#
#  RiskExtrapolation
#
#  Uncertainty assessment
#
################################################################################

#----------------------------
# Simulations from 2nd stage
#----------------------------

#----- Simulate coefficients for each city age

# Simulate meta-coefficients from multivariate normal distribution
set.seed(12345)
metacoefsim <- mvrnorm(nsim, coef(stage2res), vcov(stage2res))

# Design matrix from second stage, expand for the multivariate outcome
design_mat <- model.matrix(delete.response(terms(stage2res)), metadf) %x% 
  diag(nc)

# Generate fixed part of city-age coefficients
fixsim <- design_mat %*% t(metacoefsim)

# Generate random part from
# nearPD forces vcov matrix to be positive definite
ransim <- foreach(co = iter(rancoef, "row"), v = iter(ranvcov, "row"), 
  .combine = rbind) %do% 
{
  t(mvrnorm(nsim, unlist(co), nearPD(xpndMat(v))$mat))
}

# Add them and rearrange as array
coefsim <- fixsim + ransim
coefsim <- aperm(array(c(coefsim), dim = c(nc, na, nsim)), c(2, 1, 3))

#----- Apply attribution to all simulated
# NB: Can be combined with the computation in script 27 
#   for computational efficiency
simres <- foreach(ires = iter(cityageres, by = "row"), sim = iapply(coefsim, 1),
    .final = rbindlist) %do%
{
    
  # Extract tmean series
  tmean <- tempdata[city_code == ires$city_code,]$tmean

  # Basis
  bvar <- onebasis(tmean, fun = varfun, degree = vardegree, 
    knots = quantile(tmean, varper / 100))
  cenvec <- onebasis(ires$mmt, fun = varfun, degree = vardegree, 
    knots = quantile(tmean, varper / 100), Bound = range(tmean))
  bvarcen <- scale(bvar, center = cenvec, scale = F)
  
  # Compute daily AN
  anday <- (1 - exp(-bvarcen %*% sim)) * ires$death / length(tmean)

  # Indicator of heat days
  heatind <- tmean >= ires$mmt
  
  # Sum all and return, including the MMT
  cbind(ires[, c("city_code", "agegroup", "pop", "geozone")], sim = 1:nsim, 
    an_total = colSums(anday), an_cold = colSums(anday[!heatind,]), 
    an_heat = colSums(anday[heatind,]))
}

# Compute excess death rates
simres <- mutate(simres, across(starts_with("an_"), 
  ~ byrate * .x / pop, .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Compute CIs
cityageci <- summarise(simres, 
  across(starts_with("an_") | starts_with("excessrate"), 
    list(low = ~ quantile(.x, .025), high = ~ quantile(.x, .975))),
  .by = c(city_code, agegroup))
cityageres <- merge(cityageres, cityageci)

#----- Compute city level results

citysim <- merge(simres, isp) |>
  # Aggregate impacts
  summarise(
    # Attributable numbers: sum age specific ANs
    across(starts_with("an_"), sum), pop = sum(pop),
    # Compute standardised excess death rates
    across(starts_with("excessrate"), weighted.mean, w = ispw),
    .by = c(city_code, sim)) |>
  rename_with(~ gsub("excess", "std", .x)) |>
  # Non standardised 
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Compute CIs
cityci <- summarise(citysim, 
  across(starts_with("an_") | starts_with("excessrate"), 
    list(low = ~ quantile(.x, .025), high = ~ quantile(.x, .975))),
  .by = city_code)
cityres <- merge(cityres, cityci)

#----- Compute regional level result

# Aggregate by region, starting by summing AN and pop
regionagesim <- summarise(simres, 
  across(starts_with("an_"), sum), pop = sum(pop),
    .by = c(geozone, agegroup, sim)) |>
  # Compute excess death rates
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Compute CIs
regionageci <- summarise(regionagesim, 
  across(starts_with("an_") | starts_with("excessrate"), 
    list(low = ~ quantile(.x, .025), high = ~ quantile(.x, .975))),
  .by = geozone)
regionageres <- merge(regionageres, regionageci)

# Aggregate across age levels, starting by merging standard pop weights
regionsim <-  merge(regionagesim, isp) |>
  # Compute total AN and std death rates
  summarise(across(starts_with("an_"), sum), pop = sum(pop),
    across(starts_with("excessrate"), weighted.mean, w = ispw),
    .by = c(geozone, sim)) |>
  rename_with(~ gsub("excess", "std", .x)) |>
  # Non-standardised rates for comparison
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Compute CIs
regionci <- summarise(regionsim, 
  across(starts_with("an_") | starts_with("excessrate") | 
      starts_with("stdrate"), 
    list(low = ~ quantile(.x, .025), high = ~ quantile(.x, .975))),
  .by = geozone)
regionres <- merge(regionres, regionci)

#----------------------------
# Simulations from predicted coef as comparison
#----------------------------

benchsim <- foreach(ico = iter(predcoefs, "row"), iv = iter(predvcov, "row"), 
  ires = iter(cityageres, by = "row"), .combine = rbind) %do% 
{
  
  # Simulate coefficients
  sim <- mvrnorm(nsim, unlist(ico), xpndMat(iv))

  # Extract tmean series
  tmean <- tempdata[city_code == ires$city_code,]$tmean
  
  # Basis
  bvar <- onebasis(tmean, fun = varfun, degree = vardegree, 
    knots = quantile(tmean, varper / 100))
  cenvec <- onebasis(ires$mmt, fun = varfun, degree = vardegree, 
    knots = quantile(tmean, varper / 100), Bound = range(tmean))
  bvarcen <- scale(bvar, center = cenvec, scale = F)
  
  # Compute daily AN
  anday <- (1 - exp(-bvarcen %*% t(sim))) * ires$death / length(tmean)
  
  # Indicator of heat days
  heatind <- tmean >= ires$mmt
  
  # Sum all and return, including the MMT
  cbind(ires[, c("city_code", "agegroup", "pop", "geozone")], sim = 1:nsim, 
    an_total = colSums(anday), an_cold = colSums(anday[!heatind,]), 
    an_heat = colSums(anday[heatind,]))
}

#----- Compute region level death rates

# Aggregate by region, starting by summing AN and pop
benchra <- summarise(benchsim, 
  across(starts_with("an_"), sum), pop = sum(pop),
    .by = c(geozone, agegroup, sim)) |>
  # Compute excess death rates
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Aggregate across age levels, starting by merging standard pop weights
benchreg <-  merge(benchra, isp) |>
  # Compute total AN and std death rates
  summarise(across(starts_with("an_"), sum), pop = sum(pop),
    across(starts_with("excessrate"), weighted.mean, w = ispw),
    .by = c(geozone, sim)) |>
  rename_with(~ gsub("excess", "std", .x)) |>
  # Non-standardised rates for comparison
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Compute CIs
benchci <- summarise(benchreg, across(
    starts_with("an_") | starts_with("excessrate") | starts_with("stdrate"), 
    list(low = ~ quantile(.x, .025), high = ~ quantile(.x, .975))),
  .by = geozone)
names(benchci)[-1] <- paste(names(benchci)[-1], "old", sep = "_")
regionres <- merge(regionres, benchci, by = "geozone")

#----------------------------
# Plots
#----------------------------

#----- Summarise heterogeneity of regions by Thompson and Higgins' H

# Extract coef & vcovs
resco <- select(cityageres, starts_with("coef"))
resvcov <- select(cityageres, starts_with("vcov"))

# Split by region and loop
hetdf <- foreach(coefs = isplit(resco, cityageres$geozone), 
  vcovs = isplit(resvcov, cityageres$geozone), .combine = rbind) %do% 
{
  # Transform coefs and vcovs appropriately
  geo <- coefs$key |> unlist(use.names = F)
  coefs <- as.data.frame(t(coefs$value))
  vcovs <- apply(vcovs$value, 1, xpndMat, simplify = F)
  
  # Compute weighted average of region
  wsum <- mapply(function(b, v) solve(v) %*% b, coefs, vcovs) |> rowSums()
  vmu <- lapply(vcovs, solve) |> Reduce(f = "+") |> solve()
  wmu <- vmu %*% wsum
    
  # Compute H
  qvec <- mapply(function(b, v) t(b - wmu) %*% solve(v) %*% (b - wmu),
    coefs, vcovs)
  h <- sum(qvec) / (length(coefs) - 1)
  
  # number of cities
  ncities <- length(coefs) / length(agelabs)
  
  # Return
  data.frame(geozone = geo, ncity = ncities, h2 = h)
}

#----- Gradient interval for the distribution of std excess at regional level

# Put together the simulations
simdf <- rbind(
  subset(regionsim, select = c(geozone, sim, stdrate_total)) |> 
    mutate(type = "meta"), 
  subset(benchreg, select = c(geozone, sim, stdrate_total)) |> 
    mutate(type = "pred"))

# Data.frame for the point estimate
ptdf <- merge(regionres, 
  expand.grid(geozone = regionres$geozone, type = unique(simdf$type)))

# Plot
ggplot(simdf) + theme_classic() + 
  stat_interval(aes(x = geozone, y = stdrate_total, color = type, 
    color_ramp = after_stat(level)), 
    position = position_dodge(width = .5), point_colour = NA, 
    .width = c(.5, .8, .95, 1)) + 
  scale_color_scico_d(palette = "cork", name = "Sampling", direction = -1,
    labels = c(meta = "Meta-regression", pred = "Predictions")) +
  geom_point(aes(x = geozone, y = stdrate_total, group = type), 
    data = ptdf, size = 5, position = position_dodge2(width = .5)) + 
  geom_label(aes(x = geozone, 
    label = sprintf("%i cities \n H: %2.2f", ncity, h2)), 
    y = max(simdf$stdrate_total), data = hetdf, vjust = .5, label.size = 0) + 
  labs(x = "Region", y = "Standardized excess mortality rates (x 100,000)",
    color_ramp = "Level")

# Save
ggsave("figures/Fig8_IntervalsComparison.pdf")
