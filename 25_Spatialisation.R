################################################################################
#
#  RiskExtrapolation
#
#  Spatialisation of risk
#
################################################################################

#---------------------------
# Prepare BLUP spatial information
#---------------------------

# Extract random effect
blupres <- blup(stage2res, type = "residual")

# Add city and lat/long
blupres <- cbind(blupres, stage2_obs[,c("city", "lon", "lat")])

# Remove duplicates (same for all age groups within location)
dups <- duplicated(blupres)
blupres <- blupres[!dups,]

#---------------------------
# Fit spatial interpolation model by Kriging
#---------------------------

# Create gstat object: add all coefficients in the object for co-kriging
cokrig <- NULL
for (i in seq_len(ncol(coefs))) {
  form <- sprintf("b%i ~ 1", i)
  cokrig <- gstat(cokrig, id = sprintf("b%i", i), formula = as.formula(form), 
    data = blupres, locations = ~ lon + lat)
}

# Empirical semi-variogram
vario <- variogram(cokrig)

# Parameters for variogram model
variopars <- list(
  model = "Gau",
  nugget = NA,
  range = 1
)

# Estimate variogram model
varmod <- do.call(vgm, variopars)
spatmod <- fit.lmc(vario, cokrig, varmod, fit.method = 6, 
  correct.diagonal = 1.01)

####################################
# Plots
####################################

#------------------
# Extrapolation of BLUPs
#------------------

#----- Kriged values on a grid covering Italy

# Create grid
itaraster <- rast(nrows = 100, ncols = 100, extent = ext(italymap), 
  crs = crs(italymap))

# Predict
surfpred <- interpolate(itaraster, spatmod, xyNames=c("lon", "lat"))

# Mask sea part
surfpred <- mask(surfpred, italymap)

#----- Extract locations

# Location dataset
loc_sf <- st_as_sf(stage2df, coords = c("lon","lat"), crs = st_crs(4326))
loc_sf <- subset(loc_sf, !duplicated(city))
loc_sf$obs <- factor(loc_sf$obs, labels = c("Unobserved", "Observed"))

# Include blups for observed cities
loc_sf <- merge(loc_sf, blupres, all.x = T)
loc_sf <- subset(loc_sf, select = -c(age, agegroup, lon, lat))

#----- Interpolation plot

# Select spline coefficient to display
coefdisp <- 1

# Shapes
shp <- c(Observed = 21, Unobserved = 15)

# Limits of random coefficients
lims <- range(c(values(surfpred)[, sprintf("b%i.pred", coefdisp)], 
  loc_sf[[sprintf("b%i", coefdisp)]]), na.rm = T)

# Plot
ggpred <- ggplot() + theme_void() + 
  geom_spatraster(aes(fill = .data[[sprintf("b%i.pred", coefdisp)]]), 
    data = surfpred) + 
  geom_sf(data = italymap, fill = NA, inherit.aes = F, col = grey(.5)) + 
  geom_sf(aes(fill = .data[[sprintf("b%i", coefdisp)]], shape = obs), 
    data = loc_sf, size = 3, col = grey(.3)) + 
  scale_fill_scico(palette = "cork", midpoint = 0,
    name = "BLUP residual", limits = lims, na.value = "white") + 
  scale_shape_manual(values = shp, name = "")

#----- Variance plot

# Plot
ggvar <- ggplot() + theme_void() + 
  geom_spatraster(aes(fill = .data[[sprintf("b%i.var", coefdisp)]]), 
    data = surfpred) + 
  scale_fill_scico(palette = "oslo", name = "Variance", direction = -1,
    na.value = "white") + 
  new_scale("fill") +
  geom_sf(data = italymap, fill = NA, inherit.aes = F, col = grey(.5)) + 
  geom_sf(aes(fill = obs, shape = obs), data = loc_sf, size = 3, 
    col = grey(.3)) + 
  scale_fill_manual(values = c(Observed = "white", Unobserved = "grey50"), 
    name = "", guide = "none") + 
  scale_shape_manual(values = shp, name = "", guide = "none")

#----- Put together and save

wrap_plots(ggpred, ggvar, guides = "collect")

# Save
ggsave("figures/Fig4_extrapol_maps.pdf")
