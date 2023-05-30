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
blupext <- blup(stage2res, type = "residual")

# Remove duplicates (same for all age groups within location)
dups <- duplicated(blupext)
cityblup <- blupext[!dups,]

# Get city lat long
coords <- subset(metadata, obs, c("lon", "lat"))

# Create spatial object
blupgeo <- st_as_sf(cbind(cityblup, coords), coords = c("lon","lat"), 
  crs = st_crs(4326))

#---------------------------
# Fit spatial interpolation model by Kriging
#---------------------------

# Create gstat object: add all coefficients in the object for co-kriging
cokrig <- NULL
for (i in seq_len(ncol(coefs))) {
  form <- sprintf("b%i ~ 1", i)
  cokrig <- gstat(cokrig, id = sprintf("b%i", i), formula = as.formula(form), 
    data = blupgeo, set = list(nocheck = 1))
}

# Empirical semi-variogram
vario <- variogram(cokrig, cutoff = 300)

# Parameters for variogram model
variopars <- list(
  model = "Gau",
  nugget = NA,
  range = 100
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

# Extract map for Italy
italymap <- gisco_get_countries(year = "2020", epsg = "4326", country = "IT")

# Create grid
createraster <- raster(extent(italymap), nrow = 100, ncol = 100, 
  crs = st_crs(blupgeo))
rastcoord <- st_as_sf(as.data.frame(coordinates(createraster)), 
  coords = c("x","y"), crs = st_crs(blupgeo))

# Keep only land part
itaraster <- st_intersection(rastcoord, st_union(italymap))

# Predict on each cell
surfpred <- predict(spatmod, itaraster)

#----- Extract locations

# Location dataset
loc_sf <- st_as_sf(stage2df, coords = c("lon","lat"), crs = st_crs(4326))
loc_sf <- subset(loc_sf, !duplicated(city))
loc_sf$obs <- factor(loc_sf$obs, labels = c("Unobserved", "Observed"))

# Select spline coefficient to display
coefdisp <- 2

# Include blups for observed cities
loc_sf$blupres <- NA
loc_sf[loc_sf$obs == "Observed", "blupres"] <- cityblup[, coefdisp]

# Prepare surface
surfpred <- mutate(surfpred, pred = .data[[sprintf("b%i.pred", coefdisp)]],
  var = .data[[sprintf("b%i.var", coefdisp)]])

# Shapes
shp <- c(Observed = 21, Unobserved = 15)

#----- Interpolation plot

# Limits of random coefficients
lims <- range(c(surfpred$pred, loc_sf$blupres), na.rm = T)

# Plot
ggpred <- ggplot() + theme_void() + 
  geom_sf(aes(col = pred), data = surfpred, pch = 15) + 
  geom_sf(data = italymap, fill = NA, inherit.aes = F, col = grey(.5)) + 
  geom_sf(aes(fill = blupres, shape = obs), data = loc_sf, size = 3, 
    col = grey(.3)) + 
  scale_colour_scico(palette = "cork", midpoint = 0,
    name = "Random coefficient", limits = lims) + 
  scale_fill_scico(palette = "cork", midpoint = 0,
    name = "Random coefficient", limits = lims) + 
  scale_shape_manual(values = shp, name = "")

#----- Variance plot

# Plot
ggvar <- ggplot() + theme_void() + 
  geom_sf(aes(col = var), data = surfpred, pch = 15) + 
  geom_sf(data = italymap, fill = NA, inherit.aes = F, col = grey(.5)) + 
  geom_sf(aes(fill = obs, shape = obs), data = loc_sf, size = 3, 
    col = grey(.3)) + 
  scale_colour_scico(palette = "oslo", 
    name = "Variance", direction = -1) + 
  scale_fill_manual(values = c(Observed = "white", Unobserved = "grey50"), 
    name = "", guide = "none") + 
  scale_shape_manual(values = shp, name = "", guide = "none")

#----- Put together and save

wrap_plots(ggpred, ggvar, guides = "collect")

# Save
ggsave("figures/Fig4_extrapol_maps.pdf")
