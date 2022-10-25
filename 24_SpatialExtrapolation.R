################################################################################
#
#   New extensions
#   4. Spatialization
#
################################################################################

#---------------------------
# Prepare BLUP spatial information
#---------------------------

# Extract random effect
blupext <- blup(stage2res, type = "residual")

# Remove duplicates
dups <- duplicated(blupext)
cityblup <- blupext[!dups,]

# Get city lat long
coords <- subset(metadata, obs, c("lon", "lat"))

# Create spatial object
blupgeo <- st_as_sf(cbind(cityblup, coords), coords = c("lon","lat"), 
  crs = st_crs(4326))

#---------------------------
# Fit spatial interpolation model
#---------------------------

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
  
#----- Kriging

# Create gstat object: add all coefficients in the object for co-kriging
cokrig <- NULL
for (i in seq_len(nc)) {
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
spatmod <- fit.lmc(vario, cokrig, varmod,
  fit.method = 6, correct.diagonal = 1.01)

#---------------------------
# Predict
#---------------------------

#----- At locations (including unobserved ones)

# Predict at every city
spatpred <- predict(spatmod, as_Spatial(
  st_as_sf(stage2df, coords = c("lon","lat"), crs = st_crs(4326))))

# Extract prediction and covariance matrix
ranpred <- apply(data.matrix(spatpred@data), 1, function(x){
  fit <- x[1:nc * 2 - 1]
  vcov <- matrix(NA, nc, nc)
  diag(vcov) <- x[1:nc * 2]
  vcov[upper.tri(vcov)] <- x[-(1:(nc * 2))]
  vcov[lower.tri(vcov)] <- t(vcov)[lower.tri(vcov)]
  vcov[is.na(vcov)] <- 0 # To deal with IDW
  list(fit = fit, vcov = vcov)
})
names(ranpred) <- unique(stage2df$city)

# Add to fixed effect prediction
predcoefs <- Map(function(fix, ran) Map(`+`, fix, ran),
  newcoefs, ranpred)

#----- On a grid

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


#---------------------------
# Map interpolation
#---------------------------

# Location dataset
loc_sf <- st_as_sf(cbind(spatpred@data, stage2df),
  coords = c("lon","lat"), crs = st_crs(4326))
loc_sf <- subset(loc_sf, !duplicated(city))
loc_sf$obs <- factor(loc_sf$obs, labels = c("Unobserved", "Observed"))

# Select spline coefficient to display
coefdisp <- 5
loc_sf <- mutate(loc_sf, pred = .data[[sprintf("b%i.pred", coefdisp)]],
  var = .data[[sprintf("b%i.var", coefdisp)]])

# Prepare surface
surfpred <- mutate(surfpred, pred = .data[[sprintf("b%i.pred", coefdisp)]],
  var = .data[[sprintf("b%i.var", coefdisp)]])

# Shapes
shp <- c(Observed = 21, Unobserved = 15)

#----- Interpolation plot

# Limits of random coefficients
lims <- range(c(surfpred$pred, loc_sf$pred))

# Plot
ggplot() + theme_void() + 
  geom_sf(aes(col = pred), data = surfpred, pch = 15) + 
  geom_sf(data = italymap, fill = NA, inherit.aes = F, col = grey(.5)) + 
  geom_sf(aes(fill = pred, shape = obs), data = loc_sf, size = 3, 
    col = grey(.3)) + 
  scale_colour_scico(palette = "cork", midpoint = 0,
    name = "Random coefficient", limits = lims) + 
  scale_fill_scico(palette = "cork", midpoint = 0,
    name = "Random coefficient", limits = lims) + 
  scale_shape_manual(values = shp, name = "")

# Save
ggsave("figures/Fig4_blupresMap.pdf")

#----- Variance plot

# Plot
ggplot() + theme_void() + 
  geom_sf(aes(col = var), data = surfpred, pch = 15) + 
  geom_sf(data = italymap, fill = NA, inherit.aes = F, col = grey(.5)) + 
  geom_sf(aes(fill = var, shape = obs), data = loc_sf, size = 3, 
    col = grey(.3)) + 
  scale_colour_scico(palette = "oslo", 
    name = "Variance", direction = -1) + 
  scale_fill_scico(palette = "oslo", 
    name = "Variance", direction = -1, guide = "none") + 
  scale_shape_manual(values = shp, name = "")

# Save
ggsave("figures/Fig4_blupVarMap.pdf")