################################################################################
#
#   New extensions
#   6. Predictions at all locations, including unobserved
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

#----- Create curves

# Add to fixed and random effect predictions
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
curveind <- with(stage2df, which(agegroup == levels(agegroup)[4] & !obs))

# Palette (depends on latitude)
latcut <- cut(stage2df[curveind, "lat"], breaks = 10)
predpal <- mako(10, end = .8)[latcut]

# Initialize plot, draw grid and custom x-axis
par(mar = c(5, 4, 4, 7) + .1)
plot(NA, bty = "l", xaxt = "n", 
  xlab = "Temperature percentile", ylab = "RR",
  xlim = range(ovper), ylim = c(.95, 2.5))
abline(v = ovaxis, h = axTicks(2), lty = 2, col = "lightgrey")
axis(1, at = ovaxis, labels = axisper)

# Add curves
for (i in seq_along(curveind)){
  lines(predERF[[curveind[i]]], ptype = "overall", col = predpal[i], lwd = 2)
}

# Overlay the RR=1 line
abline(h = 1)

# Add legend
legend(x = par("usr")[2], y = par("usr")[4], legend = rev(levels(latcut)),
  fill = rev(mako(10, end = .8)), bty = "n", xpd = T, title = "Latitude")

# Save
dev.print(pdf, file = "figures/Fig5_examplePred.pdf")

#-------------------------------
# Prediction scatterplots
#-------------------------------

#----- Model with all cities to compare

# Fit model
fullform <- coefs ~ ns(age, knots = 60) + comps
refmod <- mixmeta(fullform, random = ranform, data = stage2df, S = vcovs) 

# Extract BLUP
refcoefs <- blup(refmod, vcov = T)

# BLUP curves
refERF <- lapply(refcoefs, function(x){
  uncentred <- mmtbasis %*% x$blup
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = x$blup, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

#----- Extract prediction at extreme percentiles

# Extract heat and cold for predicted
predRR <- sapply(predERF, "[[", "allRRfit")[predper %in% c(1, 99),]

# Extract heat and cold for BLUPs
refRR <- sapply(refERF, "[[", "allRRfit")[predper %in% c(1, 99),]

#----- Plot

# Color palette
pals <- list(
  scico(length(agebreaks), direction = -1, end = .8, palette = "devon"),
  scico(length(agebreaks), direction = 1, begin = .2, palette = "bilbao")
)

# Labels for percentile
perlabs <- c("Cold", "Heat")

# Loop on percentiles to display scatterplot
scatters <- lapply(seq_len(nrow(predRR)), function(i){
  df <- subset(
    data.frame(pred = predRR[i,], ref = refRR[i,], age = stage2df$agegroup),
    !stage2df$obs)
  
  ggplot(df) + theme_classic() + 
    geom_point(aes(x = ref, y = pred, fill = age), size = 3, shape = 21) + 
    geom_abline(slope = 1, linetype = 2, col = "grey", size = 1) + 
    scale_fill_manual(values = pals[[i]], guide = guide_legend(
      title.position = "left", label.position = "bottom"
    )) + 
    labs(x = "BLUP from full model", y = "Prediction", fill = "Age group") + 
    annotate("text", x = -Inf, y = Inf, hjust = -.1, vjust = 1, size = 5,
      label = sprintf("R2 = %2.0f %%", 
        summary(lm(pred ~ ref, df))$r.squared * 100)) + 
    ggtitle(sprintf("%s) %s", letters[i], perlabs[i])) + 
    lims(x = range(df[,1:2]), y = range(df[1:2])) + 
    theme(legend.position = "bottom", aspect.ratio = 1)
})

# Put together
wrap_plots(scatters, widths = 1)

# Save
ggsave("figures/Fig6_ObsPred.pdf", width = 10, height = 6.5)
