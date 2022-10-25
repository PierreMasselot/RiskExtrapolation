################################################################################
#
#   New extensions
#   3. Prediction to unobserved locations
#
################################################################################

######################################
# PREDICTING ERF AT UNOBSERVED LOCATIONS
######################################

#--------------------------
# Second stage model
#--------------------------

#----- Dimension reduction

# Number of components
kcomp <- 3

# Apply PLS
plsres <- plsr(coefs ~ metapred, scale = T, subset = stage2df$obs)

# Extract scores
comps <- predict(plsres, newdata = metapred, type = "scores",
  ncomp = 1:kcomp)

#----- Apply model

# Formula (see scripts 21 and 22 for selection)
fullform <- update(fixform, ~ . + ns(age, knots = 60) + comps)

# Fit model
stage2res <- mixmeta(fullform, random = ranform, data = stage2df, S = vcovs,
  subset = obs) 

#--------------------------
# Prediction
#--------------------------

#----- Predict at all locations, including unobserved

# Add to prediction data
newdata <- list(comps = comps, age = stage2df$age)

# Predict coefficients
newcoefs <- predict(stage2res, newdata, vcov = T)

# Reconstruct curves
newERF <- lapply(newcoefs, function(x){
  uncentred <- mmtbasis %*% x$fit
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = x$fit, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

#--------------------------
# Comparison: model with all cities
#--------------------------

# A model with all cities
refmod <- mixmeta(fullform, random = ranform, data = stage2df, S = vcovs) 

# Extract BLUP
refcoefs <- predict(refmod, vcov = T)

# BLUP curves
refERF <- lapply(refcoefs, function(x){
  uncentred <- mmtbasis %*% x$fit
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = x$fit, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

####################################
# Plots
####################################

#-------------------------------
# Just show predicted
#-------------------------------

# Indices for highest age-group
curveind <- with(stage2df, which(agegroup == levels(agegroup)[4] & !obs))

# Palette
predpal <- scico(sum(!metadata$obs), palette = "batlow")

# Initialize plot, draw grid and custom x-axis
plot(NA, bty = "l", xaxt = "n", 
  xlab = "Temperature percentile", ylab = "RR",
  xlim = range(ovper), ylim = c(.95, 2))
abline(v = ovaxis, h = axTicks(2), lty = 2, col = "lightgrey")
axis(1, at = ovaxis, labels = axisper)

# Add curves
for (i in seq_along(curveind)){
  lines(newERF[[curveind[i]]], ptype = "overall", col = predpal[i], lwd = 2)
}

# Overlay the RR=1 line
abline(h = 1)

# Save
dev.print(pdf, file = "figures/Fig3b_examplePred.pdf")

#-------------------------------
# Prediction scatterplots
#-------------------------------

#----- Extract prediction at extreme percentiles

# First-stage ERF
firststageERF <- Map(function(fit, vcov){
  uncentred <- mmtbasis %*% fit
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = fit, vcov = vcov, model.link = "log", 
    at = ovper, cen = mmt)
}, as.data.frame(t(coefs)), vcovs)

# Extract heat and cold for first-stage
firststageRR <- sapply(firststageERF, "[[", "allRRfit")[predper %in% c(1, 99),]

# Extract heat and cold for predicted
predRR <- sapply(newERF, "[[", "allRRfit")[predper %in% c(1, 99),]

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
    labs(x = "Full model", y = "Prediction", fill = "Age group") + 
    annotate("text", x = min(df$ref), y = max(df$pred), hjust = 0, size = 5,
      label = sprintf("R2 = %2.0f %%", 
        summary(lm(pred ~ ref, df))$r.squared * 100)) + 
    ggtitle(sprintf("%s) %s", letters[i], perlabs[i])) + 
    lims(x = range(df[,1:2]), y = range(df[1:2])) + 
    theme(legend.position = "bottom", aspect.ratio = 1)
})

# Put together
wrap_plots(scatters, widths = 1)

# Save
ggsave("figures/Fig3_ObsPred.pdf", width = 10, height = 7)

#-------------------------------
# Comparison of prediction with observed curves
#-------------------------------

#----- Plot side by side for comparison

# Graphical parameters
pal <- viridis(2)

# Plot dimensions
nrows <- 6
ncols <- nlevels(stage2df$agegroup)

# Plot layout
pdf("figures/predictionComparison.pdf", width = 12, height = 15)
layout(matrix(seq(nrows * ncols), nrow = nrows, byrow = T))
par(mar = c(4, 3.8, 3, 2.4), mgp = c(2.5,1,0), las = 1)

# Loop on ERF
for (i in which(!stage2df$obs)){
  
  # Plot first-stage
  plot(firststageERF[[i]], xlab = "Temperature percentile (%)", 
    main = paste(stage2df$CITY_NAME[i], stage2df$agegroup[i]), 
    ylab = "RR", col = pal[1], lwd = 2, ylim = c(.5, 3), 
    cex.main = .8, cex.lab = .8, cex.axis = .8, 
    xaxt = "n", ci.arg = list(col = adjustcolor(pal[1], .2)))
  
  # Add prediction
  lines(newERF[[i]],col = pal[2], lwd = 2, ci = "area",
    ci.arg = list(col = adjustcolor(pal[2], .2)))
  
  # MMT
  abline(v = c(firststageERF[[i]]$cen, newERF[[i]]$cen), lty = 2, col = pal)
  
  # Add axis
  axis(1, at = ovaxis, labels = axisper)
  abline(v = ovaxis, lty = 3, col = "lightgray")
  
  # Horizontal grid
  abline(h = axTicks(2), lty = 3, col = "lightgray")
  abline(h = 1)
  
  # Add legend
  if (i %% ncols == 1) legend("top", legend = c("First-stage", "Predicted"), 
    lwd = 1, col = pal, horiz = T, cex = .5, bg = "white")
}

dev.off()