################################################################################
#
#   New extensions
#   1. Age
#
################################################################################

#######################################
# ESTIMATE SECOND_STAGE MODEL WITH AGE
#######################################

#-------------------------
# Prepare age data
#-------------------------

#----- Add average age of death for each age group of each city

# Order death rates by age
deathrates <- metadata[,sort(grep("deathrate_[[:digit:]]{4}", 
  names(metadata), value = T))]

# Detail of age groups
amin <- as.numeric(substr(names(deathrates), 11, 12))
amax <- as.numeric(substr(names(deathrates), 13, 14))
amean <- (amin + amax) / 2
arange <- amax - amin

# Compute average age of death inside each group
# Mean of ages weighted by their specific death rates
groups <- cut(amax, agebreaks)
deathages <- tapply(seq_along(deathrates), groups, function(ind){
  data.matrix(deathrates[,ind]) %*% (amean[ind] * arange[ind]) / 
    (data.matrix(deathrates[,ind]) %*% arange[ind])
})

# Add life expectancy for oldest age group
deathages[[tail(agelabs, 1)]] <- metadata[[sprintf("lifexp_%2.0f", 
  tail(agebreaks, 1))]] + tail(agebreaks, 1)

# Merge with stage 2 data.frame
deathagesdf <- data.frame(agegroup = rep(agelabs, each = nrow(metadata)),
  city = rep(metadata$CITY_CODE, length(agelabs)), 
  age = unlist(deathages))
stage2df <- merge(stage2df, deathagesdf, all.x = T)


#-------------------------
# Fit model
#-------------------------

#----- Compare various models of age

# List of compared models
compared <- c(
  "agegroup", # Age group as a factor
  "age", # Age as linear
  sprintf("ns(age, df = %i, Bound = c(0, 100))", 2:4), # Splines with various df
  sprintf("ns(age, knots = %i, Bound = c(0, 100))", 
    seq(60, 80, by = 5))) # Splines with predetermined knot locations

# Create formulas
formlist <- lapply(compared, function(f) update(fixform, 
  as.formula(sprintf("~ . + %s", f))))

# Apply mixmeta (can take several minutes)
reslist <- lapply(formlist, mixmeta, random = ranform, 
  data = stage2df, S = vcovs)

# Extract AIC
aicvec <- sapply(reslist, AIC)

# Select the best model
stage2res <- reslist[[which.min(aicvec)]]

#----- Prepare results

# Prepare prediction: at each break defined earlier (can be anything)
newdata <- data.frame(age = agebreaks[-1])

# Labels
reslabs <- as.character(agebreaks[-1])

#######################################
# PLOTS
#######################################

#------------------------
# AIC comparison
#------------------------

# Create data.frame
aicdf <- data.frame(aic = aicvec, 
  model = c("Factor", "Linear", sprintf("%idf", 2:4), 
    sprintf("Knot: %i", seq(60, 80, by = 5))),
  group = c("Age group", "Continuous: linear", 
    rep("Continuous: spline df", 3), 
    rep("Continuous: spline knots", 5)))

# Transform as factors
aicdf <- mutate(aicdf,
  model = factor(model, levels = model),
  group = factor(group, levels = unique(group))
)

# Choose color palette
pal <- mako(4, end = .8)
aicmin <- aicdf$aic == min(aicdf$aic)

# Plot
ggplot(aicdf) + theme_classic() + 
  geom_point(aes(x = aic, y = model, col = aicmin, fill = group, stroke = 2), 
    size = 3, shape = 21) + 
  scale_fill_manual(values = pal, guide = guide_legend(nrow = 2)) + 
  scale_colour_manual(values = c("transparent", 2), guide = "none") + 
  labs(x = "AIC", y = "", fill = "") + 
  ggtitle("a) AIC comparison") +
  theme(panel.grid.major.y = element_line(linetype = 2, 
      colour = ifelse(aicmin, 2, "grey")),
    legend.position = "bottom",
    axis.text.y = element_text(colour = ifelse(aicmin, 2, 1)))

# Export
ggsave("figures/Fig1a_AICage.pdf", height = 7, width = 5)

#------------------------
# Curves for various age (groups)
#------------------------

#----- Construct age-group specific ERF

# Predict coefficients for each age (group)
agepreds <- predict(stage2res, newdata, vcov = T)

# Determine common MMT
uncentred <- mmtbasis %*% sapply(agepreds, "[[", "fit")
mmts <- mmtper[apply(uncentred, 2, which.min)]
mmt <- median(mmts)

# Construct centred ERFs
ageERF <- lapply(agepreds, function(x){
  crosspred(ovbasis, coef = x$fit, vcov = x$vcov, model.link = "log", 
    at = ovper, cen = mmt)
})

#----- Plot age-group specific curves

# Color palette
agepal <- mako(length(ageERF), end = .8, direction = -1)

# Initialize plot, draw grid and cusotm x-axis
plot(NA, bty = "l", xaxt = "n", 
  xlab = "Temperature percentile", ylab = "RR",
  xlim = range(ovper), ylim = c(.8, 2))
abline(v = ovaxis, h = axTicks(2), lty = 2, col = "lightgrey")
axis(1, at = ovaxis, labels = axisper)
mtext("b) ERF at various ages", at = par("usr")[1], adj = 0, line = 1, 
  cex = 1.2)

# Add age curves
for (i in seq_along(ageERF)){
  lines(ageERF[[i]], ptype = "overall", col = agepal[i], ci = "area", 
    lwd = 2, ci.arg = list(col = adjustcolor(agepal[i], .2)))
}

# Overlay the RR=1 line
abline(h = 1)

# Add legend
legpars <- list(legend = reslabs, col = agepal, cex = .8,
  lty = 1, lwd = 2, title = "Age", horiz = T, xpd = T)
legdim <- do.call(legend, c(legpars, list(x = "center", plot = F)))
do.call(legend, c(legpars, 
  list(x = mean(par("usr")[1:2]) - legdim$rect$w / 2, 
    y = par("usr")[4])))

# Save
dev.print(pdf, file = "figures/Fig1b_age.pdf", width = 7, height = 5)
