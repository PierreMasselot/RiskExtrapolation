################################################################################
#
#  Appendix 2.
#  Meta regression model selection
#
################################################################################

# Should be run after script 23

#----- Parameters

# Maximum number of components to try
maxk <- 10

# List of models for age
compared <- c(
  "agegroup", # Age group as a factor
  "age", # Age as linear
  sprintf("ns(age, knots = %i, Bound = c(0, 100))", seq(60, 80, by = 5)), 
    # Splines with predetermined knot locations
  sprintf("ns(age, df = %i, Bound = c(0, 100))", 3:4) # Splines with various df
) 

#-------------------------
# Grid search of best model
#-------------------------

# Create full factorial grid
modgrid <- expand.grid(ncomp = 0:maxk, agemod = compared)

# Initialise random formula
ranform <- ~ 1|city_code

# Prepare parallelisation
ncores <- detectCores()
cl <- makeCluster(max(1, ncores - 2))
registerDoParallel(cl)

# Loop across grid
aiclist <- foreach(k = modgrid$ncomp, a = modgrid$agemod, .combine = c, 
  .packages = c("mixmeta", "splines")) %dopar% 
{
  
  # Create formula
  if (k == 0){
    kaform <- sprintf("cbind(%s) ~ %s", 
      paste(coefvars, collapse = ", "), a)
  } else {
    kaform <- sprintf("cbind(%s) ~ %s + %s", 
      paste(coefvars, collapse = ", "), a,
      paste(colnames(comps)[1:k], collapse = " + "))
  }
  
  # Fit model
  fit <- mixmeta(as.formula(kaform), random = ranform, 
    data = metadf, S = smat, subset = conv) 
  
  # Extract AIC
  AIC(fit)
}

stopCluster(cl)

#-------------------------
# Plot results
#-------------------------

# Model labels
names(compared) <- c("Factor", "Linear", 
  sprintf("NS (knot = %i)", seq(60, 80, by = 5)),
  sprintf("NS (%i DF)", 3:4))

# Color palette
pal <- mako(length(compared), end = 1)
shppal <- rep_len(rep(15:18), length(compared))

# Plot data.frame
plotdf <- cbind(modgrid, deltaAIC = aiclist - min(aiclist))

# Plot
ggplot(plotdf, aes(x = ncomp, y = deltaAIC, group = agemod, col = agemod, 
    shape = agemod)) + 
  theme_classic() + 
  geom_line(linewidth = 1) + geom_point(size = 4) + 
  scale_color_manual(values = pal, labels = names(compared), 
    name = "Age model") + 
  scale_shape_manual(values = shppal, name = "Age model", 
    labels = names(compared)) +
  scale_x_continuous(name = "Number of components", breaks = 1:maxk) + 
  scale_y_continuous(name = "\U0394 AIC") + 
  geom_hline(yintercept = c(0, 2, 5, 10), linetype = c(1, 2, 2, 2))

# Save
dev.print(cairo_pdf, "figures/FigS1_npcSelect.pdf")

