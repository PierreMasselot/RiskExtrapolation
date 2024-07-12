################################################################################
#
#  RiskExtrapolation
#
#  Standardisation of impacts
#
################################################################################

#-----------------------------
# Initialize age groups informations
#-----------------------------

#----- Prepare demographic data

# Add info about population 
demodat <- mutate(metadata_age, agegroup = cut(age, 
    breaks = c(0, agebreaks, Inf), labels = agelabs, right = F)) |>
  merge(subset(metadf, select = c(city_code, agegroup, pop)), all.x = T)

# Compute age-specific demographic data
demodat <- mutate(demodat, pop = prop * pop / 100, deaths = pop * deathrate)

# Aggregate by age group to create results data.frame
demogrps <- summarise(demodat, 
  popprop = sum(prop), pop = sum(pop), death = sum(deaths),
  .by = c("city_code", "agegroup"))
cityageres <- merge(cityageres, demogrps, all.x = T)
    
#----- Prepare standard population

# Average population proportion in Italy
isp <- aggregate(popprop ~ agegroup, cityageres, mean)
names(isp)[2] <- "ispw"

#-----------------------------
# Compute standardized excess death rates
#-----------------------------

#----- Compute city-age level results

# Loop on city ages to compute attributable bumbers
cityageres <- foreach(ires = iter(cityageres, by = "row"), 
  .combine = rbind) %do% 
{
  
  # Extract tmean series
  tmean <- tempdata[city_code == ires$city_code,]$tmean
  
  # Extract coefficients
  coefs <- ires[,coefvars] |> unlist()
  
  # Estimate MMT
  bvar <- onebasis(tmean, fun = varfun, degree = vardegree, 
    knots = quantile(tmean, varper / 100))
  erf <- bvar %*% coefs
  inrange <- between(tmean, quantile(tmean, mmprange[1] / 100), 
    quantile(tmean, mmprange[2] / 100))
  immt <- tmean[inrange][which.min(erf[inrange])]
  
  # Centre basis
  cenvec <- onebasis(immt, fun = varfun, degree = vardegree, 
    knots = quantile(tmean, varper / 100), Bound = range(tmean))
  bvarcen <- scale(bvar, center = cenvec, scale = F)
  
  # Compute daily AF
  anday <- (1 - exp(-bvarcen %*% coefs)) * ires$death / length(tmean)
  
  # Indicator of heat days
  heatind <- tmean >= immt
  
  # Sum all and return, including the MMT
  cbind(ires, an_total = sum(anday), an_cold = sum(anday[!heatind]), 
    an_heat = sum(anday[heatind]), mmt = immt)
}

# Compute excess death rates
cityageres <- mutate(cityageres, across(starts_with("an_"), 
    ~ byrate * .x / pop, .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

#----- Compute city-level results

# Add italian standard population as weights
cityres <- merge(cityageres, isp) |>
  # Aggregate impacts
  summarise(
    # Attributable numbers: sum age specific ANs
    across(starts_with("an_"), sum), pop = sum(pop),
    # Compute standardised excess death rates
    across(starts_with("excessrate"), ~ weighted.mean(.x, w = ispw)),
    .by = city_code) |>
  rename_with(~ gsub("excess", "std", .x))

# Compute non-standardised excess rates
cityres <- mutate(cityres, across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

#----- Compute regional level excess death rates

# Aggregate by region, starting by summing AN and pop
regionageres <- summarise(cityageres, 
    across(starts_with("an_"), sum), pop = sum(pop),
    .by = c(geozone, agegroup)) |>
  # Compute excess death rates
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))

# Aggregate across age levels, starting by merging standard pop weights
regionres <-  merge(regionageres, isp) |>
  # Compute total AN and std death rates
  summarise(across(starts_with("an_"), sum), pop = sum(pop),
    across(starts_with("excessrate"), weighted.mean, w = ispw),
    .by = geozone) |>
  rename_with(~ gsub("excess", "std", .x)) |>
  # Non-standardised rates for comparison
  mutate(across(starts_with("an_"), ~ byrate * .x / pop, 
    .names = "excessrate_{.col}")) |>
  rename_with(~ gsub("_an_", "_", .x))



####################################
# Plots
####################################

# Select cities to display based in demographic structure
citysel <- cityageres |> subset(agegroup == "85+") |>
  subset(rank(popprop, ties.method = "first") %in% 
      round(seq(1, n, length.out = 5)) & !duplicated(popprop), 
    city_code, drop = T) 
citysel <- as.character(citysel[order(subset(cityageres, 
  city_code %in% citysel & agegroup == "85+", popprop, drop = T))])
namesel <- with(metadf, city_name[match(citysel, city_code)])

#----- Compare age distribution

# Prepare df
demoplot <- subset(demodat, city_code %in% citysel) |> 
  merge(unique(metadf[, c("city_code", "city_name")])) |>
  mutate(city_name = factor(city_code, levels = citysel, labels = namesel))

# Plot
figa <- ggplot(demoplot) + theme_classic() + 
  geom_smooth(aes(x = age, y = prop / 5, col = city_name), 
    linewidth = 2, se = F) + 
  scale_colour_viridis(option = "mako", end = .9, name = "", discrete = T) + 
  labs(x = "Age", y = "Population percentage (%)", title = "A") + 
  theme(legend.position = "top")

#----- Compare standardised and non-standardised

# Select cities, pivot for ggplot and relabel
excessdf <- subset(cityres, city_code %in% citysel, 
    select = c(city_code, stdrate_heat, excessrate_heat)) |>
  pivot_longer(cols = c(stdrate_heat, excessrate_heat), names_to = "type",
    values_to = "excessrate") |>
  mutate(type = factor(type, c("stdrate_heat", "excessrate_heat"), 
      labels = c("Standardised", "Non-standardised")),
    city_name = factor(city_code, levels = rev(citysel), labels = rev(namesel)))
  
# Color palette
pal <- mako(2, begin = .2, end = .7)
names(pal) <- levels(excessdf$type) |> rev()

# Plot
figb <- ggplot(excessdf) + theme_classic() + 
  geom_col(aes(x = city_name, y = excessrate, fill = type), 
    position = "dodge", width = .7) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = pal, name = "") +
  labs(y = "Excess motality rate (x 100,000)", x = "", title = "B") + 
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = 2),
    legend.position = "top") + 
  coord_flip()

#----- Put Figures together

# Wrap plots
wrap_plots(figa, figb, nrow = 1)

# Save
ggsave("figures/Fig7_AgeStandardisation.pdf", width = 11)
