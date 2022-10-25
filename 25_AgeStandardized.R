################################################################################
#
#   New extensions
#   5. Standardized impacts
#
################################################################################

#-----------------------------
# Initialize age groups informations
#-----------------------------

#----- Prepare demographic data

# Transform as long
demodf <- pivot_longer(metadata, 
  cols = matches("(prop|deathrate)_[[:digit:]]{4}"),
  names_to = c(".value", "age"), names_pattern = "(.*)_(.*)")

# Create new variables
demodf <- mutate(demodf, 
  # About age group
  agemin = as.numeric(substr(age, 1, 2)),
  agemax = as.numeric(substr(age, 3, 4)),
  agelength = agemax - agemin,
  agegroup = cut(agemax, breaks = c(agebreaks, 100), labels = agelabs),
  # Demographics
  agepop = prop * pop / 100,
  deaths = agepop * deathrate
)

# Aggregate by age group
demogrp <- demodf |> group_by(CITY_CODE, agegroup) |>
  summarise(popprop = sum(prop), agepop = sum(agepop), death = sum(deaths))

# Create result data.frame
cityageres <- merge(stage2df, demogrp, by.x = c("city", "agegroup"),
  by.y = c("CITY_CODE", "agegroup"))
    
#----- Prepare standard population

# # Define minimum age of groups (5 years bands) see ?esp2013
# espbreaks <- (seq_along(esp2013) - 1) * 5
# 
# # Create groups
# espgrps <- cut(espbreaks, c(agebreaks, 100), right = F, labels = agelabs)
# 
# # Sum standard population for each group
# esptot <- tapply(esp2013, espgrps, sum)

# Average population proportion in Italy
isp <- aggregate(popprop ~ agegroup, cityageres, mean)$popprop

#-----------------------------
# Compute standardized excess death rates
#-----------------------------

#----- Compute attributable fractions

# Loop on city ages
aflist <- foreach(pred = predcoefs, dat = dlist[cityageres$city],
  .combine = rbind) %do% {
  
  # Estimate MMT
  bvar <- onebasis(dat$tmean, fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100))
  erf <- bvar %*% pred$fit
  inrange <- between(dat$tmean, quantile(dat$tmean, mmprange[1] / 100), 
    quantile(dat$tmean, mmprange[2] / 100))
  immt <- dat$tmean[inrange][which.min(erf[inrange])]
  
  # Centre basis
  cenvec <- onebasis(immt, fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100), 
    Boundary.knots = range(dat$tmean))
  bvarcen <- scale(bvar, center = cenvec, scale = F)
  
  # Compute daily AF and AN
  afday <- (1 - exp(-bvarcen %*% pred$fit))
  
  # Indicator of heat days
  heatind <- dat$tmean >= immt
  
  # Sum all and return
  c(total = sum(afday), cold = sum(afday[!heatind]), heat = sum(afday[heatind]))
}

#----- Compute standardized excess death rates

cityageres <- mutate(cityageres,
  # Compute annual AN
  an = aflist * death / nts,
  
  # Compute death rate
  excessrate = byrate * an / agepop,
)

# Compute standardized death rates for each city
stdexcessrate <- aggregate(excessrate ~ city, cityageres, 
  weighted.mean, w = isp)

#----- Compute regional level excess death rates

# Aggregate by region
regionageres <- cityageres |> group_by(agegroup, geozone) |>
  summarise(
    an = t(colSums(an)), 
    agepop = sum(agepop), 
    excessrate = byrate * an / agepop
  )

# Computae standradized excess rates
stdexcess_region <- aggregate(excessrate ~ geozone, regionageres, 
  weighted.mean, w = isp)
colnames(stdexcess_region)[-1] <- c("total", "cold", "heat")

#----- Compute non-standardized excess death rates for comparison

# Aggregate by city
totexcessrate <- cityageres |> group_by(city) |>
  
  # Denominator is total population
  summarise(totexcessrate = byrate * colSums(an) / sum(agepop), 
    type = colnames(an)) |>
  
  # To have the excess rates as columns
  dcast(city ~ type, value.var = "totexcessrate")

#-----------------------------
# Plots
#-----------------------------

# Select cities to display
citysel <- cityageres |> subset(agegroup == "85+")|>
  subset(popprop %in% quantile(popprop, seq(0, 1, length.out = 5)) & 
      !duplicated(popprop), city, drop = T) 
citysel <- as.character(citysel[order(subset(cityageres, 
  city %in% citysel & agegroup == "85+", popprop, drop = T))])

#----- Compare age distribution

# Prepare df
demoplot <- subset(demodf, CITY_CODE %in% citysel) |> 
  mutate(CITY_CODE = factor(CITY_CODE, levels = unique(CITY_CODE), 
      labels = unique(CITY_NAME)), 
    city = reorder(CITY_CODE, prop, tail, 1))

# Plot
ggplot(demoplot) + theme_classic() + 
  geom_step(aes(x = agemin, y = prop, col = city), size = 2) + 
  scale_colour_scico_d(palette = "oslo", end = .8, name = "") + 
  labs(x = "Age", y = "Population percentage (%)")

# Save
ggsave("figures/Fig5a_popStructure.pdf")

#----- Compare standardised and non-standardised

# Bind the data.frames
excessdf <- rbind(totexcessrate, stdexcessrate)
excessdf <- mutate(excessdf,
  type = rep(c("tot", "std"), each = n))

# Select
excessplot <- subset(excessdf, city %in% citysel)
excessplot <- mutate(excessplot, 
  city = factor(as.character(city), rev(as.character(citysel))))

# Color palette
pal <- scico(2, palette = "oslo", begin = .2, end = .7)
names(pal) <- c("tot", "std")

# Plot
ggplot(excessplot) + theme_classic() + 
  geom_col(aes(x = city, y = heat, fill = type, col = city), 
    position = "dodge", width = .7, size = 2) + 
  scale_x_discrete(breaks = citysel,
    labels = with(metadata, CITY_NAME[match(citysel, CITY_CODE)])) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = pal, name = "", 
    labels = c(std = "Standardised", tot = "Non-standardised")) +
  scale_colour_scico_d(palette = "oslo", end = .8, name = "", guide = "none") +
  labs(y = "Excess motality rate (x 100,000)", x = "") + 
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = 2)) + 
  coord_flip()

# Save
ggsave("figures/Fig5b_ExcessComparison.pdf", height = 8)
