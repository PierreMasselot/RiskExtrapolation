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

# Average population proportion in Italy
isp <- aggregate(popprop ~ agegroup, cityageres, mean)$popprop

#-----------------------------
# Compute standardized excess death rates
#-----------------------------

#----- Compute attributable fractions

# Loop on city ages
aflist <- foreach(pred = predcoefs, dat = dlist[cityageres$city], 
  .combine = rbind) %do% 
{
  
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
  
  # Compute daily AF
  afday <- (1 - exp(-bvarcen %*% pred$fit))
  
  # Indicator of heat days
  heatind <- dat$tmean >= immt
  
  # Sum all and return, including the MMT
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

# Compute standardized excess rates
stdexcess_region <- aggregate(excessrate ~ geozone, regionageres, 
  weighted.mean, w = isp)
colnames(stdexcess_region)[-1] <- c("total", "cold", "heat")

#----- Compute non-standardized excess death rates for comparison

# Aggregate by city
totexcessrate <- cityageres |> group_by(city) |>
  # Denominator is total population
  summarise(totexcessrate = t(byrate * colSums(an) / sum(agepop))) |>
  do.call(what = data.frame) |>
  rename(total = "totexcessrate.1", cold = "totexcessrate.2", 
    heat = "totexcessrate.3")

####################################
# Plots
####################################

# Select cities to display based in demographic structure
citysel <- cityageres |> subset(agegroup == "85+") |>
  subset(rank(popprop, ties.method = "first") %in% 
      round(seq(1, n, length.out = 5)) & !duplicated(popprop), 
    city, drop = T) 
citysel <- as.character(citysel[order(subset(cityageres, 
  city %in% citysel & agegroup == "85+", popprop, drop = T))])

#----- Compare age distribution

# Prepare df
demoplot <- subset(demodf, CITY_CODE %in% citysel) |> 
  mutate(CITY_CODE = factor(CITY_CODE, levels = unique(CITY_CODE), 
      labels = unique(CITY_NAME)), 
    city = reorder(CITY_CODE, prop, tail, 1))

# Plot
figa <- ggplot(demoplot) + theme_classic() + 
  geom_smooth(aes(x = agemin, y = prop / 5, col = city), 
    linewidth = 2, se = F) + 
  scale_colour_viridis(option = "mako", end = .9, name = "", discrete = T) + 
  labs(x = "Age", y = "Population percentage (%)", title = "A") + 
  theme(legend.position = "top")

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
pal <- mako(2, begin = .2, end = .7)
names(pal) <- c("tot", "std")

# Plot
figb <- ggplot(excessplot) + theme_classic() + 
  geom_col(aes(x = city, y = heat, fill = type), 
    position = "dodge", width = .7) + 
  scale_x_discrete(breaks = citysel,
    labels = with(metadata, CITY_NAME[match(citysel, CITY_CODE)])) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = pal, name = "", 
    labels = c(std = "Standardised", tot = "Non-standardised")) +
  labs(y = "Excess motality rate (x 100,000)", x = "", title = "B") + 
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = 2),
    legend.position = "top") + 
  coord_flip()

#----- Put Figures together

# Wrap plots
wrap_plots(figa, figb, nrow = 1)

# Save
ggsave("figures/Fig6_AgeStandardisation.pdf", width = 11)
