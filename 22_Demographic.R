################################################################################
#
#  RiskExtrapolation
#
#  Modelling demographic differences
#
################################################################################

#######################################
# ESTIMATE SECOND_STAGE MODEL WITH AGE
#######################################

#-------------------------
# Prepare age data
#-------------------------

#----- Add average age of death for each age group of each city

# Detailed info about age groups
agegrs <-  mutate(metadf, agelow = as.numeric(substring(agegroup, 1, 2)),
    agehigh = as.numeric(substring(agegroup, 3))) |>
  subset(select = c(city_code, agegroup, agelow, agehigh))

# Extract all ages represented in each age group
ages <- subset(agegrs, !is.na(agehigh)) |>
  reframe(age = agelow:agehigh, .by = c(city_code, agegroup))
agesdemo <- subset(metadata_age, !is.na(agehigh)) |>
  reframe(age = age:agehigh, deathrate,  .by = c(city_code, age))

# Compute weighted average age
ages <- merge(ages, agesdemo, all.x = T) |>
  summarise(age = weighted.mean(age, deathrate), .by = c(city_code, agegroup))

# Life expectancy for oldest age group
oldage <- subset(agegrs, is.na(agehigh)) |>
  merge(metadata_age, by.x = c("city_code", "agelow"), 
    by.y = c("city_code", "age")) |>
  mutate(age = agelow + lifexp)
ages <- rbind(ages, oldage[, c("city_code", "agegroup", "age")]) |>
  arrange(city_code, agegroup)

# Merge with metadata
metadf <- merge(metadf, ages)
