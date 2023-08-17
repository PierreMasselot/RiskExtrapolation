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
agebreaks <- substr(agelabs, 1, 2) |> as.numeric()
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
stage2_obs <- merge(stage2_obs, deathagesdf, all.x = T)
stage2_pred <- merge(stage2_pred, deathagesdf, all.x = T)
