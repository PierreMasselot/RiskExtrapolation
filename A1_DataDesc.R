################################################################################
#
#  Appendix 1.
#  Data description
#
################################################################################

#-----------------------
# List of cities
#-----------------------

# Select variables, sort, transform and rename
city_desc <- metadf |>
  subset(select = c(city_code, city_name, geozone, obs, pop, tmean)) |>
  arrange(city_code) |>
  mutate(obs = factor(obs, levels = c(T, F), labels = c("Yes", "No")),
    pop = formatC(pop, format = "f", big.mark = ",", digits = 0),
    tmean = formatC(tmean, format = "f", digits = 1)) |>
  rename(`Eurostat code` = "city_code", Name = "city_name", Region = "geozone", 
    Observed = "obs", Population = "pop", `Mean temperature` = "tmean")

# Export
write.table(city_desc, file = "figures/TableS1_citydesc.txt", sep = ";", 
  quote = F, row.names = F)
