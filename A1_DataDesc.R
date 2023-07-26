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
city_desc <- metadata |>
  subset(select = c(CITY_CODE, LABEL, geozone, obs, pop, tmean)) |>
  arrange(CITY_CODE) |>
  mutate(obs = factor(obs, levels = c(T, F), labels = c("Yes", "No")),
    pop = formatC(pop, format = "f", big.mark = ",", digits = 0),
    tmean = formatC(tmean, format = "f", digits = 1)) |>
  rename(`Eurostat code` = "CITY_CODE", Name = "LABEL", Region = "geozone", 
    Observed = "obs", Population = "pop", `Mean temperature` = "tmean")

# Export
write.table(city_desc, file = "figures/TableS1_citydesc.txt", sep = ";", 
  quote = F, row.names = F)
