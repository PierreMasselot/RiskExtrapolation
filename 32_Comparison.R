################################################################################
#
#  RiskExtrapolation
#
#  Comparison
#
################################################################################

#-------------------------
# ERF preparation
#-------------------------

#----- Compute ERFs

# Merge all coefs together
coefdf <- rbind(
  select(metafull, city_code, agegroup, starts_with("coef"), 
    starts_with("vcov")) |> mutate(estimate = "fs"),
  select(cityageres, city_code, agegroup, starts_with("coef"), 
    starts_with("vcov")) |> mutate(estimate = "pred"),
  mutate(validblups, estimate = "blup")
  # mutate(fixedres, estimate = "fixed"),
  # mutate(ageonlyres, estimate = "ageonly"),
  # mutate(compres, estimate = "componly"),
  # mutate(nullres, estimate = "null")
)

# Loop across cities, age groups and estimates
allerf <- foreach(ires = iter(coefdf, "row")) %do% 
{
  
  # Extract coefs
  coefs <- select(ires, starts_with("coef")) |> unlist()
  vcovs <- select(ires, starts_with("vcov")) |> unlist() |> xpndMat()

  # Compute ERF
  uncentred <- mmtbasis %*% coefs
  mmt <- mmtper[which.min(uncentred)]
  crosspred(ovbasis, coef = coefs, vcov = vcovs, 
    model.link = "log", at = ovper, cen = mmt)  
}

# Extract RRs with confidence intervals
erfdf <- foreach(ires = iter(coefdf, "row"), erf = allerf, 
  .final = rbindlist) %do% 
  {
    cbind(ires[,c("city_code", "agegroup", "estimate")], per = predper,
      erf[c("predvar", "allRRfit", "allRRlow", "allRRhigh", "allse", "cen")]) |>
      suppressWarnings()
  }

# Select unobserved cities add names, and create factors for plotting
nms <- unique(metadf[, c("city_code", "city_name")])
erfdf <- subset(erfdf, !city_code %in% obs) |>
  mutate(erf = factor(estimate, 
    levels = c("fs", "null", "ageonly", "componly", "fixed", "pred", "blup"),
    labels = c("First-stage", "Null", "Age only", "Components only", 
      "Fixed only", "Predicted", "Reference")),
    city = factor(city_code, levels = nms$city_code, labels = nms$city_name))
  
#-------------------------
# Plots
#-------------------------

#----- Graphical parameters

# Age-group labels
ageplot <- agelabs
ageplot[-length(ageplot)] <- gsub(pattern = "(.{2})(.*)",
  replacement = "\\1-\\2", ageplot[-length(ageplot)])
names(ageplot) <- agelabs

# Model labels
estlabs <- c("First-stage", "Null", "Age only", "Components only", 
  "Fixed only", "Full prediction")
names(estlabs) <- c("fs", "null", "ageonly", "componly", "fixed", "pred")

# Palette for estimates
erfpal <- mako(length(unique(erfdf$erf)), direction = -1)
names(erfpal) <- unique(erfdf$erf)

#----- Plot all exposure-response functions

# plot by city and age
erfplot <- ggplot(subset(erfdf, estimate %in% c("fs", "pred", "blup"))) + 
  theme_minimal() + 
  facet_grid_paginate(city ~ agegroup, nrow = 7, ncol = 5, 
    labeller = labeller(agegroup = ageplot)) + 
  geom_ribbon(aes(x = predvar, ymin = allRRlow, ymax = allRRhigh, 
    fill = erf), alpha = .2) + 
  geom_line(aes(x = predvar, y = allRRfit, col = erf)) + 
  coord_cartesian(ylim = c(.8, 2)) +
  scale_x_continuous(name = "Temperature percentile", breaks = ovaxis, 
    labels = axisper) +
  ylab("RR") + 
  geom_hline(yintercept = 1, linewidth = .3) + 
  # scale_color_manual(values = erfpal, name = "Exposure-response\nfunction") +
  # scale_fill_manual(values = erfpal, name = "Exposure-response\nfunction") +
  scale_color_viridis_d(option = "G", direction = -1, end = .8,
    name = "Exposure-response\nfunction") +
  scale_fill_viridis_d(option = "G", direction = -1, alpha = .2, end = .8,
    name = "Exposure-response\nfunction") +
  theme(panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linetype = 2, color = "grey", 
      linewidth = .1),
    panel.spacing = unit(1, "lines"),
    axis.text = element_text(size = 6))

# Export
for (i in 1:n_pages(erfplot)) {
  pi <- erfplot + 
    facet_grid_paginate(city ~ agegroup, nrow = 7, ncol = 5, page = i, 
      labeller = labeller(agegroup = ageplot))
  ggsave(sprintf("figures/SFig_ERFcomparison_%i.pdf", i), pi, height = 8)
}

#----- Plot specific features of ERFs
# 
# # Extract 1st and 99th percentiles and remove extreme ones
# erfsum <- subset(erfdf, per %in% c(1, 99)) |> 
#   subset(allRRfit < 3) |>
#   mutate(res = allRRfit)
# 
# # Extract MMT
# mmtdf <- subset(erfdf, predvar == cen) |> mutate(res = cen, per = "mmt")
# erfsum <- rbind(erfsum, mmtdf) |>
#   mutate(per = factor(per, c("1", "99", "mmt"), 
#     labels = c("Cold", "Heat", "MMT")))
# 
# # Palettes
# coldpal <- mako(2, begin = .3, end = .7, direction = -1)
# heatpal <- rocket(2, begin = .3, end = .7, direction = -1)
# mmtpal <- grey.colors(2)

# # Plot
# ggplot(erfsum) + theme_minimal() + 
#   geom_point(aes(x = city, y = res, col = per, shape = erf, alpha = erf)) + 
#   geom_line(aes(x = city, y = res, col = per, group = erf, alpha = erf)) + 
#   # facet_wrap(~ agegroup + per, scales = "free_y", ncol = 3) + 
#   ggh4x::facet_grid2(rows = vars(agegroup), cols = vars(per), scales = "free",
#     independent = "y") + 
#   scale_color_manual(values = c("Cold" = 4, "Heat" = 2, "MMT" = 1), 
#     guide = "none") + 
#   scale_alpha_manual(name = "Estimate", 
#     values = c("First-stage" = .4, "BLUP" = .7, "Extrapolated" = 1)) + 
#   scale_shape(name = "Estimate") + 
#   labs(x = "", y = "") + 
#   theme(panel.grid.minor = element_blank(), 
#     panel.grid.major.x = element_blank(),
#     axis.text.x.bottom = element_text(angle = 90, hjust = 1, size = 5),
#     strip.text = element_text(size = 12),
#     panel.spacing = unit(1, "lines"))


#------ Alternative

# # Pivot
# erfsum <- pivot_wider(erfsum,
#   id_cols = all_of(c("city_code", "agegroup", "per")),
#   names_from = "erf", values_from = "res")
# 
# # Compute regression results
# lmdata <- summarise(erfsum, .by = c("per", "agegroup"), 
#     r2 = summary(lm(Extrapolated ~ BLUP))$r.squared) |>
#   mutate(r2 = sprintf("%2.0f%%", r2 * 100))
# 
# # Plot
# ggplot(erfsum) + theme_minimal() + 
#   ggh4x::facet_grid2(rows = vars(agegroup), cols = vars(per), scales = "free",
#     independent = "y", labeller = labeller(agegroup = ageplot)) + 
#   geom_point(aes(x = BLUP, y = Extrapolated, col = per)) + 
#   geom_smooth(aes(x = BLUP, y = Extrapolated, col = per, fill = per),
#     method = "lm", alpha = .2) + 
#   scale_color_manual(values = c("Cold" = 4, "Heat" = 2, "MMT" = 1), 
#     guide = "none") + 
#   scale_fill_manual(values = c("Cold" = 4, "Heat" = 2, "MMT" = 1), 
#     guide = "none") + 
#   geom_text(aes(x = Inf, y = -Inf, vjust = -.3, hjust = 1.05, label = r2), 
#     data = lmdata, fontface = "bold") + 
#   theme(panel.grid.minor = element_blank(),
#     strip.text = element_text(size = 12),
#     panel.spacing = unit(1, "lines"),
#     panel.border = element_rect(fill = NA))
#   
# 
# ggsave("figures/SFig_resComparison.pdf", width = 10)

#----- Plot RMSE

# Compute squared error
sqerrdf <- pivot_wider(erfdf, id_cols = c("city_code", "agegroup", "per"),
    names_from = estimate, values_from = c("allRRfit", "allse")) |>
  mutate(across(starts_with("allRRfit"), ~ (log(.x) - log(allRRfit_blup))^2)) |>
  pivot_longer(starts_with("allRRfit") | starts_with("allse"), 
    names_to = c(".value", "estimate"), names_sep = "_") |>
  subset(estimate != "blup") |>
  rename(sqerr = "allRRfit")

# Overall RMSE
rmsedf <- mutate(sqerrdf, w = 1/(allse + 1e-8)^2) |>
  summarise(total = sqrt(weighted.mean(sqerr, w)), 
    cold = sqrt(weighted.mean(sqerr[per <= 5], w[per <= 5])),
    heat = sqrt(weighted.mean(sqerr[per >= 95], w[per >= 95])),
    .by = c("agegroup", "estimate"))

# Refactor
rmsedf <- pivot_longer(rmsedf, c("total", "cold", "heat"),
    names_to = "range", values_to = "rmse")

# Plot
ggplot(rmsedf) + theme_minimal() + 
  facet_grid(rows = vars(agegroup), cols = vars(range),
    labeller = labeller(agegroup = ageplot, 
      range = c(cold = "Cold", heat = "Heat", total = "Total"))) + 
  geom_col(aes(x = estimate, y = 100 * (exp(rmse) - 1), fill = range, 
    alpha = estimate)) + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = c("cold" = 4, "heat" = 2, "total" = 1), 
    guide = "none") + 
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  scale_x_discrete(labels = estlabs) +
  labs(x = "", y = "RMSE (RR ratio %)") + 
  theme(axis.text.x.bottom = element_text(angle = -45, vjust = .5, hjust = 0),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = .1),
    axis.text.y = element_text(size = 6))

# Save
ggsave("figures/Fig9_RMSE.pdf", width = 7)

# #----- Plot RMSE curves
# 
# # Compute RMSE for each percentile
# rmsecurves <- summarise(sqerrdf, rmse = sqrt(mean(sqerr)), 
#     .by = c("estimate", "per")) |>
#   mutate(estimate = factor(estimate, 
#     levels = c("fs", "null", "ageonly", "componly", "fixed", "pred"),
#     labels = c("First-stage", "Null", "Age only", "Components only", 
#       "Fixed only", "Full prediction")))
# 
# # Plot
# ggplot(subset(rmsecurves, estimate != "First-stage")) + theme_minimal() +
#   geom_line(aes(x = per, y = rmse, col = estimate, linetype = estimate)) + 
#   scale_color_viridis_d(option = "G", direction = -1, end = .8)
  