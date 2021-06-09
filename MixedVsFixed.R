##########################################################################
#
#                     RISK EXTRAPOLATION
#                   Exploration Mixed vs fixed
#
##########################################################################

library(mixmeta)
library(pls)
library(rsample)
library(ggplot2)

load("data/MCC-EUcityTRM.RData")

#----------------------------
# Parameters
#----------------------------

# Number of PLS components
ncomp <- 7

# Number of folds for CV
nfold <- 10

#----------------------------
# Prepare data
#----------------------------

# Coefs
coefs <- t(sapply(stage1res, "[[", "coef"))

# Vcovs
vcovs <- lapply(stage1res, "[[", "vcov")

# Prepare PLS variables
metapls <- scale(metapred)
plsres <- plsr(coefs ~ metapls)
plsvar <- predict(plsres, ncomp = 1:ncomp, type = "scores")
colnames(plsvar) <- sprintf("pls%i", 1:ncomp)

# Split data
splitinds <- vfold_cv(data.frame(1:nrow(coefs)), v = nfold)

#----------------------------
# CV on different levels of mixed effects
#----------------------------

rmselist <- list()

# Loop on folds
for (i in seq_len(nfold)){
  
  cat(i)
  
  # Get training and validation indices
  train <- unlist(analysis(splitinds$splits[[i]]))
  valid <- unlist(assessment(splitinds$splits[[i]]))
  
  #----- Basic linear model
  lmres <- lm(coefs ~ plsvar, subset = train)
  lmerr <- coefs[valid,] - 
    predict(lmres, newdata = model.frame(coefs ~ plsvar)[valid,])
  rmselist$lm[i] <- sqrt(mean(lmerr^2))
  
  #----- Fixed effects
  fixedres <- mixmeta(coefs ~ plsvar, S = vcovs, method = "fixed",
    subset = train)
  fixederr <- coefs[valid,] - 
    predict(fixedres, newdata = model.frame(coefs ~ plsvar)[valid,])
  rmselist$fixed[i] <- sqrt(mean(fixederr^2))
  
  #----- City level random effects
  
  cityranres <- mixmeta(coefs ~ plsvar, data = cities, S = vcovs, 
    random = ~ 1|city, subset = train)
  cityerr <- coefs[valid,] - 
    predict(cityranres, newdata = model.frame(coefs ~ plsvar)[valid,])
  rmselist$cityran[i] <- sqrt(mean(cityerr^2))
  
  #----- City and country level random effects
  
  countryranres <- mixmeta(coefs ~ plsvar, data = cities, S = vcovs, 
    random = ~ 1|country/city, control = list(maxiter = 10), subset = train)
  countryerr <- coefs[valid,] - 
    predict(countryranres, newdata = model.frame(coefs ~ plsvar)[valid,])
  rmselist$countryran[i] <- sqrt(mean(countryerr^2))
}

#----- Compute CV scores
cvm <- sapply(rmselist, mean)
cvsd <- sapply(rmselist, sd) / nfold

#----------------------------
# Plot
#----------------------------

ggplot(data.frame(mod = names(rmselist), cvm, cvsd), aes(x = mod)) + 
  theme_classic() + 
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd)) + 
  geom_point(aes(y = cvm), pch = 16, size = 5) + 
  ylab("Cross-validated RMSE") + xlab("Model") + 
  scale_x_discrete(labels = c("lm" = "Linear model", "fixed" = "Fixed",
    "cityran" = "Mixed city", "countryran" = "Mixed city & country"),
    limits = names(rmselist))

ggsave("figures/FixedvsRandom.pdf", device = "pdf")
