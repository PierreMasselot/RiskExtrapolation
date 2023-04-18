################################################################################
#
#  Appendix 1.
#  Alternative ways to construct composite indices of vulnerability
#
################################################################################

# Needs the metapred and metapred_obs object created
# Replace the "Compute PLS" part of the code in 23 by the code below

#---------------------
# PCA
#---------------------
# Unsupervised dimension reduction useful in any context

# Compute PCA
pcares <- prcomp(metapred_obs, center = TRUE, scale = TRUE)

# Extract PCs
comps <- scale(metapred, center = pcares$center, scale = pcares$scale) %*% 
  pcares$rotation

# Correlation with original metapredictors
coords <- t(apply(pcares$rotation, 1, "*", pcares$sdev))

#---------------------
# CCA
#---------------------
# Mostly useful when the response (ERF) is multidimensional

# Compute CCA
ccares <- cancor(metapred_obs, coefs_obs)

# Create new metavariables
comps <- scale(metapred, center = ccares$xcenter, scale = F) %*% ccares$xcoef
colnames(comps) <- sprintf("Comp%i", 1:ncol(metapred_obs))

# Correlation with original metapredictors
coords <- cor(metapred, comps)
