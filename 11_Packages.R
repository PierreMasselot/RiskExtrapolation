################################################################################
#
#   Risk Extrapolation project
#   
#   Packages
#
################################################################################

#------------------------
# Packages
#------------------------

#----- Data management
library(eurostat) # To download Eurostat datasets
library(giscoR) # For geographic information in Europe
library(sf) # GIS functions
library(doParallel) # To parallelize
library(doSNOW)
library(data.table) # Efficient tables and function `between`
library(dplyr) # data.frame management
library(terra) # For rasters
library(PHEindicatormethods) # Includes standard European population
library(Matrix) # For nearPD function
library(reshape2) # For long-to-wide
library(tidyr)
library(zen4R) # To download data from Zenodo
library(R.utils) # To read gz files

#----- Analysis
library(dlnm) # Functions to perform DLNM on first stage
library(splines) # For access to B-splines
library(mixmeta) # Second-stage meta-analysis
library(pls) # To create indices of vulnerability
library(gstat) # Kriging / IDW
library(MASS) # For multivariate normal in simulations

#----- Vizualization
library(viridis)
library(corrplot)
library(ggplot2)
library(scales)
library(scico)
library(patchwork)
library(ggdist)
library(fields)
library(tidyterra)
library(ggnewscale) # To have several fill scales on the same plot
