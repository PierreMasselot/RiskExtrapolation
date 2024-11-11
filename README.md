# A new modelling framework for multi-location studies in environmental epidemiology
 
Fully reproducible code and data for the illustrative application in the paper *A new modelling framework for multi-location studies in environmental epidemiology*. 

All components of the dataset are provided in the *data* subdirectory. The subdirectory also contains the script `00_DownloadData.RData` that shows how the data were obtained and pre-processed. The only exception is *data/tmean.csv.gz*, that was extracted using Python from ERA5Land. :warning: This script should not be executed, it is here mainly to document the data gathering process.

Execution of the code should start from *11_Packages.R*. 

1. Scripts *11_Packages.R* to *13_PrepData.R* setup the modelling and load data.

2. Scripts *21_LocationAgeSpecific.R* to *28_Uncertainty.R* implement each corresponding package (sub-section) in section 3 of the paper. Each one rely on previous ones, which means they should be ran in order.

3. Scripts *A1_DataDesc.R* to *A4_AlternativePC.R* provides additional results corresponding to the supplementary materials.

Finally, the *figures* subdirectory includes all Figures and Tables of the paper.
