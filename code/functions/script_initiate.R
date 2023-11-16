# script aim: This script loads required packages as well as all functions used in this study 
# load packages for analysis
source(paste0(packages.folder,'packages_to_load.R'))

# load functions used in many scripts
source(paste0(functions.folder,'functions.R'))

# you might face issue on installing libraries like PCPhelpers and pcpr as well as "esri2sf". Please use following code to install them from github repositories

# library(devtools)
# install_github("Columbia-PRIME/PCPhelpers")
#devtools::install_github("https://github.com/Columbia-PRIME/pcpr/tree/dev")
#install_github("yonghah/esri2sf") 
