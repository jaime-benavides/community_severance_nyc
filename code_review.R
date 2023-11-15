## Code Review 
## Development of a Community Severance Index for urban areas in the United States: A case study in New York City
## Vijay Kumar 
## Updated: Nov 15, 2023

# script: line:
# comment: Jaime, can we add chunk of code or details for installing the PCPhelpers and pcpr packages?
#comment:The way I installed instead of downloading them
library(devtools)
install_github("Columbia-PRIME/PCPhelpers") #worked fine 
install_github("Columbia-PRIME/pcpr") #having error, also tried Lawrence suggested way but have same error 
install.packages("pcpr-master", repos = NULL, type="source")
#installing *source* package ‘pcpr’ ...
# ** using staged installation
# Then try reinstalling 'pcpr'
remotes::install_github("Columbia-PRIME/pcpr")

# Jaime; In general I have two issues so far; 1- I coudn't install pcpr package, and 2- couldn't clone this repository to R; maybe you need to grant me persmission!
#fixed packages issues but we need to add details on packages like line number 33-35
#checking "a_01_preproc_smart_location_dta.R" packages issue "esri2sf", "rgeos" and "pcpr"
install_github("yonghah/esri2sf") #worked!
#other two packages still have problem but I could run a_01_preproc_smart_location_dta.R file and generate the data

# added on line 15 of code "#Selection of urban spatial variables, details in table 1 of manuscript"

#checking "a_02_barrier_factor_prep.R" no data folder and file roads/osm_driving_network_northeast.rds at line #36
#to fix this, data directory has been changed to generated data in script at line#

#
#data exploration script b_01_community_serverance_index; changed line 53 to correct data directory faf5_network <- sf::read_sf(paste0(geometry.data.folder, "FAF5Network.gdb"))
#working fine and understandable. 

#after talking to Jaime, fixed pcpr package installation issue; here is the way as suggested by Lawrence at packages github link but dev tree
#install.packages("devtools") # if you don't already have devtools
#devtools::install_github("https://github.com/Columbia-PRIME/pcpr/tree/dev") # you can also replace dev with any other branch you'd like to install

#c_01_estimate_community_sev_index_nyc worked fine and is clearly understandable

#c_02_estimate_community_sev_index_nyc_sens_anal is also working fine and understandable

#still aims of each scripts are missing. I think you are right person to add 

# I don't see table 2 generation code; Jaime, its straight forward summary function but maybe I missed it in scripts

#Note: I had issues in running pre-processing scripts but mostly issues were with directories, which is very small fix.
