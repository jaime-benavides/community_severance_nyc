## Code Review 
## Development of a Community Severance Index for urban areas in the United States: A case study in New York City
## Vijay Kumar 
## Updated: Nov 16, 2023

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

#packages fixed! #after talking to Jaime, fixed pcpr package installation issue; here is the way as suggested by Lawrence at packages github link but dev tree
#install.packages("devtools") # if you don't already have devtools
#devtools::install_github("https://github.com/Columbia-PRIME/pcpr/tree/dev") # you can also replace dev with any other branch you'd like to install

#Comment on each RScript

#1. functions.R: 
#script line 1: for consistency script aim is missing (I know it's already clear from name of file but I just added a line), 
#script line 3: also function number is added

#2. script_initiate.R
#Script line 1: script aim: This script loads required packages as well as all functions used in this study 

#also added notes about packages installations importantly pcpr

#script line 8-13: you might face issue on installing libraries like PCPhelpers and pcpr as well as "esri2sf". Please use following code to install them from github repositories

# library(devtools)
# install_github("Columbia-PRIME/PCPhelpers")
#devtools::install_github("https://github.com/Columbia-PRIME/pcpr/tree/dev")
#install_github("yonghah/esri2sf") 

#3. init_directory_structure.R
#script line 1: #script aim: this script creates directories required for the code structure of this study
#script line 34: we might need a path for reading car crashes file in c_03 file line 21 "crashes <- readr::read_csv(paste0(traffic.data.folder, "Motor_Vehicle_Collisions_-_Crashes_20231113.csv"))"

#4. a_01_preproc_smart_location_dta.R
#script line 1: # script aim: This script is preparing spatial data sets for this study from Smart Location Database (Jaime, please add or correct the aims/objective)
#script line 15: "Selection of urban spatial variables, details in table 1 of manuscript (Jaime, you can add more details or delete if it's incorrect)"

#5. a_02_barrier_factor_prep.R:
#no data folder and file roads/osm_driving_network_northeast.rds at line #36
#script line 36: to fix this, data directory has been changed to generated data in script
# I am not changing it into script, assuming it can be fixed easily

#from here, I am using Jaime's preprocessed data for analysis.


#6. b_01_community_serverance_index.R
#script line 1: # script aim is missing
#script line 53: corrected directory by changing code faf5_network <- sf::read_sf(paste0(geometry.data.folder, "FAF5Network.gdb"))
# I didn't change in code as directory can be fixed easily, otherwise code working fine and understandable. 


#7. c_01_estimate_community_sev_index_nyc.R
#worked fine and is clearly understandable

#8. c_02_estimate_community_sev_index_nyc_sens_anal.R
#is also working fine and understandable

#9. c_03_build_csi_models.R
#script line 13: notes on data of car crashes  
#you might need to download car crash data from  nyc open data and change folder name accordingly
 
# I don't see table 2 generation code; Jaime, its straight forward summary function but maybe I missed it in scripts

#Note: I had issues in running pre-processing scripts but mostly issues were with directories, which is very small fix.

#Jaime, aims of each script are missing. I think you are right person to add and also check what I have added.
