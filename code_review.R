## Code Review 
## Development of a Community Severance Index for urban areas in the United States: A case study in New York City
## Vijay Kumar 
## Updated: Nov 6, 2023

# script: line:
# comment: Jaime, can we add chunk of code or details for installing the PCPhelpers and pcpr packages?
#comment:The way I installed instead of downloading them
library(devtools)
install_github("Columbia-PRIME/PCPhelpers") #worked fine 
install_github("Columbia-PRIME/pcpr") #having error, also tried Lawrence suggested way but have same error 
install.packages("pcpr-master", repos = NULL, type="source")
#installing *source* package ‘pcpr’ ...
# ** using staged installation
# ** R
# ** byte-compile and prepare package for lazy loading
# Error: object ‘multiprocess’ is not exported by 'namespace:future'
# Execution halted
# ERROR: lazy loading failed for package ‘pcpr’
# * removing ‘/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/pcpr’
# Warning message:
#   In i.p(...) :
#   installation of package ‘/var/folders/g8/dk9kfnvs4g126rq6lwy_pjhr0000gq/T//Rtmp9sRd3c/filedb4bbdc388a/pcpr_1.0.0.tar.gz’ had non-zero exit status
# I tried following but looks like multipurpose package is not on R
# Reinstall 'future' package
install.packages("future")

# Update 'multiprocess' package if available
install.packages("multiprocess")

# Then try reinstalling 'pcpr'
remotes::install_github("Columbia-PRIME/pcpr")

# Jaime; In general I have two issues so far; 1- I coudn't install pcpr package, and 2- couldn't clone this repository to R; maybe you need to grant me persmission!

#checking "a_01_preproc_smart_location_dta.R" packages issue "esri2sf", "rgeos" and "pcpr"
install_github("yonghah/esri2sf") #worked!
#other two packages still have problem but I could run a_01_preproc_smart_location_dta.R file and generate the data

# added on line 15 of code "#Selection of urban spatial variables, details in table 1 of manuscript"
