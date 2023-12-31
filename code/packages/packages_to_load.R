

# list of packages to use
list.of.packages = c('PCPhelpers', 'pcpr', 'ggridges', 'foreach', 'tictoc', 'progressr', 'httr', 'RCurl', 'vtable', 'osmextract', 'esri2sf', 'automap', 'tidyverse', 'haven','ggplot2','GGally','caret','readr',
                     'fBasics','raster', 'nortest','MASS','extrafont', 'forestplot',
                     'ggsn','survival','grid','gridExtra', 'cowplot','reshape2', 'ggrepel',
                     'naniar','rpart','rpart.plot','pROC', 'ggfortify', 'psych', 'GPArotation',
                     'car','olsrr','EnvStats','janitor','tableone', 'reshape', 'kableExtra',
                     'sf', 'areal', 'tmap', 'mgcv','gamm4', 'plotly', 'viridis', 'rgeos', 
                     'spdep', 'gstat', 'splm', 'readxl', 'plm', 'itsadug')

# check if list of packages is installed. If not, it will install ones not yet installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
lapply(list.of.packages, require, character.only = TRUE)




