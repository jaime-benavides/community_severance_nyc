# script aim: 
# First step to load packages etc.
# 1a Declare root directory, folder locations and load essential stuff
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'init_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# traffic aadt esri
url <- "https://demographics5.arcgis.com/arcgis/rest/services/USA_Traffic_Counts/MapServer/0"
tok <- "0E0uMkcA0SwWVbeHW5BNKNh55fkPUkiHHy8e6F7TL-DlIkVGxPYHfMb0mKESeAWPhlEW73kNLAR9PyLtL7xgWeoYji8m7dRCTTckdmuPSBZNEB0iSg6PhFe9EebI18Int2_A4tc6trEZBvu2HVkcDw.."
df <- esri2sf::esri2sf(url, token = tok)
saveRDS(df, paste0(generated.data.folder, "traffic_counts_esri.rds"))

# fhwa traffic segment aadt
url_ny <- "https://geo.dot.gov/server/rest/services/Hosted/HPMS_FULL_NY_2019/FeatureServer/0"
df_ny <- esri2sf(url_ny)
saveRDS(df_ny, paste0(paste0(generated.data.folder, "aadt_ny_2019.rds")))
