# script aim: generate road infrastructure input data
# First step to load packages etc.
# 1a Declare root directory, folder locations and load essential stuff
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'init_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# set coordinate reference system
crs <- 2163

# load grids
# load grids
sld_us_loc <- readRDS(paste0(generated.data.folder, "smart_location_data_subset.rds"))
grid <- sf::st_centroid(sld_us_loc[,c("GEOID20")]) %>%
  sf::st_transform(crs)
# load testing areas
# obtain a spatial context for an example
boroughs <- sf::st_read(paste0(demography.data.folder, "nybb.shp")) %>%
  sf::st_transform(crs)

nyc_boundaries <- sf::st_union(boroughs)
spatial_context <- nyc_boundaries

# read urban grid
grid_id_cntxt <- sapply(sf::st_intersects(grid, spatial_context),function(x){length(x)>0})
grid_contxt <- grid[grid_id_cntxt, ]

# road infrastructure from faf 5 model - Freight Analysis Framework 5.0 Model Network Database
faf5_network <- sf::read_sf(paste0(geometry.data.folder, "FAF5Network.gdb"))
faf5_highways <- faf5_network[which(faf5_network$F_Class %in% c(1,2,3)),] # Interstate, Principal Arterial - Other Freeways and Expressways, and Principal Arterial - Other)
faf5_highways <- faf5_highways %>%
  sf::st_transform(crs)
faf5_highways_nyc <- faf5_highways[spatial_context,]
# road infrastructure generated at community_severance_main_us
osm_driving_network <- readRDS(paste0(generated.data.folder, "osm_driving_network_northeast.rds"))
# extract driving network following Larkin et al., (2017) Major roads were derived from OSM motorways, motorway links, trunks, trunk links, primary and secondary roads and links
os_highway_roads <- c(
  #Major roads were derived from OSM motorways, motorway links, trunks, trunk links, primary and secondary roads and links.
  'motorway', 'motorway_link', 'trunk', 'trunk_link', 'primary', 'secondary', 'primary_link', 'secondary_link',
  #Minor roads were derived from OSM tertiary roads and tertiary road links.
  'tertiary', 'tertiary_link', 'unclassified',
  #Residential roads were derived from OSM residential roads and residential road links
  'residential'
)
roads_ne_us <- osm_driving_network[which(osm_driving_network$highway %in% os_highway_roads),]
roads_ne_us <- sf::st_transform(roads_ne_us, config$crs)
roads_ne_us_nyc <- roads_ne_us[spatial_context,]
# motorway dist
cs_category <- c("motorway", "motorway", 'trunk', 'trunk', "primary", "primary", "secondary", "secondary", "tertiary", "tertiary", "residential", "residential", "residential", "residential")
osm_category <- c("motorway", "motorway_link", 'trunk', 'trunk_link', "primary", "primary_link", "secondary", "secondary_link", 'tertiary', 'tertiary_link', "residential", "residential_link", 'unclassified', 'unclassified_link')
road_types <- data.frame(cs_category = cs_category, osm_category = osm_category)
cats <- unique(road_types$cs_category)
for(c in 1:length(cats)){
cat <- cats[c]
osm_types <- road_types[which(road_types$cs_category == cat), "osm_category"]
roads_ne_us_loc <- roads_ne_us[which(roads_ne_us$highway %in% osm_types),]
if(nrow(roads_ne_us_loc) > 0){
roads_ne_us_loc_id_cntxt <- sapply(sf::st_intersects(roads_ne_us_loc, spatial_context),function(x){length(x)>0})

roads_ne_us_loc_test <- roads_ne_us_loc[roads_ne_us_loc_id_cntxt, ]

dist <- sf::st_distance(grid_contxt, roads_ne_us_loc_test)
var_name <- paste0(cat, "_dist")
grid_contxt[,var_name] <- apply(dist, 1, FUN = min, na.rm = TRUE)
} else {
grid_contxt[,var_name] <- NA
}
}

# distance to major roads according to faf 5
# FHWA highway functional class designation:
#   1 - Interstate
# 2 - Principal Arterial - Other Freeways and Expressways
# 3 - Principal Arterial - Other
cs_category <- c("interstate_highway", "freeways_expressways", 'other_princ_arter')
fhwa_category <- c(1, 2, 3)
road_types <- data.frame(cs_category = cs_category, fhwa_category = fhwa_category)
cats <- unique(road_types$cs_category)
for(c in 1:length(cats)){
  cat <- cats[c]
  fhwa_types <- road_types[which(road_types$cs_category == cat), "fhwa_category"]
  roads_loc <- faf5_highways[which(faf5_highways$F_Class %in% fhwa_types),]
  if(nrow(roads_loc) > 0){
    roads_loc_id_cntxt <- sapply(sf::st_intersects(roads_loc, spatial_context),function(x){length(x)>0})
    
    roads_loc_test <- roads_loc[roads_loc_id_cntxt, ]
    
    dist <- sf::st_distance(grid_contxt, roads_loc_test)
    var_name <- paste0(cat, "_dist")
    grid_contxt[,var_name] <- apply(dist, 1, FUN = min, na.rm = TRUE)
  } else {
    grid_contxt[,var_name] <- NA
  }
}


# distance metrics to proximity (reverse direction)
max_motorway <- max(data_in_cs$motorway_dist) 
data_in_cs$motorway_prox <- (max_motorway - data_in_cs$motorway_dist) / max_motorway

# checking if abouve method is consistent with normalization below
normalize <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
test <- data_in_cs
test$motorway_prox_2 <- 1- normalize(test$motorway_dist)
summary(test[,c("motorway_prox", "motorway_prox_2")])
# yes, it is,same result

max_primary <- max(grid_contxt$primary_dist) 
grid_contxt$primary_prox <- (max_primary - grid_contxt$primary_dist) / max_primary

max_secondary <- max(grid_contxt$secondary_dist) 
grid_contxt$secondary_prox <- (max_secondary - grid_contxt$secondary_dist) / max_secondary

max_tertiary <- max(grid_contxt$tertiary_dist) 
grid_contxt$tertiary_prox <- (max_tertiary - grid_contxt$tertiary_dist) / max_tertiary

max_residential <- max(grid_contxt$residential_dist) 
grid_contxt$residential_prox <- (max_residential - grid_contxt$residential_dist) / max_residential

max_trunk <- max(grid_contxt$trunk_dist) 
grid_contxt$trunk_prox <- (max_trunk - grid_contxt$trunk_dist) / max_trunk

max_interstate_highway <- max(grid_contxt$interstate_highway_dist) 
grid_contxt$interstate_highway_prox <- (max_interstate_highway - grid_contxt$interstate_highway_dist) / max_interstate_highway

max_freeways_expressways <- max(grid_contxt$freeways_expressways_dist) 
grid_contxt$freeways_expressways_prox <- (max_freeways_expressways - grid_contxt$freeways_expressways_dist) / max_freeways_expressways

max_other_princ_arter <- max(grid_contxt$other_princ_arter_dist)
grid_contxt$other_princ_arter_prox <- (max_other_princ_arter - grid_contxt$other_princ_arter_dist) / max_other_princ_arter


grid_contxt <- grid_contxt[,-which(colnames(grid_contxt) %in% c("motorway_dist", "primary_dist", "secondary_dist", "tertiary_dist", "residential_dist", "trunk_dist",
                                                             "interstate_highway_dist", "freeways_expressways_dist", "other_princ_arter_dist"))]

saveRDS(grid_contxt, paste0(generated.data.folder, "road_inf_dist_2_grid_sld_nyc.rds"))

