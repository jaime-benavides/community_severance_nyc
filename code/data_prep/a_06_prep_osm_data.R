# script aim: 
# First step to load packages etc.
# 1a Declare root directory, folder locations and load essential stuff
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'init_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# set coordinate reference system
crs <- 2163

### load data

us_cont_bbox <- sf::st_as_sfc(sf::st_bbox(c(xmin = -130, xmax = -60, ymax = 51, ymin = 21), crs = 4326)) 
test_bbox <- sf::st_as_sfc(sf::st_bbox(c(xmin = -73.8, xmax = -74.1, ymax = 41.2, ymin = 40.4), crs = 4326)) 
us_boundaries <- sf::st_read(paste0(config$demography_path, "us/nation/cb_2018_us_nation_20m.shp")) %>%
  sf::st_transform(config$crs) %>%
  sf::st_intersection(us_cont_bbox) 

boroughs <- sf::st_read(paste0(demography.data.folder, "nybb.shp")) %>%
  sf::st_transform(crs) %>%
  sf::st_union()

## road characteristics (make query for osm data)

driving_network_v_t_opts <- c(
  "-where", "
    (highway IS NOT NULL)
    AND
    (highway NOT IN (
    'abandoned', 'bus_guideway', 'byway', 'construction', 'corridor', 'elevator',
    'fixme', 'escalator', 'gallop', 'historic', 'no', 'planned', 'platform',
    'proposed', 'cycleway', 'pedestrian', 'bridleway', 'path', 'footway',
    'steps'
    ))
    AND
    (access NOT IN ('private', 'no'))
    AND
    (service NOT ILIKE 'private%')
    "
)
driving_network_ext_tgs <- c("lanes", "maxspeed", "access", "service", "barrier", "surface", "tiger:cfcc", "parking:lane:both", "parking:lane:left", "parking:lane:right")

pbf = file.path(paste0(geometry.data.folder, "us-northeast-latest.osm.pbf"))

# road network and parking data call 
osmextract::oe_vectortranslate(
  pbf,
  layer = "lines",
  vectortranslate_options = driving_network_v_t_opts,
  osmconf_ini = NULL,
  extra_tags = driving_network_ext_tgs,
  force_vectortranslate = TRUE,
  never_skip_vectortranslate = FALSE,
  boundary = NULL,
  boundary_type = c("spat", "clipsrc"),
  quiet = FALSE
)
# storing the just saved file in a known path and delete unused variables
osm_driving_network <- osmextract::oe_read(paste0(raw.data.folder, "us-northeast-latest.gpkg"))
unused_vars_ind <- which(colnames(osm_driving_network) %in% c("waterway", "aerialway", "man_made"))
osm_driving_network <- osm_driving_network[,-unused_vars_ind]
osm_driving_network <- saveRDS(osm_driving_network, paste0(generated.data.folder, "osm_driving_network_northeast.rds"))