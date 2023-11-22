# script aim: interpolate traffic intensity data to census block groups
# First step to load packages etc.
# 1a Declare root directory, folder locations and load essential stuff
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'init_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# set coordinate reference system
crs <- 2163

# load traffic 
aadt_nys_2019 <- readRDS(paste0(paste0(generated.data.folder, "aadt_ny_2019.rds")))
aadt_segments <- aadt_nys_2019 
# load grids
sld_us_loc <- readRDS(paste0(generated.data.folder, "smart_location_data_subset.rds"))
grid <- sf::st_centroid(sld_us_loc[,c("GEOID20")]) %>%
  sf::st_transform(crs)

boroughs <- sf::st_read(paste0(demography.data.folder, "nybb.shp")) %>%
  sf::st_transform(crs)

nyc_boundaries <- sf::st_union(boroughs)
spatial_context <- nyc_boundaries

grid_id_cntxt <- sapply(sf::st_intersects(grid, spatial_context),function(x){length(x)>0})
grid_contxt <- grid[grid_id_cntxt, ]


aadt_segments <- sf::st_transform(aadt_segments, crs = crs)
aadt_segments_id_cntxt <- sapply(sf::st_intersects(aadt_segments, spatial_context),function(x){length(x)>0})
aadt_segments_contxt <- aadt_segments[aadt_segments_id_cntxt, ]
aadt_segments_contxt <- sf::st_cast(aadt_segments_contxt, "LINESTRING")

# select three points per line to represent the street segment in interpolation
n_points <- 3
road_points <- sf::st_transform(aadt_segments_contxt, crs) %>%
  sf::st_line_sample(n = n_points, type = "regular") %>%
  sf::st_cast("POINT") 
aadt_segments_p <- st_sf(aadt = rep(aadt_segments_contxt$aadt, each =n_points), geom = road_points)

# kriging from mid point 
# adapted from Criado et al. (2022) https://earth.bsc.es/gitlab/es/universalkriging/-/blob/production/general/UK_mean.R
UK_mean_uniform_ok <- regrid_ok(non_uniform_data = sf::as_Spatial(aadt_segments_p), # traffic_esri_cntxt
                                target_grid = sf::as_Spatial(grid_contxt), crs_sim = crs) # grid_contxt
colnames(UK_mean_uniform_ok)[1] <- "aadt"
UK_mean_uniform_ok$GEOID20 <- grid_contxt$GEOID20
saveRDS(UK_mean_uniform_ok, paste0(generated.data.folder, "traffic_segment_2_grid_sld_nyc.rds"))
