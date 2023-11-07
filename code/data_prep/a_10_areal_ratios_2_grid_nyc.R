# script aim: 
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
rm(sld_us_loc)
# load testing areas
# obtain a spatial context for an example
boroughs <- sf::st_read(paste0(demography.data.folder, "nybb.shp")) %>%
  sf::st_transform(crs)

nyc_boundaries <- sf::st_union(boroughs)

spatial_context <- nyc_boundaries
rm(nyc_boundaries, boroughs)

# limit to spatial context
grid_id_cntxt <- sapply(sf::st_intersects(grid, spatial_context),function(x){length(x)>0})

grid_contxt <- grid[grid_id_cntxt, ]
rm(grid)
grid_contxt$id_local <- 1:nrow(grid_contxt)

# load sidewalks 

edges <- sf::read_sf(paste0(geometry.data.folder, "edges.shp"))
# leave out incorrect entries
edges$width[which(edges$width < 0)] <- NA 
edges$width[which(edges$width > 10)] <- NA
edges_nyc <- edges[-which(is.na(edges$width)),]
edges_nyc <- edges_nyc %>%
  sf::st_transform(crs)
# estimate area of sidewalk 
edges_nyc$sidewalk_area_m <- edges_nyc$width * edges_nyc$length

# load nyc roads in order to estimate area occupied by roads 
street_nyc <- sf::read_sf(paste0(geometry.data.folder, 'geo_export_a18d2619-e5c3-4105-bcf5-99ad8f35058e.shp'))
roads_contxt <- street_nyc %>%
  sf::st_transform(4326) %>% # intermediate step needed to make it work :)
  sf::st_transform(crs)
# transform measures to meters
roads_contxt$shape_leng_m <- roads_contxt$shape_leng * 0.3048 
roads_contxt$st_width_m <- roads_contxt$st_width * 0.3048 
roads_contxt$st_area <- roads_contxt$shape_leng_m * roads_contxt$st_width_m
roads_contxt$bike_lane <- as.numeric(roads_contxt$bike_lane)
# assuming 6 feet width bike lane
roads_contxt$bike_lane_area <- roads_contxt$bike_lane * 6 * 0.3048

# running a 10% each time for making the run easier
n_groups <- 10
group_1 <- 1:as.integer(nrow(grid_contxt)/10)
group_2 <- group_1 + as.integer(nrow(grid_contxt)/10)
group_3 <- group_2 + as.integer(nrow(grid_contxt)/10)
group_4 <- group_3 + as.integer(nrow(grid_contxt)/10)
group_5 <- group_4 + as.integer(nrow(grid_contxt)/10)
group_6 <- group_5 + as.integer(nrow(grid_contxt)/10)
group_7 <- group_6 + as.integer(nrow(grid_contxt)/10)
group_8 <- group_7 + as.integer(nrow(grid_contxt)/10)
group_9 <- group_8 + as.integer(nrow(grid_contxt)/10)
group_10 <- c(group_9 + as.integer(nrow(grid_contxt)/10), nrow(grid_contxt))

areal_road_ped_ratio_coll <- data.frame()
for(grp_n in 1:10){
grp <- get(paste0("group_", grp_n))
print(paste0("working on group ", grp_n))
sp_unit_df <- as.data.frame(grid_contxt[grp,])
sp_unit_df$sp_unit_pos <- grp
sp_unit_ap <- as.array(unlist(sp_unit_df[,c("sp_unit_pos", "id_local")]))
dim(sp_unit_ap) <- c(sp_unit_pos = nrow(sp_unit_df), col = 2)
# function estim_areal_road_ped_ratio can be found at functions.R
areal_road_ped_ratio <- multiApply::Apply(list(sp_unit_ap), target_dims  = 'col',
                                       fun = estim_areal_road_ped_ratio, ncores = config$geom_prec_n_cores)$output1
areal_road_ped_ratio <- t(areal_road_ped_ratio)
colnames(areal_road_ped_ratio) <- c("id_local", "streets_local_area", "roads_local_area", "bikes_local_area", "sidewalk_local_area", "sidewalk_to_road", "sidewalk_to_street","road_to_street") # "barrier_sp_units",
areal_road_ped_ratio_coll <- rbind(areal_road_ped_ratio_coll, areal_road_ped_ratio)
if (grp_n < 10){
saveRDS(areal_road_ped_ratio_coll, paste0(generated.data.folder, "areal_ratio_partial_", "group_", grp_n, ".rds"))
} else {
  saveRDS(areal_road_ped_ratio_coll, paste0(generated.data.folder, "areal_ratio_all.rds"))
}
rm(areal_road_ped_ratio, sp_unit_ap)
}
rm(roads_contxt)

areal_ratio_nyc  <- readRDS(paste0(generated.data.folder, "areal_ratio_all", ".rds"))
areal_ratio_nyc$road_to_street[which(areal_ratio_nyc$road_to_street < 0)] <- NA 

sld_us_loc <- readRDS(paste0(generated.data.folder, "smart_location_data_subset.rds"))
grid <- sf::st_centroid(sld_us_loc[,c("GEOID20")]) %>%
  sf::st_transform(config$crs)
rm(sld_us_loc)
# load testing areas
# obtain a spatial context for an example
boroughs <- sf::st_read(paste0(config$demography_path, "us/nyc/nybb_21a/nybb.shp")) %>%
  sf::st_transform(config$crs)

nyc_boundaries <- sf::st_union(boroughs)

spatial_context <- nyc_boundaries
rm(nyc_boundaries, boroughs)
# subset to spatial context
grid_id_cntxt <- sapply(sf::st_intersects(grid, spatial_context),function(x){length(x)>0})
grid_contxt <- grid[grid_id_cntxt, ]
grid_contxt$id_local <- 1:nrow(grid_contxt)

grid_contxt_df <- grid_contxt
sf::st_geometry(grid_contxt_df) <- NULL


areal_ratio_nyc <- dplyr::left_join(areal_ratio_nyc, grid_contxt_df[, c("id_local", "GEOID20")])
areal_ratio_nyc <- areal_ratio_nyc[,c("GEOID20", "sidewalk_to_road", "sidewalk_local_area", "roads_local_area", "road_to_street")]
saveRDS(areal_ratio_nyc, paste0(generated.data.folder, "areal_ratio_nyc", ".rds"))
