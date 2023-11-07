# script aim: 
# First step to load packages etc.
# 1a Declare root directory, folder locations and load essential stuff
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'init_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# set coordinate reference system
crs <- 2163


# read data inputs
traffic_count_2_grid_sld_nyc <- readRDS(paste0(generated.data.folder, "traffic_count_2_grid_sld_nyc.rds"))
traffic_segment_2_grid_sld_nyc <- readRDS(paste0(generated.data.folder, "traffic_segment_2_grid_sld_nyc.rds"))
road_inf_dist_2_grid_sld_nyc <- readRDS(paste0(generated.data.folder, "road_inf_dist_2_grid_sld_nyc.rds"))
barrier_factor_nyc <- readRDS(paste0(generated.data.folder, "barrier_sp_units_nyc_for_use.rds"))
traffic_co2_emis_nyc <- readRDS(paste0(generated.data.folder, "traffic_co2_emis_nyc.rds"))

# load grids and smart location dataset subset
data_desc <- readRDS(paste0(generated.data.folder, "smart_location_data_subset_desc.rds"))
sld_us_loc <- readRDS(paste0(generated.data.folder, "smart_location_data_subset.rds"))
grid <- sf::st_centroid(sld_us_loc) %>%
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
sld_us_loc_df <- grid_contxt
sf::st_geometry(sld_us_loc_df) <- NULL
# subset data for community severance index estimation
sld_us_loc_df <- sld_us_loc_df[,c("GEOID20",  
                                  "D3AAO", "D3APO", 
                                   "D3B", "D3BAO", "NatWalkInd")]
colnames(sld_us_loc_df)[c(2:5)] <- c("autom_netw_dens", "pedest_netw_dens", "street_no_autom_inters_dens", "autom_inters_dens")

sld_us_loc_df <- as.data.frame(sld_us_loc_df)
sld_us_loc_df <- dplyr::na_if(sld_us_loc_df, -99999)

# join data
traffic_count_2_grid_sld_nyc_df <- traffic_count_2_grid_sld_nyc
sf::st_geometry(traffic_count_2_grid_sld_nyc_df) <- NULL
traffic_count_2_grid_sld_nyc_df <- traffic_count_2_grid_sld_nyc_df[,c("GEOID20", "aadt")]
colnames(traffic_count_2_grid_sld_nyc_df) <- c("GEOID20", "aadt_esri_point")

traffic_segment_2_grid_sld_nyc_df <- traffic_segment_2_grid_sld_nyc
sf::st_geometry(traffic_segment_2_grid_sld_nyc_df) <- NULL
traffic_segment_2_grid_sld_nyc_df <- traffic_segment_2_grid_sld_nyc_df[,c("GEOID20", "aadt")]
colnames(traffic_segment_2_grid_sld_nyc_df) <- c("GEOID20", "aadt_fhwa_segm")


barrier_factor_nyc <- barrier_factor_nyc[,c("GEOID20", "barrier_factor_osm", "barrier_factor_fhwa")]

road_inf_dist_2_grid_sld_nyc_df <- road_inf_dist_2_grid_sld_nyc
sf::st_geometry(road_inf_dist_2_grid_sld_nyc_df) <- NULL

data_in_cs <- dplyr::left_join(sld_us_loc_df, traffic_count_2_grid_sld_nyc_df, by = "GEOID20") %>%
  dplyr::left_join(traffic_segment_2_grid_sld_nyc_df, by = "GEOID20") %>%
  dplyr::left_join(road_inf_dist_2_grid_sld_nyc_df, by = "GEOID20") %>%
  dplyr::left_join(traffic_co2_emis_nyc, by = "GEOID20") %>%
  dplyr::left_join(barrier_factor_nyc, by = "GEOID20")

# manuscript Table 2
# explore summary descriptive
summ_data_in_cs <- sumtable(data_in_cs, out = "return")
colnames(data_desc)[1] <- "Variable"
summ_data_in_cs <- dplyr::left_join(summ_data_in_cs, data_desc, by = "Variable")

# homogenize data ranges by scaling by standard deviation
data_in_cs_id <- data_in_cs[,c("GEOID20")]
dta <- data_in_cs[ , -c(which(colnames(data_in_cs) == "GEOID20"))]
dta_scaled = as.data.frame(apply(dta, 2, function(a) a/sd(a, na.rm = T)))
dta_prep <- cbind(data_in_cs_id, dta_scaled)
colnames(dta_prep)[1] <- "GEOID20"

# build dataframe for distributional characteristics (paper table)
# delete variables distance to road and add proximity based

vtable::st(dta, add.median = T,fit.page = '\\textwidth', digits = 2, out = 'latex')


# save data
dta_prep <- cbind(data_in_cs_id, dta_scaled)
colnames(dta_prep)[1] <- "GEOID20"
saveRDS(dta_prep, paste0(generated.data.folder, "community_severance_nyc_input_data.rds"))
