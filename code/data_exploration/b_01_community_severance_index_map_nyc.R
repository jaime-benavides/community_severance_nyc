# script aim: 
# First step to load packages etc.
# 1a Declare root directory, folder locations and load essential stuff
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'init_directory_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# set coordinate reference system
crs <- 2163

sld_us_loc <- readRDS(paste0(generated.data.folder, "smart_location_data_subset.rds"))
dat_scores <- readRDS(paste0(generated.data.folder, "comm_sev_fa_scores_nyc_dta_us.rds"))

sld_us_loc <- sld_us_loc%>%
  sf::st_transform(crs)

grid <- sf::st_centroid(sld_us_loc) 

# load spatial context (NYC)
boroughs <- sf::st_read(paste0(demography.data.folder, "nybb.shp")) %>%
  sf::st_transform(crs)
nyc_boundaries <- sf::st_union(boroughs)
spatial_context <- nyc_boundaries

# create grid for nyc
grid_id_cntxt <- sapply(sf::st_intersects(grid, spatial_context),function(x){length(x)>0})
grid_contxt <- grid[grid_id_cntxt, ]
sld_us_loc_df <- grid_contxt
sf::st_geometry(sld_us_loc_df) <- NULL

# add community severance index
grid_contxt <- dplyr::left_join(grid_contxt[,c("GEOID20", "TotPop", "NatWalkInd")], dat_scores[,c("GEOID20","MR1_norm")], by = "GEOID20")
grid_contxt <- grid_contxt[,c("GEOID20", "TotPop", "NatWalkInd", "MR1_norm")]
colnames(grid_contxt)[which(colnames(grid_contxt) == "MR1_norm")] <- "community_severance_index"

# back to dataframe
grid_contxt_df <- grid_contxt
sf::st_geometry(grid_contxt_df) <- NULL

# 1 rank census block groups from lower to higher community severance index
comm_sev <- grid_contxt_df

comm_sev_sf <- dplyr::left_join(comm_sev, sld_us_loc[,"GEOID20"], by = "GEOID20")
comm_sev_sf <- sf::st_as_sf(comm_sev_sf)
comm_sev_sf <- sf::st_transform(comm_sev_sf, crs) 

# leave out areas over water and parks for plotting
comm_sev_sf <- sf::st_intersection(comm_sev_sf,spatial_context)
park_geiods <- sld_us_loc_df[which(sld_us_loc_df$TotPop < 20),"GEOID20"]
comm_sev_sf_p <- comm_sev_sf[-which(comm_sev_sf$GEOID20 %in% park_geiods$GEOID20),]

# read road infrastructure for plotting
faf5_network <- sf::read_sf(paste0(geometry.data.folder, "FAF5_Model_Highway_Network/Networks/geodatabase_format/FAF5Network.gdb"))
faf5_highways <- faf5_network[which(faf5_network$F_Class %in% c(1,2,3)),]
faf5_highways <- faf5_highways %>%
  sf::st_transform(crs)

# intersect with spatial context
roads_contxt_id_fhwa <- sapply(sf::st_intersects(faf5_highways, spatial_context),function(x){length(x)>0})
roads_contxt_fhwa <- faf5_highways[roads_contxt_id_fhwa, ]
rm(faf5_highways)

# prepare map of community severance index
# manuscript Figure 5
tmap_mode(mode = "view")
tmap_options(check.and.fix = TRUE)
map_community_severance_ctxt <-  tm_shape(sf::st_make_valid(spatial_context)) +
  tm_borders(alpha = 0.1, col = "black") +
  tm_shape(comm_sev_sf_p) +
  tm_polygons("community_severance_index", palette = 'Reds', scale=2, style = "order",
              title = "Community Severance Index") +
  tm_shape(roads_contxt_fhwa) +
  tm_lines(alpha = 0.5)
# 
tmap_save(map_community_severance_ctxt, paste0(output.folder, "map_community_severance_index_nyc.html"))

# read barrier factors for plotting

barrier_factors <- readRDS(paste0(generated.data.folder, "barrier_sp_units_nyc_for_use.rds"))
comm_sev_sf_p <- dplyr::left_join(comm_sev_sf_p, barrier_factors)

# manuscript Figure 2
map_barrier_factor_ctxt <-  tm_shape(sf::st_make_valid(spatial_context)) +
  tm_borders(alpha = 0.1, col = "black") +
  tm_shape(comm_sev_sf_p) +
  tm_polygons("barrier_factor_fhwa", palette = 'Oranges', scale=2, style = "order",
              title = "Barrier Factor") +
  tm_shape(roads_contxt_fhwa) +
  tm_lines(alpha = 0.5)

tmap_save(map_barrier_factor_ctxt, paste0(output.folder, "map_barrier_factor_nyc.html"))

