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

# add covariates
# ses - Area Deprivation Index (ADI)
cbg_adi <- readr::read_csv(demography.data.folder, "US_2020_ADI_Census_Block_Group_v4_0_1.csv")
colnames(cbg_adi)[2] <- "GEOID20"

# motor vehicle collisions

crashes <- readr::read_csv(paste0(crashes_path, "Motor_Vehicle_Collisions_-_Crashes.csv"))
# subset to those with coordinates
crashes_loc <- crashes[which(complete.cases(crashes[,c("LONGITUDE", "LATITUDE")])), ]
# dates
crashes_loc$`CRASH DATE` <- as.Date(crashes_loc$`CRASH DATE`, "%m/%d/%Y")
# subset crashes from 2019
colnames(crashes_loc)[1] <- "date"
crashes_loc <- crashes_loc[crashes_loc$date >= as.Date("2019-01-01") & crashes_loc$date <= as.Date("2019-12-31"), ]
# build spatial object with crashes
crashes_loc_sf <- sf::st_as_sf(crashes_loc, coords = c("LONGITUDE","LATITUDE"  ), crs = 4326)  %>%
  sf::st_transform(crs) 

# create subset of crashes with cyclists or pedestrians involved
crashes_loc_byc_ped_sf <- crashes_loc_sf[which(
  crashes_loc_sf$`NUMBER OF PEDESTRIANS INJURED`>0  | 
    crashes_loc_sf$`NUMBER OF PEDESTRIANS KILLED` > 0 |
    crashes_loc_sf$`NUMBER OF CYCLIST INJURED` > 0  |
    crashes_loc_sf$`NUMBER OF CYCLIST KILLED` > 0),]

# prepare smart location dataset
sld_us_loc <- sld_us_loc%>%
  sf::st_transform(crs)


# intersect crashes with census block groups
crashes_loc_sf <- sf::st_intersection(crashes_loc_sf, sld_us_loc[,"GEOID20"])
# back to dataframe
crashes_loc_df <- crashes_loc_sf
sf::st_geometry(crashes_loc_df) <- NULL
# count crashes
crashes_bg <- crashes_loc_df %>% dplyr::count(GEOID20)
# count crashes for pedestrians/cyclists injured or killed
crashes_loc_byc_ped_sf <- sf::st_intersection(crashes_loc_byc_ped_sf, sld_us_loc[,"GEOID20"])
crashes_loc_byc_ped_df <- crashes_loc_byc_ped_sf
sf::st_geometry(crashes_loc_byc_ped_df) <- NULL
crashes_byc_ped <- crashes_loc_byc_ped_df %>% dplyr::count(GEOID20)

# obtain centroids
grid <- sf::st_centroid(sld_us_loc[,c("GEOID20", "TotPop", "Ac_Total", "NatWalkInd")]) 
# spatial context
boroughs <- sf::st_read(paste0(demography.data.folder, "nybb.shp")) %>%
  sf::st_transform(crs)
nyc_boundaries <- sf::st_union(boroughs)
spatial_context <- nyc_boundaries


# prepare dataset for regression models
# intersect data with spatial context (nyc)
grid_id_cntxt <- sapply(sf::st_intersects(grid, spatial_context),function(x){length(x)>0})
grid_contxt <- grid[grid_id_cntxt, ]
sld_us_loc_df <- grid_contxt
sf::st_geometry(sld_us_loc_df) <- NULL
# put together data
grid_contxt <- dplyr::left_join(grid_contxt[,c("GEOID20", "TotPop", "Ac_Total", "NatWalkInd")], dat_scores[,c("GEOID20","MR1_norm")], by = "GEOID20")
grid_contxt <- grid_contxt[,c("GEOID20", "TotPop", "Ac_Total","NatWalkInd", "MR1_norm")]
# name community severance index
colnames(grid_contxt)[which(colnames(grid_contxt) == "MR1_norm")] <- "community_severance_index"
# add area deprivation index
grid_contxt <- dplyr::left_join(grid_contxt, cbg_adi[,c("GEOID20", "ADI_NATRANK")], by = "GEOID20")
# get coordinates for assessing spatial autocorrelation
grid_contxt_cent <- data.frame(GEOID20 = grid_contxt$GEOID20, sf::st_coordinates(sf::st_transform(grid_contxt, "WGS84")))
colnames(grid_contxt_cent)[2:3] <- c("lon", "lat")

# compare community severance to walkability index
grid_contxt_df <- grid_contxt
sf::st_geometry(grid_contxt_df) <- NULL
# add population density
grid_contxt_df$pop_dens <- grid_contxt_df$TotPop/grid_contxt_df$Ac_Total
# add total crashes
grid_contxt_df <- dplyr::left_join(grid_contxt_df, crashes_bg, by = "GEOID20")
colnames(grid_contxt_df)[which(colnames(grid_contxt_df)== "n")] <- "traffic_incidents"
# add crashes involving pedestrians and cyclists
grid_contxt_df <- dplyr::left_join(grid_contxt_df, crashes_byc_ped, by = "GEOID20")
colnames(grid_contxt_df)[which(colnames(grid_contxt_df)== "n")] <- "traffic_incidents_ped_byc"

# scale variables by standard deviation
grid_contxt_df$ADI_NATRANK <- as.numeric(grid_contxt_df$ADI_NATRANK)
grid_contxt_df <- dplyr::left_join(grid_contxt_df, grid_contxt_cent, by = "GEOID20")
scaled <- grid_contxt_df %>% 
  dplyr::select(pop_dens, ADI_NATRANK, lon, lat) %>% 
  purrr::map_df(~ scale(.x, center = FALSE, scale = sd(.x, na.rm = TRUE)))
# dataset for models
data_scale <- cbind(grid_contxt_df[, c("GEOID20", "community_severance_index", "NatWalkInd", "traffic_incidents", "traffic_incidents_ped_byc")],
                    scaled[,c("pop_dens", "ADI_NATRANK", "lon", "lat")])

# choosing optimal number of knots

traffic_incidents_csi_k5 <- mgcv::gam(traffic_incidents ~ s(community_severance_index, k=5, fx=TRUE) + 
                                        pop_dens  + ADI_NATRANK, family = nb(1), data = data_scale)
traffic_incidents_csi_k4 <- mgcv::gam(traffic_incidents ~ s(community_severance_index, k=4, fx=TRUE) + 
                                        pop_dens  + ADI_NATRANK, family = nb(1), data = data_scale)
traffic_incidents_csi_k3 <- mgcv::gam(traffic_incidents ~ s(community_severance_index, k=3, fx=TRUE) + 
                                        pop_dens  + ADI_NATRANK, family = nb(1), data = data_scale)
AIC(traffic_incidents_csi_k5, traffic_incidents_csi_k4)
AIC(traffic_incidents_csi_k3, traffic_incidents_csi_k4)

# change reference level (from average to minimum) for plotting regression model results
# need to use the crossbasis from the dlnm package as an input in the gam. and then use crosspredict afterwards!
# estimate minimum community severance index
min_csi <- min(data_scale$community_severance_index)
# create crossbasis for community severance using the spline from above result (df=2)
cb.csi <- dlnm::crossbasis(data_scale$community_severance_index,
                        lag=0,
                        argvar=list(fun = "ns", df = 2))
# check that AIC is consistent when running this adapted model 
traffic_incidents_csi_k3_dlnm <- mgcv::gam(traffic_incidents ~ cb.csi + 
                                        pop_dens  + ADI_NATRANK,  family = nb(1), data = data_scale)
AIC(traffic_incidents_csi_k3, traffic_incidents_csi_k3_dlnm)

# get predictions using this crossbasis and adding the minimum as the reference level
pred.csi <- dlnm::crosspred(
  basis =cb.csi,   
  model=traffic_incidents_csi_k3_dlnm, 
  model.link = "log",
  at = data_scale$community_severance_index, 
  cen = min(data_scale$community_severance_index))


mod_dta <- format_pred_mod(pred.csi)

# manuscript Figure 6
name_plot <- "vehicle_collisions"
png(paste0(output.folder, name_plot,"_model_res.png"), 900, 460)
plot_pred_mod(mod_dta)
dev.off()

# sensitivity analysis to including walkability index in the model as covariate
traffic_incidents_csi_k3_dlnm_walk <- mgcv::gam(traffic_incidents ~ cb.csi + 
                                             pop_dens + NatWalkInd + ADI_NATRANK,  family = nb(1), data = data_scale)

pred.csi_walk <- dlnm::crosspred(
  basis =cb.csi,   
  model=traffic_incidents_csi_k3_dlnm_walk, 
  model.link = "log",
  at = data_scale$community_severance_index, 
  cen = min(data_scale$community_severance_index))

mod_dta <- format_pred_mod(pred.csi_walk)

# manuscript Figure S3
name_plot <- "vehicle_collisions_walk"
png(paste0(output.folder, name_plot,"_model_res.png"), 900, 460)
plot_pred_mod(mod_dta)
dev.off()

# sensitivity analysis to traffic incidents involving pedestrians cyclists 

traffic_incidents_ped_byc_csi_k3_dlnm <- mgcv::gam(traffic_incidents_ped_byc ~ cb.csi + 
                                                          pop_dens + ADI_NATRANK,  family = nb(1), data = data_scale)

pred.ped_byc_csi <- dlnm::crosspred(
  basis =cb.csi,   
  model=traffic_incidents_ped_byc_csi_k3_dlnm, 
  model.link = "log",
  at = data_scale$community_severance_index, 
  cen = min(data_scale$community_severance_index))

mod_dta <- format_pred_mod(pred.ped_byc_csi)

# manuscript Figure S4
name_plot <- "vehicle_collisions_ped_byc"
png(paste0(output.folder, name_plot,"_model_res.png"), 900, 460)
plot_pred_mod(mod_dta)
dev.off()

# sensitivity analysis to traffic incidents involving pedestrians cyclists and walkability index as covariate


traffic_incidents_ped_byc_csi_k3_dlnm_walk <- mgcv::gam(traffic_incidents_ped_byc ~ cb.csi + 
                                                          pop_dens + NatWalkInd + ADI_NATRANK,  family = nb(1), data = data_scale)

pred.ped_byc_csi_walk <- dlnm::crosspred(
  basis =cb.csi,   
  model=traffic_incidents_ped_byc_csi_k3_dlnm_walk, 
  model.link = "log",
  at = data_scale$community_severance_index, 
  cen = min(data_scale$community_severance_index))

mod_dta <- format_pred_mod(pred.ped_byc_csi_walk)

# manuscript Figure S5
name_plot <- "vehicle_collisions_ped_byc_walk"
png(paste0(output.folder, name_plot,"_model_res.png"), 900, 460)
plot_pred_mod(mod_dta)
dev.off()


# adding lat and long as covariates to check spatial autocorrelation

traffic_incidents_csi_k3_dlnm_sp <- mgcv::gam(traffic_incidents ~ cb.csi + 
                                                  pop_dens + ADI_NATRANK + t2(lat,lon),  family = nb(1), data = data_scale)

pred.csi_sp <- dlnm::crosspred(
  basis =cb.csi,   
  model=traffic_incidents_csi_k3_dlnm_sp, 
  model.link = "log",
  at = data_scale$community_severance_index, 
  cen = min(data_scale$community_severance_index))

mod_dta <- format_pred_mod(pred.csi_sp)

# manuscript Figure S6
name_plot <- "vehicle_collisions_sp"
png(paste0(output.folder, name_plot,"_model_res.png"), 900, 460)
plot_pred_mod(mod_dta)
dev.off()