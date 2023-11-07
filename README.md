# community_severance_nyc
Repository for reviewing the code of the project Development of a Community Severance Index for urban areas in the United States: A case study in New York City

note: please run init_directory_structure.R first to create folders and add the necessary input data to the raw folder as detailed at the end of this readme in order to run the project. 

# note on packages
need to download PCPhelpers (https://github.com/Columbia-PRIME/PCPhelpers) and pcpr (https://github.com/Columbia-PRIME/pcpr) from github before running the models

## Code and data generated stored at data/generated/ (file name - short description)

### Data preparation (data_prep) list including tables/figures:

a_01_preproc_smart_location_dta.R - preprocess smart location database

- smart_location_data_subset.rds (750Mb) - smart location data subset for each census block group

- smart_location_data_subset_desc.rds - description of smart location data subset for use in this project 

a_02_barrier_factor_prep.R (time consuming: divided in 2 processes) - barrier factor for both osm and faf5 data

- barrier_sp_units_nyc_for_use.rds - barrier factor from both sources for each census block group

# traffic counts from ESRI source 

a_03_prep_traffic.R - download to local machine traffic count data from esri

- traffic_counts_esri.rds - object containing traffic counts from esri

a_04_traffic_count_esri_to_cbg.R - interpolate traffic count esri data to census block group centroids

- traffic_count_2_grid_sld_nyc.rds - object containing traffic count for each census block group centroids

a_05_traffic_segment_hpms_to_cbg.R - interpolate traffic segment hpms data to census block group centroids

- traffic_segment_2_grid_sld_nyc - object containing traffic segment for each census block group centroids

a_06_prep_osm_data.R

- osm_driving_network_northeast.rds - osm roads 

a_07_road_infrastructure_dist_to_cbg.R - estimate distance from each type of road infrastructure to census block group and estimate proximity metric

- road_inf_dist_2_grid_sld_nyc.rds - proximity to different type of road infrastructure for each census block group

a_08_traffic_co2_emissions_to_cbg.R - estimate traffic co2 emissions for each census block group and estimate co2 emissions / area

- traffic_co2_emis_nyc - traffic co2 emissions for each census block group

- a_09_put_together_inputs_to_csi.R - create dataframe containing all the necessary inputs to estimate community severance index

- community_severance_nyc_input_data.rds - input data for community severance index estimation

- Table 2

### Data exploration (data_exploration) list including tables/figures:

b_01_community_severance_index_map_nyc.R - exploratory data analysis of noise complaints

- Figure 2, Figure 5


### Model running (models) list including tables/figures:

c_01_estimate_community_sev_index_nyc.R - estimating community severance index

Table S1, Figure 3a, Figure 3b, Figure 4, Figure S1

c_02_estimate_community_sev_index_nyc_sens_anal.R - estimating community severance index adding as input variable sidewalks

Figure S2

c_03_build_csi_models.R - build regression models

Figure 6, Figure S3, Figure S4, Figure S5, Figure S6

## Data (data) list:

### Raw (description - file name - link to source - where should be stored)

#### smart location database

The Smart Location Database is a nationwide geographic data resource for measuring location efficiency - SmartLocationDatabase.gdb - https://www.epa.gov/smartgrowth/smart-location-mapping#SLD - data/raw/geometry/

#### demography

NYC boroughs - nybb.shp -  https://www1.nyc.gov/site/planning/data-maps/open-data/districts-download-metadata.page - data/raw/demography/

Area deprivation index US - https://www.neighborhoodatlas.medicine.wisc.edu/ US_2020_ADI_Census_Block_Group_v4_0_1.csv - data/raw/demography/

#### road infrastructure

- OSM data 

us-northeast-latest.osm.pbf - downloaded from geofabrik (search nyc) https://download.geofabrik.de/north-america.html - data/raw/geometry/ (1.4Gb)

- FAF 5 Network

FAF5Network.gdb - downloaded from faf 5 model - Freight Analysis Framework 5.0 Model Network Database - https://geodata.bts.gov/datasets/9343414b46794fb8be9867db2d1ccb75/about - data/raw/geometry

#### traffic activity

ESRI traffic counts

traffic intensity data from ESRI database - https://demographics5.arcgis.com/arcgis/rest/services/USA_Traffic_Counts/MapServer/0 - downloaded at a_03_prep_traffic.R and stored at data/generated/ as traffic_counts_esri.rds

FHWA traffic intensity - https://geo.dot.gov/server/rest/services/Hosted/HPMS_FULL_NY_2019/FeatureServer/0 - downloaded at a_03_prep_traffic.R and stored at data/generated

Traffic co2 emissions at census block group from darte - DARTE_v2.gdb - downloaded from https://daac.ornl.gov/CMS/guides/CMS_DARTE_V2.html - stored at


#### geometry

sidewalks for nyc from Rhoads et al. (2021) A sustainable strategy for Open Streets in (post)pandemic cities  https://www.nature.com/articles/s42005-021-00688-z#data-availability - need to store edges shapefile in the data/raw/geometry folder

road centerline from nyc open data portal - downloaded from https://data.cityofnewyork.us/City-Government/NYC-Street-Centerline-CSCL-/exjm-f27b file name geo_export_a18d2619-e5c3-4105-bcf5-99ad8f35058e.shp - data/raw/geometry (1Gb)


#### motor vehicle collisions

motor vehicle collisions crashes NYC open data https://data.cityofnewyork.us/Public-Safety/Motor-Vehicle-Collisions-Crashes/h9gi-nx95 - Motor_Vehicle_Collisions_-_Crashes.csv - store at data/raw/traffic




