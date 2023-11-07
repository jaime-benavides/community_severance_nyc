
# kriging from mid point 
# adapted from Criado et al. (2022) https://earth.bsc.es/gitlab/es/universalkriging/-/blob/production/general/UK_mean.R
regrid_ok <- function(non_uniform_data, target_grid,crs_sim = "+proj=utm +zone=31 +ellps=intl +units=m +no_defs"){
  if(isTRUE(class(non_uniform_data) != "SpatialPointsDataFrame")){
    non_uniform_data <- sp::SpatialPointsDataFrame(non_uniform_data[,c("X", "Y")],non_uniform_data)
    sp::proj4string(non_uniform_data) <- crs_sim
  }
  vf_ok      <- automap::autofitVariogram(aadt ~ 1, non_uniform_data)
  ok_regular <- gstat(formula = aadt ~ 1, data = non_uniform_data, model = vf_ok$var_model, nmax = 20) 
  regular <- predict(ok_regular, target_grid)
  regular_sf <- sf::st_as_sf(regular)
  # regular <- sp::spTransform(regular, sp::CRS("+proj=longlat"))
  # pixels <- sp::SpatialPixelsDataFrame(regular,tolerance = 0.99, as.data.frame(regular[,"var1.pred"]))
  # mean_raster <- raster::raster(pixels[,'var1.pred'])
  # return(mean_raster)}
  return(regular_sf)}



print_patterns_loc <- function (pats, colgroups = NULL, n = 1:6, pat_type = "pat", 
                                title = "", size_line = 1, size_point = 1) 
{
  if (!is.null(colgroups)) {
    colgroups <- colgroups %>% dplyr::rename(chem = !!names(colgroups)[1])
  }
  else {
    colgroups <- data.frame(chem = rownames(pats), group = "1")
  }
  if (n > ncol(pats)) 
    n <- ncol(pats)
  grouping <- names(colgroups)[2]
  colnames(pats) <- paste0(pat_type, stringr::str_pad(1:ncol(pats), 
                                                      width = 2, pad = "0", side = "left"))
  pats.df <- pats %>% tibble::as_tibble() %>% dplyr::mutate(chem = colgroups[[1]]) %>% 
    tidyr::pivot_longer(-chem, names_to = "pattern", values_to = "loading") %>% 
    dplyr::right_join(., colgroups, by = "chem")
  pats.df$chem <- factor(as.character(pats.df$chem), levels = unique(as.character(pats.df$chem)))
  loadings <- pats.df %>% dplyr::filter(pattern %in% paste0(pat_type, 
                                                            stringr::str_pad(n, width = 2, pad = "0", side = "left"))) %>% 
    ggplot(aes(x = chem, y = loading, color = !!sym(grouping))) + 
    geom_point(size = size_point) + geom_segment(aes(yend = 0, xend = chem), size = size_line) + 
    facet_wrap(~pattern) + theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size=12), legend.title = element_text(size=14),
                                              axis.text.x = element_text(angle = 45, hjust = 1, size = 14), strip.background = element_rect(fill = "white"), 
                                              axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    geom_hline(yintercept = 0, size = 0.2) #+ ggtitle(title)
  loadings
}


format_pred_mod <- function(pred.csi) {
  fit.table.csi <- as.data.frame(pred.csi$matRRfit)  
  colnames(fit.table.csi) <- paste0("rr_", colnames(fit.table.csi))
  fit.table.csi <-fit.table.csi %>% mutate(csi = as.numeric(row.names(fit.table.csi)))
  # 3d.ii Extract 95% CI  
  lci.table.csi <- as.data.frame(pred.csi$matRRlow)  
  colnames(lci.table.csi) <- paste0("lci_", colnames(lci.table.csi))
  uci.table.csi <- as.data.frame(pred.csi$matRRhigh)  
  colnames(uci.table.csi) <- paste0("uci_", colnames(uci.table.csi))
  ## plot 
  plot.csi_mod <- data.frame(fit.table.csi, lci.table.csi, uci.table.csi)
  names(plot.csi_mod) <- c("pe", "csi", "lci", "uci")
  return(plot.csi_mod)
}

plot_pred_mod <- function(mod_dta) {
  mynamestheme <- theme(
    plot.title = element_text(family = "Helvetica", hjust = 0.5, size = (18)),
    axis.title = element_text(family = "Helvetica", size = (25)),
    axis.text = element_text(family = "Helvetica", size = (18))
  )
  
  p <- ggplot() +
    geom_rug(aes(x = csi),
             data = mod_dta,
             sides = "b", length = grid::unit(0.02, "npc")) +
    geom_ribbon(aes(ymin = lci, ymax = uci, x = csi), fill = "blue",
                alpha = 0.2, data = mod_dta) + 
    geom_line(aes(x = csi, y = pe), lwd = 1, data = mod_dta) +
    geom_hline(yintercept=1, linetype="dashed", color = "grey", size = 1.5) + 
    labs(y = "Motor vehicle collisions", x = "Community severance index") + 
    theme_bw() +
    mynamestheme
  return(p)
}


estim_areal_road_ped_ratio  <- function(sp_unit_pos) {
  sp_unit <- grid_contxt[sp_unit_pos[1],]
  influence_area <- sf::st_buffer(sp_unit, dist = 804.672)
  roads_buffer_id <- sapply(sf::st_intersects(roads_contxt, influence_area),function(x){length(x)>0})
  roads_local <- roads_contxt[roads_buffer_id, ]
  sidewalk_buffer_id <- sapply(sf::st_intersects(edges_nyc, influence_area),function(x){length(x)>0})
  sidewalk_local <- edges_nyc[sidewalk_buffer_id, ]
  if(nrow(roads_local) > 0){
    roads_local_area <- sum(roads_local$st_area, na.rm = T)
    bikes_local_area <- sum(roads_local$bike_lane_area, na.rm = T)
    sidewalk_local_area <- sum(sidewalk_local$sidewalk_area_m, na.rm = T)
    roads_local_area <- roads_local_area - bikes_local_area
    streets_local_area <- roads_local_area + bikes_local_area + sidewalk_local_area
    sidewalk_to_road <- sidewalk_local_area / roads_local_area
    sidewalk_to_street <- sidewalk_local_area / streets_local_area
    road_to_street <- roads_local_area / streets_local_area
  } else {
    streets_local_area <- NA
    roads_local_area <- NA
    bikes_local_area <- NA
    sidewalk_local_area <- NA
    sidewalk_to_road <- NA
    sidewalk_to_street <- NA
    road_to_street <- NA
  }
  areal_ratio <- c(sp_unit_pos[2], streets_local_area, roads_local_area, bikes_local_area, sidewalk_local_area, sidewalk_to_road, sidewalk_to_street, road_to_street) 
  return(areal_ratio)
}


 


