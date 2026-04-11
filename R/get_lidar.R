
library(dplyr)
library(httr)
library(sf)
library(raster)
library(ggplot2)
library(raster)
library(sf)
library(leaflet)
library(osmdata)
library(mapview)
library(htmlwidgets)
library(ggmap)
library(viridis)
library(stringr)


# using the package gblidar import from EAs lidar catalouge
get_gb_lidar = function(domain,product = c("DSM", "DTM"),product_type = c("elevation", "hillshade"), res = c(1,2)){
  
  #> Linking to GEOS 3.12.1, GDAL 3.8.3, PROJ 9.3.1; sf_use_s2() is TRUE
  if (rlang::is_installed("terra")) {
    library(terra)
    options(gblidar.out_raster_type = "SpatRaster")
  }
  #> terra 1.7.71
  
  options(gblidar.progress = FALSE) # for readability in this example.
  
  domain_os = domain |> 
    st_transform(27700) 
  
  # search_box <- st_point(c(532054, 181145)) |>
  #   st_buffer(1000) |>
  #   st_sfc() |>
  #   st_set_crs(27700)
  
  # the elevation data for the first return DSM
  r = eng_composite(domain_os, product = product)
  
  return(r)
  
}


# get nl lidar using the WCS service
get_nl_lidar = function(domain,product = c("DSM", "DTM"),res = c(0.5,5)){
  
  domain_rdnew = domain |> 
    st_transform(28992) |> 
    st_bbox()
  
  ahn = rAHNextract::ahn_area(
    bbox = domain_rdnew,
    resolution = res,
    dem = product
  )
  
  return(ahn)
  
}



## define an area
osm_bb <- st_bbox(c(xmin = 4.80, xmax = 4.9, ymax = 52.36, ymin = 52.32), crs = st_crs(4326))
mod_area <- st_as_sf(st_as_sfc(osm_bb))

# get osm land use using osmdata
get_osm_landuse = function(domain){
  
  osm_bb = st_bbox(domain)
  
  ##download landuse from osm
  x <- opq(bbox = osm_bb) %>%
    add_osm_feature(key = c('landuse')) %>%
    osmdata_sf()
  
  ##extract polygons
  landuse <- x$osm_polygons
  landuse <- select(landuse, osm_id, landuse)
  
}






# not sure about this
combine_rasters = function(rasters){
  
  rast_comb <- Reduce(function(x,y)mosaic(x,y, tolerance=1, fun=mean), lapply(rasters, raster))
  
}



# downloads NL lidar using the KBN tile codes
# if merge is false the tiles are written out
get_nl_lidar_tile = function(domain = NULL,point = NULL,product = c("DSM", "DTM"), merge = FALSE,
                             output_dir = "./"){
  
  # convert to dutch coords
  domain_rdnew = domain |> 
    st_transform(28992)
  
  # load the kaartblad grid
  KBN_grid = rAHNextract::ahn_sheets_info(AHN = "AHN4",dem = product,resolution = 0.5)
  
  grids_in = KBN_grid[domain_rdnew,]
  
  for (g in grids_in$kaartbladNr){
    tryCatch({
      
      if(product == "DSM"){   
        
        r = rast(paste0("https://service.pdok.nl/rws/ahn4/atom/downloads/dsm_05m/R_", g, ".tif"))
      }
      
      if(product == "DTM"){ 
        
        r = rast(paste0("https://service.pdok.nl/rws/ahn/atom/downloads/dtm_05m/M_", g, ".tif"))
        
      }   
      
      if(isTRUE(merge)){
        
        rast_comb <- Reduce(function(x,y)mosaic(x,y, tolerance=1, fun=mean), lapply(rasters, raster))
        
        return(rast_comb)
        
      } else {
        
        writeRaster(r,filename = paste0(output_dir,g,"_",product,".TIF"))
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
}

#r2 = get_gb_lidar(domain = domain, product = "DSM",res = 2)


#**The frontal area index problem**

#  The honest caveat is that λf properly requires **directional projection** of building faces — you want the area of building facades that are perpendicular to the wind for each sector. This is the hardest part. With a 1m nDSM (DSM minus DTM) you can approximate it by:

#  - Taking a transect through the wedge in the wind direction
#- Counting upward steps in height (leading edges of buildings)
#- Summing their heights

#It's not trivial but doable in R with `terra`. Alternatively the Grimmond & Oke method sidesteps this by using λp alone with an empirical correction, which is less accurate but much easier to compute.




library(gblidar)
library(sf)
library(raster)
library(terra)
library(mapview)
library(mapedit)
library(osmdata)
library(dplyr)
library(tidyr)
library(birk)
library(leaflet)
library(htmlwidgets)
library(tmap)

make_wedge = function(site_vect, direction_deg, half_angle = 11.25, 
                       radius = 1000) {
  # site_vect: terra SpatVector point (projected CRS, metres)
  # direction_deg: upwind direction to look (meteorological: direction wind comes FROM)
  # We want to look UPWIND, so add 180 to get the upwind direction
  
  #upwind_dir = (direction_deg- 180) %% 360
  
  upwind_dir <- direction_deg %% 360
  
  #upwind_dir = (90 - upwind_dir) %% 360
  
  # Convert to radians
  to_rad <- function(d) d * pi / 180
  
  site_vect = vect(site_vect)
  
  # Centre point coordinates
  cx <- crds(site_vect)[1, 1]
  cy <- crds(site_vect)[1, 2]
  
  # Generate arc points for the wedge boundary
  angles <- seq(
    to_rad(upwind_dir - half_angle),
    to_rad(upwind_dir + half_angle),
    length.out = 30
  )
  
  # Arc points at fetch radius
  arc_x <- cx + radius * sin(angles)
  arc_y <- cy + radius * cos(angles)
  
  # Close the polygon back through the site point
  coords <- rbind(
    c(cx, cy),
    cbind(arc_x, arc_y),
    c(cx, cy)
  )
  
  # Build sf polygon then convert to terra
  wedge_sf <- st_polygon(list(coords)) |>
    st_sfc(crs = crs(site_vect, proj = TRUE)) |>
    st_sf()
  
  wedge_vect <- vect(wedge_sf)
  return(wedge_vect)
}


options(gblidar.progress = TRUE) # for readability in this example.

pnt = st_point(c(532054, 181145)) |> 
  st_sfc() |>
  st_set_crs(27700) 



z0_wind_sectors_from_dsm = function(centre_point, buffer_m, n_sectors, country = c("ENG","NL")){
  
  if(country == "NL"){
    
    # create domain using default location in the function
    domain <- centre_point |>
      st_transform(28992) |> 
      st_buffer(buffer_m)
    
    dsm = get_nl_lidar(domain = domain,product = "DSM",res = 0.5)
    dtm = get_nl_lidar(domain = domain,product = "DTM",res = 0.5)
    
    dsm = dsm-dtm
    # shouldn't be any negative values, but is possible small discrepancies
    #dsm = dsm[dsm>=0]
    
  }
  
  if(country == "ENG"){
    
    # create domain using default location in the function
    domain <- centre_point |>
      st_transform(27700) |>
      st_buffer(buffer_m)
    
    if (rlang::is_installed("terra")) {
      library(terra)
      options(gblidar.out_raster_type = "SpatRaster")
    }
    
    # download DSM data for the domain
    dsm = get_gb_lidar(domain,res = 1, product = "dsm", product_type = "elevation")
    dtm = get_gb_lidar(domain,res = 1, product = "dtm", product_type = "elevation")
    
    dsm = dsm-dtm
    # shouldn't be any negative values, but is possible small discrepancies
    #dsm = dsm[dsm>=0]
    
  }
  
  # Define sectors
  sectors <- seq(0, 337.5, by=22.5)  # 16 sectors
  fetch <- buffer_m  # metres upwind
  
  results <- data.frame(sector=sectors, z0=NA, d=NA)
  wedges = list()
  for(i in seq_along(sectors)){
    
    direction <- sectors[i]
    
    # Create upwind wedge (±11.25 degrees, fetch distance)
    wedge <- make_wedge(centre_point, direction+11.25, half_angle=11.25, 
                        radius=fetch)  # custom function
    
    # Extract LiDAR heights within wedge
    heights <- terra::extract(dsm, wedge)
    
    # get rid of any NAs
    heights <- heights[!is.na(heights)]
    
    # Ground threshold - subtract DTM to get above-ground heights
    # (use DSM - DTM = nDSM/CHM)
    
    # Morphometric calculations
    H_mean <- mean(heights[heights > 2])  # >2m = buildings/trees
    lambda_p <- sum(heights > 2) / length(heights)
    
    # Macdonald z0 and d
    alpha <- 4.43
    d_h <- 1 + alpha^(-lambda_p) * (lambda_p - 1)
    d_val <- H_mean * d_h
    
    Cd = 1.2; k = 0.4; beta = 1.0
    lambda_f = lambda_p * 0.5  # rough approximation without directional projection
    
    z0_val <- H_mean * (1 - d_h) * exp(
      -(0.5 * beta * (Cd/k^2) * (1-d_h) * lambda_f)^(-0.5))
    
    wedge_out = wedge |> 
      st_as_sf() |> 
      mutate(sector = direction,
             z0 = z0_val,
             d = d_val)
    
    #results$z0[i] <- z0_val
    #results$d[i] <- d_val
    wedges[[i]] = wedge_out
  }
  
  wedge_sf = do.call(rbind,wedges)
  
  return(wedge_sf)
  
}

z0_from_dsm = function(domain, country = c("ENG","NL")){

if(country == "NL"){
  
  # create domain using default location in the function
  domain_nl <- domain |>
    st_transform(28992) 
  
  dsm = get_nl_lidar(domain = domain_nl,product = "DSM",res = 0.5)
  dtm = get_nl_lidar(domain = domain_nl,product = "DTM",res = 0.5)
  
  dsm = dsm-dtm
  # shouldn't be any negative values, but is possible small discrepancies
  #dsm = dsm[dsm>=0]
  
  heights <- terra::extract(dsm, domain_nl)
  
}

if(country == "ENG"){
  
  # create domain using default location in the function
  domain_eng <- domain |>
    st_transform(27700)
  
  if (rlang::is_installed("terra")) {
    library(terra)
    options(gblidar.out_raster_type = "SpatRaster")
  }
  
  print("downloading DSM")
  
  # download DSM data for the domain
  dsm = get_gb_lidar(domain_eng,res = 1, product = "dsm", product_type = "elevation")
  print("downloading DTM")
  dtm = get_gb_lidar(domain_eng,res = 1, product = "dtm", product_type = "elevation")
  print("subtracting DTM from DSM to get absolute heights")
  dsm = dsm-dtm
  # shouldn't be any negative values, but is possible small discrepancies
  #dsm = dsm[dsm>=0]
  
  # Extract LiDAR heights within wedge
  heights <- terra::extract(dsm, domain_eng)
  
}


  # get rid of any NAs
  heights <- heights[!is.na(heights)]
  
  # Ground threshold - subtract DTM to get above-ground heights
  # (use DSM - DTM = nDSM/CHM)
  
  # Morphometric calculations
  H_mean <- mean(heights[heights > 2])  # >2m = buildings/trees
  lambda_p <- sum(heights > 2) / length(heights)
  print("running the calcs")
  # Macdonald z0 and d
  alpha <- 4.43
  d_h <- 1 + alpha^(-lambda_p) * (lambda_p - 1)
  d_val <- H_mean * d_h
  
  Cd = 1.2; k = 0.4; beta = 1.0
  lambda_f = lambda_p * 0.5  # rough approximation without directional projection
  
  z0_val <- H_mean * (1 - d_h) * exp(
    -(0.5 * beta * (Cd/k^2) * (1-d_h) * lambda_f)^(-0.5))
  
  return(z0_val)
  
}


