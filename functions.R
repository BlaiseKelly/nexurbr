library(dplyr)
library(ncdf4)
library(stringr)
library(birk)
library(openair)
library(lubridate)
library(purrr)
library(sf)
library(raster)
# # Install saqtrendr
# remotes::install_github("skgrange/saqtrendr")
# # Development version
# remotes::install_github("skgrange/threadr")
# # Install sspatialr
# remotes::install_github("skgrange/sspatialr")
library(saqtrendr)
library(threadr)
library(sspatialr)

# this function writes openair plots to png files
openair2png <- function(plot, path, plot_name, width = 10500, height = 10000, units="px", res = 800){

  filename <- paste0(path, plot_name,".png")
  png(filename, width=width, height=height, units=units, res=res)
  print(plot)
  dev.off()

}


# this function imports the historic nc files from ICOS and makes a single dataframe for all the species and heights which
# matches the original format that most of the scripts work off
combine_historic_nc <- function(dir_location,
                                species,
                                heights){

  mol_conv <- data.frame(species = c("co2", "ch4", "n2o"),
                         mol_pp = c(10^6,10^9,10^9),
                         final_unit = c("ppm", "ppb", "ppb"))

  all_dat <- list()
  for (s in species){

    # find nc files with species that match
    filez <- list.files(scratch, s, full.names = TRUE)

    conv <- filter(mol_conv, species == s)

    # loop through heights
    for (h in heights){

      # pick out the species file for the specific height
      filez_h <- filez[grepl(paste0("-", h), filez)]

      # open the nc file
      icos <- nc_open(paste0(filez_h))

      # import time dimension
      TIME <- ncvar_get(icos, "time")

      # for some reason the ICOS data is at half hourly time steps, with the first 01:30. Assume it should be 01:00
      TIME_corrected <- TIME-3600/2

      # find the base date
      time_since <- gsub("T", " ", str_sub(icos$dim$time$units, 15, -2))

      # create formatted times
      d8_time <- lubridate::ymd_hms(time_since) + lubridate::seconds(TIME_corrected)

      # create a data frame of the dates. For some reason the dates are on the half past, so subtract 1800 seconds
      df <- data.frame(date = d8_time)

      ##extract variables and convert to ppm (divide by 10^6)
      df$value <- ncvar_get(icos, "value")*conv$mol_pp

      # assign the species name and height
      names(df) <- c("date", paste0(toupper(s), "_", h))

      # write to list
      all_dat[[paste0(s,h)]] <- df

      # let the user know where up to
      print(paste("imported",s,h))
    }

  }

  # left join them all by date to create one dataframe
  big_df <- purrr::reduce(all_dat, dplyr::left_join, by = c('date'))

  return(big_df)

}
# dat_dir = "DAT/"
# # this function extracts a variable from a ecmwf meteo file for specific coordinates
# import_ecmwf <- function(dat_dir,
#                          lat,
#                          lon){
# 
#   dat_dirz <- list.files(dat_dir,
#                          pattern = ".nc", full.names = TRUE)
# 
#   all_years <- list()
# 
#   for (d in dat_dirz){
# 
#     ecmwf <- nc_open(d)
# 
#     species = names(ecmwf$var)[3]
# 
#     TIME <- ncvar_get(ecmwf, "valid_time")
# 
#     time_since <- str_sub(ecmwf$dim$valid_time$units, 15, -1)
# 
#     d8_time <- lubridate::ymd(time_since) + lubridate::seconds(TIME)
# 
#     met_lon <- ncvar_get(ecmwf, 'longitude')
#     met_lat <- ncvar_get(ecmwf, 'latitude')
# 
#     lon_ind <- which.closest(met_lon, lon)
#     lat_ind <- which.closest(met_lat, lat)
# 
#     ##extract no2 variables for entire domain for 1 time step and altitude
#     b_met <- data.frame(date = d8_time,
#                         nam = ncvar_get(ecmwf, species, start = c(lon_ind,lat_ind,1), count = c(1,1,NROW(d8_time))))
# 
#     names(b_met) = c("date", species)
# 
#     all_years[[species]] <- b_met
#     print(d)
#   }
# 
#   param_df <- reduce(all_years, left_join)
# 
#   row.names(param_df) <- NULL
# 
#   param_df = openair::cutData(param_df,type = "daylight", latitude = lat, longitude = lon)
# 
#   return(param_df)
# 
# }

snap_to_era5_grid = function(domain, buffer_m = 10000, grid_deg = 0.25) {
  
  # Get bounding box in WGS84
  bbox = domain |> 
    st_transform(4326) |> 
    st_bbox()
  
  # Snap outward to nearest 0.25 degree grid
  north = ceiling(bbox["ymax"] / grid_deg) * grid_deg
  south = floor(bbox["ymin"] / grid_deg) * grid_deg
  east  = ceiling(bbox["xmax"] / grid_deg) * grid_deg
  west  = floor(bbox["xmin"] / grid_deg) * grid_deg
  
  # Check minimum span — must cover at least 2 grid cells in each direction
  # otherwise CDS sometimes returns nothing
  if ((north - south) < 2 * grid_deg) {
    north = north + grid_deg
    south = south - grid_deg
  }
  if ((east - west) < 2 * grid_deg) {
    east = east + grid_deg
    west = west - grid_deg
  }
  
  paste0(north, "/", west, "/", south, "/", east)
}



get_era_aermod = function(domain,buffer_m = 0, folder_name, start_date,end_date, folder_name){
  
 
# find the coords and snap to grids so it downloads properly
cds_area = snap_to_era5_grid(domain)

path <- "./"

d8s = unique(str_sub(as.character(seq(start_date,end_date)),1,-4))

for (d in d8s){

yrz <- str_sub(d,1,-4)

mth = as.numeric(str_sub(d,6,-1))

# get the variables from ecmwf for download
variables_single <- c("friction_velocity",
                      "boundary_layer_height",
                      "forecast_surface_roughness",
                      "total_cloud_cover")

vzs <- get_era(bb = ecmwf_format_bb,location = folder_name,months = mth,
               variables = variables_single, dataset = "single-level", api_key = "c95fc5f7-4e83-4e08-88ff-f02d03baa66e", path_out = "DAT/meteo/", years = yrz)


# get the variables from ecmwf for download
variables_land <- c("2m_dewpoint_temperature",
                    "2m_temperature",
                    "10m_u_component_of_wind",
                    "10m_v_component_of_wind",
                    "surface_pressure",
                    "surface_sensible_heat_flux",
                    "total_precipitation")

get_era(bb = ecmwf_format_bb,location = folder_name,months = mth,
        variables = variables_land, dataset = "land", api_key = "c95fc5f7-4e83-4e08-88ff-f02d03baa66e", path_out = "DAT/meteo/", years = yrz)

}
}

get_worldmet_aermod = function(){}

create_aermod_meteo = function(domain_centre, country_z0){
  
  z0_sectors = z0_wind_sectors_from_dsm(centre_point = domain_centre,buffer_m = 500,country = country_z0) |> 
    mutate(sector = as.character(sector)) |> 
    st_set_geometry(NULL)
  
  coords2get = st_coordinates(st_transform(domain_centre,4326))
  
  met_df = import_ecmwf(dat_dir = paste0("DAT/meteo/",Project_name,"/"), lat = coords2get[2],lon = coords2get[1])
  
  # height of wind measurements
  z_anem = 10
  
  met = met_df |>
    arrange(date) |>
    mutate(
      sshf_hourly = sshf - lag(sshf, default = 0),
      # Reset at midnight - when lag crosses day boundary, use value as-is
      new_day = as.Date(date) != as.Date(lag(date, default = date[1])),
      sshf_hourly = ifelse(new_day, sshf, sshf_hourly),
      sshf = sshf_hourly/3600,
      sshf_aermod = -sshf,
      tp = ifelse(new_day, tp, tp - lag(tp, default = 0)
      )) |> 
    mutate(yr = year(date),
           m = month(date),
           md = mday(date),
           jd = yday(date),
           h = hour(date),
           pres = sp/100,
           temp = t2m,
           pamt = tp*1000,
           ccvr = tcc*10,
           ws = sqrt(u10^2 + v10^2),
           wd = (270 - atan2(-v10, -u10) * 180/pi) %% 360, # wind direction from
           es_T =  611.2 * exp(17.67 * (t2m - 273.15) / (t2m - 29.65)),
           es_d = 611.2 * exp(17.67 * (d2m - 273.15) / (d2m - 29.65))) |>
    mutate(wd_bin = cut(wd, breaks = seq(0, 360, by = 22.5),
                        labels = seq(0, 337.5, by = 22.5),
                        include.lowest = FALSE, right = FALSE)) |> 
    mutate(wd_bin = as.character(wd_bin)) |> 
    left_join(z0_sectors, by = c("wd_bin" = "sector")) |> 
    mutate(rh = 100 * (es_d / es_T),
           rho = sp/(287.05*t2m),
           rho_cp = rho * 1005.0,
           ustar = (ws * 0.4) / log((z_anem - d) / z0)) |> 
    mutate(sshf_aermod = ifelse(
      daylight == "nighttime" & sshf_aermod > 0,
      -sshf_aermod,  # or just a small negative value like -5
      sshf_aermod
    )) |> 
    mutate(
      # L = ifelse(abs(sshf_aermod) < 2,
      #                 ifelse(sshf_aermod >= 0, 99999, -99999),
      #                 -(rho_cp * t2m * zust^3) / (0.4 * 9.80665 * sshf_aermod)),
      L = -(rho_cp * t2m * ustar^3) / (0.4 * 9.80665 * sshf_aermod),
      wstar = ifelse(sshf_aermod > 0 & daylight == "daylight",
                     (9.80665/t2m * (sshf_aermod/rho_cp) * blh)^(1/3),
                     -999),
      VPTG = (0.011 * t2m) / 9.80665,
      Zic = ifelse(daylight == "nighttime", -999, blh),  # convective
      Zim = ifelse(daylight == "nighttime", blh, -999)) |>
    mutate(wstar = ifelse(daylight == "nighttime" | sshf_aermod < 0 , -999, wstar)) |>
    mutate(sigma_theta = ifelse(sshf_aermod > 0,
                                (0.7 * wstar + 1.5 * ustar) / ws,  # unstable
                                1.5 * ustar / ws)) |>                    # stable/neutral, in radians -> degrees
    mutate(Wdnn = sigma_theta * (180/pi)) |> 
    mutate(Wsnn = case_when(
      L < -10   ~ ustar * 1.6 * (1 - 10/L)^(1/3),
      L > 10    ~ ustar * 1.3 * (1 + 2.0 * 10/L)^(-1/3),
      TRUE      ~ 1.3 * ustar  # near-neutral band
    )) 
  
  return(met)
}


write_meteo = function(meteo_dat){
  
  meteo_df = meteo_dat
  
  ##format for outputs
  SFC <- data.frame(YEAR = str_sub(met$yr,3,-1),
                    MNTH = met$m,
                    MDAY = met$md,
                    JDAY = met$jd,
                    HOUR = met$h+1,
                    HFLX = gsub("[[:space:]]", "", format(round(met$sshf, 1), nsmall = 1)), ## sensible heat flux (W/m2)
                    USTR = gsub("[[:space:]]", "", format(round(met$ustar, 3), nsmall = 3)), ## surface friction velocity (m/s)
                    COVS = gsub("[[:space:]]", "", format(round(met$wstar, 3), nsmall = 3)), ## convective velocity scale (m/s)
                    UALR = gsub("[[:space:]]", "", format(round(met$VPTG, 3), nsmall = 3)), ## lapse rate above mixing height (K/m)
                    COMH = gsub("[[:space:]]", "", format(round(met$Zic, 4), nsmall = 4)), ## convective mixing height (m)
                    MEMH = gsub("[[:space:]]", "", format(round(met$Zim, 4), nsmall = 4)), ## mechanical mixing height (m)
                    MONI = gsub("[[:space:]]", "", format(round(met$L, 4), nsmall = 4)), ## ## Monin-Obukhov length (m)
                    ZOHT = gsub("[[:space:]]", "", format(round(model_df$z0, 2), nsmall = 2)), ## surface roughness length (m)
                    BORA = gsub("[[:space:]]", "", format(round(met$B0, 2), nsmall = 2)), ## bowen ratio
                    ALBE = gsub("[[:space:]]", "", format(round(0, 2), nsmall = 2)), ## albedo
                    WSPD = gsub("[[:space:]]", "", format(round(met$ws, 2), nsmall = 2)), ## wind speed (m/s)
                    WDIR = gsub("[[:space:]]", "", format(round(met$wd, 1), nsmall = 1)), ## wind direction (degrees)
                    REWH = gsub("[[:space:]]", "", format(round(met$zref, 1), nsmall = 1)), ## reference wind height (m)
                    TTnn = gsub("[[:space:]]", "", format(round(met$temp, 1), nsmall = 1)), ## Ambient temperature (K)
                    RETH = gsub("[[:space:]]", "", format(round(met$ztemp, 1), nsmall = 1)), ## reference temperature height (m)
                    IPCO = gsub("[[:space:]]", "", round(met$ipcode)), ## preceipitation code (also called ipcode)
                    PAMT = gsub("[[:space:]]", "", round(met$pamt)), ## precipitation (mm)
                    RHUM = gsub("[[:space:]]", "", round(met$rh)), ## relative humidity (%)
                    PRES = gsub("[[:space:]]", "", round(met$pres)), ## surface pressure (mb)
                    TSKC = gsub("[[:space:]]", "", round(met$ccvr)), ## sky cover in tenths
                    WSAD = "NAD-OS", ## wind speed adjustment
                    BLAH = "NoSubs") ## Not sure
  
  SFC <- arrange(SFC, YEAR, MNTH, MDAY, HOUR)
  
  PFL <- data.frame(YEAR = met$yr,
                    MNTH = met$m,
                    MDAY = met$md,
                    HOUR = met$h+1,
                    HEIG = format(round(10.0, 1), nsmall = 1),
                    TFLG = round(1),
                    WDIR = gsub("[[:space:]]", "", format(round(met$wd, 1), nsmall = 1)), ## wind direction (degrees)
                    WSPD = gsub("[[:space:]]", "", format(round(met$ws, 2), nsmall = 2)), ## wind speed (m/s)
                    TTnn = gsub("[[:space:]]", "", format(round(met$temp, 2), nsmall = 2)), ## Ambient temperature (K)
                    SANN = gsub("[[:space:]]", "", format(round(met$Wdnn, 2), nsmall = 2)), ## standard deviation of wind direction (degrees)
                    SWNN = gsub("[[:space:]]", "", format(round(met$Wsnn, 2), nsmall = 2)) ## standard deviation of wind speed (m/s)
  )
  
  PFL <- arrange(PFL, YEAR, MNTH, MDAY, HOUR)
  
  ##write the meteo outputs
  ##specify output file name
  out = paste0("meteo/", Project,".SFC")
  
  #write a RECEPTOR input file (not currently used in the model, but an option)
  write.table(SFC,
              file = out, col.names = FALSE, row.names = FALSE,
              sep = " ", quote = FALSE)
  
  #add the header lines to the created RECEPTOR file
  fConn <- file(out, "r+")
  Lines <- readLines(fConn)
  
  
  NS <- ifelse(lon>0, "N", "S")
  EW <- ifelse(lat>0, "E", "W")
  
  header <- paste0("    ", abs(round(lon, 2)), NS, "    ", abs(round(lat, 2)), EW, "          ", "UA_ID:", " 00014735  SF_ID:    14737  OS_ID:   000001     VERSION: 18081")
  
  writeLines(c(
    header,
    Lines)
    , con = fConn)
  close(fConn)
  
  ##Profile file (.PFL)
  out = paste0("meteo/", Project,".PFL")
  
  #write a RECEPTOR input file (not currently used in the model, but an option)
  write.table(PFL,
              file = out, col.names = FALSE, row.names = FALSE,
              sep = " ", quote = FALSE)
  
  
}

#' Convert LE Meteorological NetCDF Data to Raster Bricks
#'
#' This function processes meteorological data from LOTOS-EUROS and converts it into raster bricks.
#'
#' @param date_folder Date folder for processing.
#' @param get_forecasts Logical, whether to include forecast data.
#' @param region The geographical region of interest.
#' @param LE_output_folder Path to LOTOS-EUROS output directory.
#' @param brick_outputs Boolean, whether to save outputs as bricks.
#' @param variables Vector of variable names (e.g., c("temper", "rain")).
#' @param output Type of output required ("hourly", "daily").
#' @param output_dir Directory where output files will be saved.
#' @return None (writes raster files to disk).
#' @export
LE_meteo2brick <- function(date_folder, get_forecasts, region, LE_output_folder, brick_outputs, variables, output, output_dir) {
  # Function implementation here


  ## test variables
  # date_to_get = tday
  # get_forecasts = FALSE
  # LE_output_folder = LEOPARD
  # variables = c("temper", "rain", "usurf", "vsurf")
  # output = c("hourly", "daily")
  # output_dir = output_folder


  if(get_forecasts == TRUE){
    ## pick days to include - yesterday, today and tomorrow
    d8s_in <- as.character(c(date_folder-1, date_folder, date_folder+1))
  } else{
    d8s_in <- as.character(date_folder)
  }

  ## loop through each day and output bricks
  for (d8 in d8s_in){

    #day_number <- paste0( "D", which(d8 == d8s_in))


    d8 <- ymd(d8)
    if(get_forecasts == TRUE){
      date2get <- gsub("-", "", as.character(d8))
    } else {
    date2get <- gsub("-", "", as.character(d8-1))
    }

    # how many hours are in the nc file, only 6 for the last day
    if(which(d8 == d8s_in) < 3){
      n_hrs = 24
    } else {
      n_hrs = 6
    }

    # want to get the future days model as it uses historic meteo rather than the forecast
    y <- year(d8)
    m <- sprintf("%02d", month(d8))
    d_1 <- sprintf("%02d", day(d8))


    ## READ IN METEO
    ## create dir name
    le_files <- paste0(LE_output_folder, y, "/", m, "/", d_1, "/lotos-euros/C71-oper-Sectors-",region, "/output/LE_C71-oper-Sectors_meteo_", date2get,".nc")
    ## open the nc file
    lemet <- nc_open(le_files)
    ## get the lat and lon centroid data
    lon_ext <- ncvar_get(lemet, "longitude_bnds")
    lat_ext <- ncvar_get(lemet, "latitude_bnds")

    ## find the extent to generate a domain
    x_min <- min(lon_ext)
    x_max <- max(lon_ext)
    y_min <- min(lat_ext)
    y_max <- max(lat_ext)

    ##create a data frame to setup polygon generation
    df <- data.frame(X = c(x_min, x_max, x_max, x_min),
                     Y = c(y_max, y_max, y_min, y_min))

    ##generate a polygon of the area
    domain_ll <- df %>%
      st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
      dplyr::summarise(data = st_combine(geometry)) %>%
      st_cast("POLYGON")


    # varz <- data.frame(le_vars = c("temper", "rain", "usurf", "vsurf"),
    #                    dims = c(4,3, 3,3))

    for (v in variables){

      #var_df <- filter(varz, le_vars == v)

      tryCatch({

        lons <- ncvar_get(lemet, "longitude")
        lats <- ncvar_get(lemet, "latitude")

        nlon <- NROW(lons)
        nlat <- NROW(lats)

        var_dim <- NROW(dim(ncvar_get(lemet, v)))

        if(var_dim == 4){
          var_in <- ncvar_get(lemet, v, start = c(1,1,1,1), count = c(nlon, nlat, 1, n_hrs))
        }

        if(var_dim == 3){
          var_in <- ncvar_get(lemet, v, start = c(1,1,1), count = c(nlon, nlat, n_hrs))
        }


        TIME <- ncvar_get(lemet, "time", start = c(1), count = c(n_hrs))

        time_since <- str_sub(lemet$dim$time$units, 15, -14)

        d8_time <- lubridate::ymd(time_since) + lubridate::seconds(TIME)
        d8_char <- format(d8_time, "%Y-%m-%d %H:%M:%S")
        d8_epoch <- as.integer(d8_time)

        ##generate raster brick from it
        var_met <- brick(var_in[1:nlon,1:nlat,1:24])

        r_met <- t(var_met)
        ##define the extent
        crs(r_met) <- 4326
        bb <- extent(domain_ll)
        extent(r_met) <- bb
        r1 <- flip(r_met, direction = 'y')

        # write out bricks
        if(any(output == "daily")){

          dir.create(paste0(output_dir, "meteo/", r, "/daily/"), recursive = TRUE)

          r_mean <- calc(r1,mean)

          terra::writeRaster(r_mean, filename=paste0(output_dir, "meteo/", region, "/daily/",v, "_", date2get,".TIF"), overwrite = TRUE)

        }

        # write out bricks
        if(any(output == "hourly")){

          dir.create(paste0(output_dir, "meteo/", r, "/hourly/"), recursive = TRUE)

          names(r1) <- d8_char

          terra::writeRaster(r1, filename=paste0(output_dir,"meteo/", region, "/hourly/", v, "_", date2get, ".TIF"), overwrite = TRUE)

        }

        # show where up to
        print(paste0(v, " ", date2get))


      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    }

  }

}


#' This function processes concentration data from LOTOS-EUROS labelled netcdfs and converts it into raster bricks.
#'
#' @param date_folder Date folder for processing.
#' @param get_forecasts Logical, whether to include forecast data.
#' @param region The geographical region of interest.
#' @param LE_output_folder Path to LOTOS-EUROS output directory.
#' @param brick_outputs Boolean, whether to save outputs as bricks.
#' @param variables Vector of variable names (e.g., c("no2", "o3")).
#' @param output Type of output required ("hourly", "daily").
#' @param output_dir Directory where output files will be saved.
#' @return None (writes raster files to disk).
#' @export
LE_lbl_concs2brick <- function(date_folder,
                           get_forecasts,
                           region,
                           LE_output_folder,
                           brick_outputs,
                           variables,
                           output,
                           output_dir){


  # date_folder = Sys.Date()
  # get_forecasts = TRUE
  # region = "NL"
  # LE_output_folder = "//tsn.tno.nl/RA-Data/Express/ra_express_modasuscratch_unix/projects/EU/CAMS/C71/Sectors/LEOPARD/"
  # variables = c("no2")
  # output = c("hourly", "daily")
  # output_dir = "//tsn.tno.nl/RA-Data/Express/ra_express_modasuscratch_unix/projects/KIP/2025/EZK_dashboard/outputs/"

  ## look up table of species names for LE, the factor convert from mol mol to ug/m3 and the title for tmap
  varient_df <- data.frame(le_vars = c("o3", "no2", "no", "nh3", "so2", 'hno3', "co", "tpm25", "tpm10","ec_f", "ec_c"),
                           omrekenfactor = c(1e12*0.048/24.4,  1e12*0.046005/24.4, 1e12*0.03001/24.4, 1e12*0.017031/24.4, 1e12*0.0640628/24.4, 1e12*0.06301/24.4,1e12*0.02801/24.4, 1e9, 1e9,1e9, 1e9)
                           )


  tryCatch({



  if(get_forecasts == TRUE){
    ## pick days to include - yesterday, today and tomorrow
    d8s_in <- as.character(c(date_folder-1, date_folder, date_folder+1))
  } else{
    d8s_in <- as.character(date_folder)
  }

    # want to get the future days model as it uses historic meteo rather than the forecast
    y <- year(date_folder)
    m <- sprintf("%02d", month(date_folder))
    d_1 <- sprintf("%02d", day(date_folder))

  ## loop through each day and output bricks
  for (d8 in d8s_in){

    #day_number <- paste0( "D", which(d8 == d8s_in))

    d8 <- ymd(d8)


    d8 <- ymd(d8)
    if(get_forecasts == TRUE){
      date2get <- gsub("-", "", as.character(d8))
    } else {
      date2get <- gsub("-", "", as.character(d8-1))
    }

    # how many hours are in the nc file, only 6 for the last day
    if(which(d8 == d8s_in) < 3){
      n_hrs = 24
    } else {
      n_hrs = 6
    }


    ## READ IN METEO
    le_files <- paste0(LE_output_folder, y, "/", m, "/", d_1, "/lotos-euros/C71-oper-Sectors-",region,"/output/LE_C71-oper-Sectors_labelled-conc-sfc_", date2get,".nc")

    lelab <- nc_open(le_files)

    lon_ext <- ncvar_get(lelab, "longitude_bnds")
    lat_ext <- ncvar_get(lelab, "latitude_bnds")

    ## find the extent to generate a domain
    x_min <- min(lon_ext)
    x_max <- max(lon_ext)
    y_min <- min(lat_ext)
    y_max <- max(lat_ext)

    ##create a data frame to setup polygon generation
    df <- data.frame(X = c(x_min, x_max, x_max, x_min),
                     Y = c(y_max, y_max, y_min, y_min))

    ##generate a polygon of the area
    domain_ll <- df %>%
      st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
      dplyr::summarise(data = st_combine(geometry)) %>%
      st_cast("POLYGON")


    lab_names <- ncvar_get(lelab, "labelnames")

    for (v in variables){

      tryCatch({

        var_info <- filter(varient_df, le_vars == v)

        lons <- ncvar_get(lelab, "longitude")
        lats <- ncvar_get(lelab, "latitude")

        nlon <- NROW(lons)
        nlat <- NROW(lats)

        var_in <- ncvar_get(lelab, v, start = c(1,1,1,1,1), count = c(nlon, nlat, 1, n_hrs, 18))

        TIME <- ncvar_get(lelab, "time", start = c(1), count = c(n_hrs))

        time_since <- str_sub(lelab$dim$time$units, 15, -14)

        d8_time <- lubridate::ymd(time_since) + lubridate::seconds(TIME)
        d8_char <- format(d8_time, "%Y-%m-%d %H:%M:%S")
        d8_epoch <- as.integer(d8_time)
        l <- lab_names[1]
        for (l in lab_names){

          lab_num <- which(l == lab_names)

          tryCatch({
            ##generate raster from it

            r_rd <- t(brick(var_in[1:nlon,1:nlat, 1:n_hrs, lab_num]))

            ##define the extent
            crs(r_rd) <- 4326
            bb <- extent(domain_ll)
            extent(r_rd) <- bb
            r1 <- flip(r_rd, direction = 'y')
            r1_ug <- r1*var_info$omrekenfactor
            names(r1_ug) <- d8_char

            nam_lab <- paste0(v, "_", date2get, "_", sprintf("%04d", lab_num))

            # write out bricks
            if(any(output == "daily")){

              dir.create(paste0(output_dir, "concs/", region, "/daily/"), recursive = TRUE)

              r_mean <- calc(r1,mean)

              terra::writeRaster(r_mean, filename=paste0(output_dir, "concs/", region, "/daily/",nam_lab, ".TIF"), overwrite = TRUE)

            }

            # write out bricks
            if(any(output == "hourly")){

              dir.create(paste0(output_dir, "concs/", region, "/hourly/"), recursive = TRUE)

              names(r1) <- d8_char

              terra::writeRaster(r1, filename=paste0(output_dir,"concs/", region, "/hourly/", nam_lab, ".TIF"), overwrite = TRUE)

            }

            # show where up to
            print(paste0("reading from folder ", date_folder, " and downloading ", nam_lab))

          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})



        }

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }

  }





      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

#' This function processes concentration data from non labelled LOTOS-EUROS and converts it into raster bricks.
#'
#' @param date_folder Date folder for processing.
#' @param get_forecasts Logical, whether to include forecast data.
#' @param region The geographical region of interest.
#' @param LE_output_folder Path to LOTOS-EUROS output directory.
#' @param brick_outputs Boolean, whether to save outputs as bricks.
#' @param variables Vector of variable names (e.g., c("no2", "o3")).
#' @param output Type of output required ("hourly", "daily").
#' @param output_dir Directory where output files will be saved.
#' @return None (writes raster files to disk).
#' @export
LE_nolbl_concs2brick <- function(date_folder,
                                 get_forecasts,
                                 region,
                                 LE_output_folder,
                                 brick_outputs,
                                 variables,
                                 output,
                                 output_dir){

  # date_folder = Sys.Date()-3
  # get_forecasts = TRUE
  # region = c("NL")
  # LE_output_folder = LEOPARD
  # variables = no_lab_varz2get
  # output = c("hourly")
  # output_dir = output_folder


  tryCatch({

    if(get_forecasts == TRUE){
      ## pick days to include - yesterday, today and tomorrow
      d8s_in <- as.character(c(date_folder-1, date_folder, date_folder+1, date_folder+2))
    } else{
      d8s_in <- as.character(date_folder)
    }

    ## loop through each day and output bricks
    for (d8 in d8s_in){
      tryCatch({
        #day_number <- paste0( "D", which(d8 == d8s_in))

        d8 <- ymd(d8)


        d8 <- ymd(d8)
        if(get_forecasts == TRUE){
          date2get <- gsub("-", "", as.character(d8))
        } else {
          date2get <- gsub("-", "", as.character(d8-1))
        }

        # how many hours are in the nc file, only 6 for the last day
        if(which(d8 == d8s_in) < 4){
          n_hrs = 24
        } else {
          n_hrs = 6
        }

        # want to get the future days model as it uses historic meteo rather than the forecast
        y <- year(d8)
        m <- sprintf("%02d", month(d8))
        d_1 <- sprintf("%02d", day(d8))


        ## READ IN METEO
        le_files <- paste0(LEOPARD, y, "/", m, "/", d_1, "/lotos-euros/C71-oper-Sectors-",region,"/output/LE_C71-oper-Sectors_conc-sfc_", date2get,".nc")

        lenolab <- nc_open(le_files)

        lon_ext <- ncvar_get(lenolab, "longitude_bnds")
        lat_ext <- ncvar_get(lenolab, "latitude_bnds")

        ## find the extent to generate a domain
        x_min <- min(lon_ext)
        x_max <- max(lon_ext)
        y_min <- min(lat_ext)
        y_max <- max(lat_ext)

        ##create a data frame to setup polygon generation
        df <- data.frame(X = c(x_min, x_max, x_max, x_min),
                         Y = c(y_max, y_max, y_min, y_min))

        ##generate a polygon of the area
        domain_ll <- df %>%
          st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
          dplyr::summarise(data = st_combine(geometry)) %>%
          st_cast("POLYGON")


        for (v in variables){

          tryCatch({

            var_info <- filter(varient_df, le_vars == v)

            lons <- ncvar_get(lenolab, "longitude")
            lats <- ncvar_get(lenolab, "latitude")

            nlon <- NROW(lons)
            nlat <- NROW(lats)

            var_in <- ncvar_get(lenolab, v, start = c(1,1,1,1), count = c(nlon, nlat, 1, n_hrs))

            TIME <- ncvar_get(lenolab, "time", start = c(1), count = c(n_hrs))

            time_since <- str_sub(lenolab$dim$time$units, 15, -14)

            d8_time <- lubridate::ymd(time_since) + lubridate::seconds(TIME)
            d8_char <- format(d8_time, "%Y-%m-%d %H:%M:%S")

            r_rd <- t(brick(var_in[1:nlon,1:nlat, 1:n_hrs]))

            ##define the extent
            crs(r_rd) <- 4326
            bb <- extent(domain_ll)
            extent(r_rd) <- bb
            r1 <- flip(r_rd, direction = 'y')
            r1_ug <- r1*var_info$omrekenfactor


            nam_r <- paste0(v, "_", date2get, "_ALL")
            # write out bricks
            if(any(output == "daily")){

              dir.create(paste0(output_dir, "concs/", region, "/daily/"))

              r_mean <- calc(r1_ug,mean)

              terra::writeRaster(r_mean, filename=paste0(output_dir, "concs/", region, "/daily/",nam_r,".TIF"), overwrite = TRUE)

            }

            # write out bricks
            if(any(output == "hourly")){

              dir.create(paste0(output_dir, "concs/", region, "/hourly/"))

              names(r1_ug) <- d8_char

              terra::writeRaster(r1_ug, filename=paste0(output_dir,"concs/", region, "/hourly/", nam_r, "_", date2get, ".TIF"), overwrite = TRUE)

            }


          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

        }



        print(paste0(nam_r, " ", date2get))

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})



    }

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#' Convert CH4 LOTOS EUROS Concentration NetCDF Data to Raster Bricks nad write to disk
#'
#' @param date_folder Date folder for processing.
#' @param get_forecasts Logical, whether to include forecast data.
#' @param LE_output_folder Path to LOTOS-EUROS output directory.
#' @param brick_outputs Boolean, whether to save outputs as bricks.
#' @param output Type of output required ("hourly", "daily").
#' @param output_dir Directory where output files will be saved.
#' @param country Optional country filter.
#' @param sector Optional sector filter.
#' @return None (writes raster files to disk).
#' @export
LE_lbl_CH4_concs2brick <- function(date_folder,
                               get_forecasts,
                               LE_output_folder,
                               brick_outputs,
                               output,
                               output_dir,
                               country = FALSE,
                               sector = FALSE){


  ## look up table of species names for LE, the factor convert from mol mol to ug/m3 and the title for tmap
  varient_df <- data.frame(le_vars = c("o3", "no2", "no", "nh3", "so2", 'hno3', "co", "ch4", "tpm25", "tpm10"),
                     omrekenfactor = c(1e12*0.048/24.4,  1e12*0.046005/24.4, 1e12*0.03001/24.4, 1e12*0.017031/24.4, 1e12*0.0640628/24.4, 1e12*0.06301/24.4,1e12*0.02801/24.4,1e12*0.01604/24.4, 1e9, 1e9),
                     title = c("'O'[3]*' ('* mu*'g/m'^3*')'", "'NO'[2]*' ('* mu*'g/m'^3*')'","'NO ('* mu*'g/m'^3*')'",
                               "'NH'[3]*' ('* mu*'g/m'^3*')'","'SO'[2]*' ('* mu*'g/m'^3*')'","'HNO'[3]*' ('* mu*'g/m'^3*')'",
                               "'CO'*' ('* mu*'g/m'^3*')'", "'CH'[4]*' ('* mu*'g/m'^3*')'","'PM'[2.5]*' ('* mu*'g/m'^3*')'",
                               "'PM'[10]*' ('* mu*'g/m'^3*')'"))

  v <- "ch4"
  region <- "EU"

  var_info <- filter(varient_df, le_vars == v)


    if(get_forecasts == TRUE){
      ## pick days to include - yesterday, today and tomorrow
      d8s_in <- as.character(c(date_folder-1, date_folder, date_folder+1))
    } else{
      d8s_in <- as.character(date_folder)
    }

    # want to get the future days model as it uses historic meteo rather than the forecast
    y <- year(date_folder)
    m <- sprintf("%02d", month(date_folder))
    d_1 <- sprintf("%02d", day(date_folder))

    ## loop through each day and output bricks
    for (d8 in d8s_in){
      tryCatch({
      #day_number <- paste0( "D", which(d8 == d8s_in))

      d8 <- ymd(d8)


      d8 <- ymd(d8)
      if(get_forecasts == TRUE){
        date2get <- gsub("-", "", as.character(d8))
      } else {
        date2get <- gsub("-", "", as.character(d8-1))
      }

      # how many hours are in the nc file, only 6 for the last day
      if(which(d8 == d8s_in) < 3){
        n_hrs = 24
      } else {
        n_hrs = 6
      }


      ## READ IN NC
      le_files <- paste0(LE_output_folder, y, "/", m, "/", d_1, "/lotos-euros/C71-CH4-oper-Countries-Sectors-EU/output/LE_C71-CH4-oper-Countries-Sectors_labelled-conc-3d_", date2get,".nc")

      lelab <- nc_open(le_files)

      lons <- ncvar_get(lelab, "longitude")
      lats <- ncvar_get(lelab, "latitude")

      nlon <- NROW(lons)
      nlat <- NROW(lats)

      lon_ext <- ncvar_get(lelab, "longitude_bnds")
      lat_ext <- ncvar_get(lelab, "latitude_bnds")

      ## find the extent to generate a domain
      x_min <- min(lon_ext)
      x_max <- max(lon_ext)
      y_min <- min(lat_ext)
      y_max <- max(lat_ext)

      # create a data frame to setup polygon generation
      df <- data.frame(X = c(x_min, x_max, x_max, x_min),
                       Y = c(y_max, y_max, y_min, y_min))

      # generate a polygon of the area
      domain_ll <- df %>%
        st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
        dplyr::summarise(data = st_combine(geometry)) %>%
        st_cast("POLYGON")

      # split into country and sector
      sector_names <- colsplit(ncvar_get(lelab, "labelnames"), " ", c("country", "sector"))

      # count
      nlabs <- NROW(sector_names)

      countries <- unique(sector_names$country)
      sectors <- unique(sector_names$sector)


      TIME <- ncvar_get(lelab, "time", start = c(1), count = c(n_hrs))

      time_since <- str_sub(lelab$dim$time$units, 15, -14)

      d8_time <- lubridate::ymd(time_since) + lubridate::seconds(TIME)
      d8_char <- format(d8_time, "%Y-%m-%d %H:%M:%S")
      #d8_epoch <- as.integer(d8_time)

      var_in <- ncvar_get(lelab, v, start = c(1,1,1,1,1), count = c(nlon, nlat, 1, n_hrs, nlabs))

      if(any(!country == FALSE) & any(sector == FALSE)){

        for (c in country){

        ind <- which(c == sector_names$country)

      # pick out labels
      var_subbed <- var_in[ , , ,ind]

      # sum labels selected
      var_summed <- apply(var_subbed, c(1, 2, 3), sum)

      # convert to raster brick
      r_brick <- t(brick(var_summed))

      # assign crs
      crs(r_brick) <- 4326

      # define extent
      extent(r_brick) <- extent(domain_ll)

      # flip along y axis
      r1 <- flip(r_brick, direction = 'y')

      # convert to concentration
      r1_ug <- r1*var_info$omrekenfactor

      # name the layers
      names(r1_ug) <- d8_char

      # write out bricks
      if(any(output == "daily")){

      dir.create(paste0(output_dir, "concs_ch4/countries/daily/"))

      r_mean <- calc(r1_ug,mean)

      terra::writeRaster(r_mean, filename=paste0(output_dir, "concs_ch4/countries/daily/",c,"_", date2get, ".TIF"), overwrite = TRUE)

      }

      # write out bricks
      if(any(output == "hourly")){

        dir.create(paste0(output_dir, "concs_ch4/countries/hourly/"))

      names(r1_ug) <- d8_char

      terra::writeRaster(r1_ug, filename=paste0(output_dir, "concs_ch4/countries/hourly/",c,"_", date2get, ".TIF"), overwrite = TRUE)

      }
      print(c)
        }

      }

      if(any(country == FALSE) & any(!sector == FALSE)){

        for (s in sector){

        ind <- which(s == sector_names$sector)

        # pick out labels
        var_subbed <- var_in[ , , ,ind]

        # sum labels selected
        var_summed <- apply(var_subbed, c(1, 2, 3), sum)

        # convert to raster brick
        r_brick <- t(brick(var_summed))

        # assign crs
        crs(r_brick) <- 4326

        # define extent
        extent(r_brick) <- extent(domain_ll)

        # flip along y axis
        r1 <- flip(r_brick, direction = 'y')

        # convert to concentration
        r1_ug <- r1*var_info$omrekenfactor

        # name the layers
        names(r1_ug) <- d8_char

        # write out bricks
        if(any(output == "daily")){

          dir.create(paste0(output_dir, "concs_ch4/sectors/daily/"), recursive = TRUE)

          r_mean <- calc(r1_ug,mean)

          terra::writeRaster(r_mean, filename=paste0(output_dir, "concs_ch4/sectors/daily/", s,"_", date2get, ".TIF"), overwrite = TRUE)

        }

        # write out bricks
        if(any(output == "hourly")){

          dir.create(paste0(output_dir, "concs_ch4/sectors/hourly/"), recursive = TRUE)

          names(r1_ug) <- d8_char

          terra::writeRaster(r1_ug, filename=paste0(output_dir, "concs_ch4/sectors/hourly/", s,"_", date2get, ".TIF"), overwrite = TRUE)

        }
        print(s)
        }

      }

      # show where up to
      print(paste0("reading from folder ", date_folder, " and downloading ", date2get))


      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }

    }



## Calc Z0
calc_z0 <- function(model_domain,
                    land_ids_file,
                    depac_file,
                    KBN_file,
                    lu_dir){
  
  # read in top10 key
  t10 <- read.csv(land_ids_file)
  
  # read in depac key
  depac <- read.csv(depac_file) %>% 
    arrange(long_name)
  
  # join TOP10 land ids with depac
  t10_depac_z0 <- left_join(t10, depac, c("depac" = "short"))
  
  ## get KBN grid
  KBN_grid <- readRDS(KBN_file)
  
  # make sure model domain is in same format as KBN_grid (rdnew)
  model_domain <- st_transform(model_domain, 28992)
  
  ## find which grid cells are within modelled area
  KBN_in <- KBN_grid[model_domain,]
  
  ## define top10 files to be read in
  rds_files <- paste0(lu_dir, 'top10/outputs/RDS/', unique(str_sub(KBN_in$KBN, 1,-4)), '.RDS')
  
  ## read in KBN gridded Z0 data for modelled area
  top10_in <- map_dfr(rds_files, readRDS)
  
  ## intersect 
  z0_poly <- st_intersection(top10_in,model_domain) %>% 
    left_join(t10_depac_z0, by = c('top10_id')) %>% 
    select(Z0, long_name, colour, geometry)
  
  ## define kappa and zr
  kappa = 0.41
  zr = 60
  
  ## find the unique land use parameters in the modelled area (domain)
  all_lu <- unique(z0_poly$long_name)
  ## define coefficient of drag to start with
  
  cd_ave = 0.0
  for (l in all_lu){
    
    lu <- filter(z0_poly, long_name == l)
    if(NROW(lu)>0){
      lu_frac <- as.numeric(sum(st_area(lu))/st_area(model_domain))
      cd = (kappa/log(zr/lu$Z0[1]))^2
      cd_ave = cd_ave+lu_frac*cd
    }
    
  }
  
  ## average roughness length over domain
  z0 = zr*exp(-kappa/sqrt(cd_ave))
  
  return(z0)
  
}

# library(tmap)
# plot_z0 <- function(model_domain,
#                     land_ids_file,
#                     depac_file,
#                     KBN_file,
#                     lu_dir,
#                     output_folder){
# 
#   model_domain = domain_polygon
#   land_ids_file = "DAT/defs/land_ids.csv"
#   depac_file = "DAT/defs/depac_z0.csv"
#   KBN_file = "DAT/defs/KBN_grid.RDS"
#   lu_dir = paste0(tsn, "RA-Data/Express/ra_express_modasuscratch_unix/projects/KIP/2023/land-use/")
#   output_folder <- "OUT/land_use/"
# 
#   t10 <- read.csv(land_ids_file)
# 
#   depac <- read.csv(depac_file) %>%
#     arrange(long_name)
# 
#   t10_depac_z0 <- left_join(t10, depac, c("depac" = "short"))
# 
#   ## get KBN grid
#   KBN_grid <- readRDS(KBN_file)
# 
#   model_domain <- st_transform(model_domain, 28992)
# 
#   ## find which grid cells are within modelled area
#   KBN_in <- KBN_grid[model_domain,]
# 
#   ## define files to be read in
#   rds_files <- paste0(lu_dir, 'top10/outputs/RDS/', unique(str_sub(KBN_in$KBN, 1,-4)), '.RDS')
# 
#   ## read in KBN gridded Z0 data for modelled area
#   top10_in <- map_dfr(rds_files, readRDS)
# 
#   ## intersect
#   z0_poly <- st_intersection(top10_in,model_domain) %>%
#     left_join(t10_depac_z0, by = c('top10_id')) %>%
#     select(Z0, long_name, colour, geometry)
# 
#   dir.create(output_folder, recursive = TRUE)
# 
#   ## simple plots of Z0 and LU classes associated
#   png(paste0(output_folder, 'depac_z0.png'))
#   plot(z0_poly['Z0'], main = paste0('average Z0 = ', mean(z0_poly$Z0)))
#   dev.off()
# 
#   bks <- seq(0,round(max(z0_poly$Z0),1),0.2)
# 
#   # create plot using tmap
#   tm_mod <- tm_shape(bg)+
#     tm_rgb()+
#     tm_shape(z0_poly)+
#     tm_polygons('Z0', alpha = 0.5, breaks = bks, lwd = 0.1, title = "Z0")+
#     tm_layout(legend.outside = FALSE, title.size = 1, legend.title.size = 0.001,
#               panel.label.color = 'white', panel.label.bg.color = 'black',bg.color = "white")+
#     tm_legend(position = c("left", "bottom"))
# 
#   tmap_save(tm_mod, paste0(output_folder, "01_Z0.png"), width = 2000, height = 2000, dpi = 300)
# 
#   ## load background map for plot for domain of data
#   bg <- basemaps::basemap_raster(ext=z0_poly, map_service = "carto", map_type = "light")
# 
#   all_lu <- unique(z0_poly$long_name)
# 
#   ## define kappa and zr
#   kappa = 0.41
#   zr = 60
# 
#   cd_ave = 0.0
#   for (l in all_lu){
# 
#     lu <- filter(z0_poly, long_name == l)
# 
#     if(NROW(lu)>0){
#       lu_frac <- as.numeric(sum(st_area(lu))/st_area(model_domain))
#       cd = (kappa/log(zr/lu$Z0[1]))^2
#       cd_ave = cd_ave+lu_frac*cd
#     }
# 
#     # create plot using tmap
#     tm_mod <- tm_shape(bg)+
#       tm_rgb()+
#       tm_shape(lu)+
#       tm_polygons('long_name', palette = lu$colour, lwd = 0.1, alpha = 0.5)+
#       tm_layout(legend.outside = FALSE, title.size = 1,
#                 panel.label.color = 'white', panel.label.bg.color = 'black',bg.color = "white", title = paste0('mean cd = ', round(cd_ave,5),
#                                                                                                                ' land use fraction = ', round(lu_frac,5)))+
#       tm_legend(position = c("left", "bottom"))
# 
#     assign(gsub(' ', '_',l), tm_mod)
# 
#     tmap_save(tm_mod, paste0(output_folder, l,".png"), width = 2000, height = 2000, dpi = 300)
# 
#   }
# 
#   tm_all <- tmap_arrange(Aquatic, Arable_Land, Broadleaf_deciduous_forest, Coniferous_evergreen_forest, Desert, Grassland, Permanent_crops, `Semi-natural_vegetation`, Urban, nrow = 4)
# 
#   tmap_save(tm_all, paste0(output_folder, 'all_lu.png'), width = 2000, height = 3200, dpi = 300)
# 
# 
# }


library(dplyr)
library(ncdf4)
library(stringr)
library(birk)
library(openair)
library(lubridate)
library(purrr)
library(sf)
library(raster)

rmw_plot_normalised <- function(df, colour = "#6B186EFF") {
  
  # Create plot
  plot <- ggplot2::ggplot(data = df)
  
  # Add se if present in data frame too
  if ("se" %in% names(df)) {
    
    # Plot the ribbon geom first
    plot <- plot +
      ggplot2::geom_ribbon(
        ggplot2::aes(x = date, ymin = value_predict - se, ymax = value_predict + se),
        alpha = 0.3,
        fill = "grey"
      )
    
  }
  
  # Overlay line
  plot <- plot +
    ggplot2::geom_line(
      ggplot2::aes(date, value_predict), colour = colour
    ) +
    ggplot2::geom_smooth(ggplot2::aes(date, value_predict), method = "loess", se = FALSE, colour = "black") +  # Smooth line
    ggplot2::theme_minimal() +
    ggplot2::ylab("Meteorologically normalised value") +
    ggplot2::xlab("Date")
  
  return(plot)
  
}

# this function writes openair plots to png files
openair2png <- function(plot, path, plot_name, width = 10500, height = 10000, units="px", res = 800){
  
  filename <- paste0(path, plot_name,".png")
  png(filename, width=width, height=height, units=units, res=res)
  print(plot)
  dev.off()
  
}


# this function imports the historic nc files from ICOS and makes a single dataframe for all the species and heights which
# matches the original format that most of the scripts work off
combine_historic_nc <- function(dir_location,
                                species,
                                heights){
  
  mol_conv <- data.frame(species = c("co2", "ch4", "n2o"),
                         mol_pp = c(10^6,10^9,10^9),
                         final_unit = c("ppm", "ppb", "ppb"))
  
  all_dat <- list()
  for (s in species){
    
    # find nc files with species that match
    filez <- list.files(scratch, s, full.names = TRUE)
    
    conv <- filter(mol_conv, species == s)
    
    # loop through heights
    for (h in heights){
      
      # pick out the species file for the specific height
      filez_h <- filez[grepl(paste0("-", h), filez)]
      
      # open the nc file
      icos <- nc_open(paste0(filez_h))
      
      # import time dimension
      TIME <- ncvar_get(icos, "time")
      
      # for some reason the ICOS data is at half hourly time steps, with the first 01:30. Assume it should be 01:00
      TIME_corrected <- TIME-3600/2
      
      # find the base date
      time_since <- gsub("T", " ", str_sub(icos$dim$time$units, 15, -2))
      
      # create formatted times
      d8_time <- lubridate::ymd_hms(time_since) + lubridate::seconds(TIME_corrected)
      
      # create a data frame of the dates. For some reason the dates are on the half past, so subtract 1800 seconds
      df <- data.frame(date = d8_time)
      
      ##extract variables and convert to ppm (divide by 10^6)
      df$value <- ncvar_get(icos, "value")*conv$mol_pp
      
      # assign the species name and height
      names(df) <- c("date", paste0(toupper(s), "_", h))
      
      # write to list
      all_dat[[paste0(s,h)]] <- df
      
      # let the user know where up to
      print(paste("imported",s,h))
    }
    
  }
  
  # left join them all by date to create one dataframe
  big_df <- purrr::reduce(all_dat, dplyr::left_join, by = c('date'))
  
  return(big_df)
  
}
dat_dir ="DAT/meteo/UTS_CT6/"
# this function extracts a variable from a ecmwf meteo file for specific coordinates
import_ecmwf <- function(dat_dir,
                         lat,
                         lon){
  
  dat_dirz <- list.files(dat_dir,
                         pattern = ".nc", full.names = FALSE)
  
  species = unique(str_sub(dat_dirz,1,-12))
  
  
  all_species = list()
  for (s in species){
    
    s_files = dat_dirz[grep(s, dat_dirz)]
    all_d8s <- list()
    for (d in s_files){
    
    ecmwf <- nc_open(paste0(dat_dir,d))
    
    species = names(ecmwf$var)[3]
    
    TIME <- ncvar_get(ecmwf, "valid_time")
    
    time_since <- str_sub(ecmwf$dim$valid_time$units, 15, -1)
    
    d8_time <- lubridate::ymd(time_since) + lubridate::seconds(TIME)
    
    met_lon <- ncvar_get(ecmwf, 'longitude')
    met_lat <- ncvar_get(ecmwf, 'latitude')
    
    lon_ind <- which.closest(met_lon, lon)
    lat_ind <- which.closest(met_lat, lat)
    
    ##extract no2 variables for entire domain for 1 time step and altitude
    b_met <- data.frame(date = d8_time,
                        nam = ncvar_get(ecmwf, species, start = c(lon_ind,lat_ind,1), count = c(1,1,NROW(d8_time))))
    
    names(b_met) = c("date", species)
    
    all_d8s[[d]] <- b_met
    #print(species)
    }
    
    all_species[[species]] = do.call(rbind,all_d8s)
    
  }
  
  param_df <- reduce(all_species, left_join)|>
    mutate(
      sun =  suncalc::getSunlightPosition(date = date, lat = lat, lon = lon),
      daylight = ifelse(sun$altitude > 0, "daylight", "nighttime")
    ) |> 
    select(-sun)
  
  row.names(param_df) <- NULL
  
  return(param_df)
  
  
}

get_country_domain <- function(countries = c('Slovenia', 'Greece'), return_format = 'sf'){
  
  dir.create("temp/countries", recursive = TRUE)
  download.file("https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/download/?format=geojson&timezone=Europe/Berlin&lang=en",
                destfile = "temp/countries/world_map.geojson")
  
  countries <- tolower(countries)
  
  countries <- (st_read("temp/countries/world_map.geojson")) %>% 
    select(name, geometry) %>% 
    mutate(name = tolower(name)) %>% 
    filter(name %in% c(countries))
  #st_geometry(countries) <- countries$geometry
  bb <- st_bbox(st_union(rbind(countries)))
  sf_dom <- bb %>% st_as_sfc() %>% st_as_sf() %>% mutate(lon_min = bb$xmin, lon_max = bb$xmax, lat_min = bb$ymin, lat_max = bb$ymax)
  
  if(return_format == 'bb'){
    return(bb)
  }
  if(return_format == 'sf'){
    return(sf_dom)
  }
  
}


find_noaa_domain <- function(domain, html_out = FALSE){
  #domain <- sub_domain
  ## import all meteo sites
  met_info <- getMeta(plot = FALSE)
  ## geo reference them
  met_sf <- st_as_sf(met_info, coords = c("longitude", "latitude"), crs = 4326)
  
  met_in <- met_sf[domain,]
  
  met_in$start_year <- year(met_in$begin)
  met_in$end_year <- year(met_in$end)
  
  if(html_out == TRUE){
    
    dir.create('plots/noaa', recursive = TRUE)
    ## set palette for site types
    pal_g <- colorFactor("Set3", reverse = FALSE, domain = met_in$end_year)
    
    lon <- mean(st_coordinates(sub_domain)[,1])
    lat <- mean(st_coordinates(sub_domain)[,2])
    
    ## create base base
    m <- leaflet() %>% 
      setView(lon, lat, zoom = 6) %>% 
      addProviderTiles('CartoDB.Positron')
    
    
    m <- m %>% addCircleMarkers(data = met_in, fillColor = ~pal_g(end_year), color = 'black', weight = 1,
                                opacity = 1.0, fillOpacity = 1.0,
                                popup = paste("station code:", met_in$code, "<br>",
                                              "station name:", met_in$station, "<br>",
                                              "country:", met_in$ctry, "<br>",
                                              "date start:", met_in$begin, "<br>",
                                              "date end:", met_in$end, "<br>",
                                              "elevation (m):", met_in$`elev(m)`, "<br>",
                                              "usaf:", met_in$usaf, "<br>"))
    
    
    m <- m %>% addLegend("bottomleft", pal=pal_g, values=met_in$end_year, opacity=1, title = "Latest year")
    
    print(paste0("saving html to: plots/noaa/",round(lon,1), "_", round(lat,1) ,".html"))
    ## use html widget saveWidget function to make a standalone html
    withr::with_dir('./', saveWidget(m, file = paste0("plots/noaa/obs_domain_",round(lon,1), "_", round(lat,1) ,".html")))
    
  }
  
  
  return(met_in)
  
}


find_noaa_sites <- function(sites, start_date, end_date){
  
  ## import all meteo sites
  met_info <- getMeta()
  ## geo reference them
  met_info <- st_as_sf(met_info, coords = c("longitude", "latitude"), crs = 4326)
  
  ## find nearest meteo station to site
  sites$nearest_NOAA <- met_info$code[st_nearest_feature(sites, met_info)]
  
  all_met_sites <- unique(sites$nearest_NOAA)
  
  ## find data range
  d8s <- seq(start_date, end_date)
  all_met <- list()
  
  for (m in all_met_sites){
    tryCatch({
      data_met <- importNOAA(m, d8s) %>% 
        select(code, station, date, latitude, longitude, elev, ws, wd, air_temp, atmos_pres,RH, ceil_hgt)
      
      all_met[[m]] <- data_met
      print(m)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  obs_met <- do.call(rbind, all_met)
  
  
}



ws_wd_plot <- function(ws_rast = NA, wd_rast = NA, start_hr = 1, end_hr = 4){
  
  
  # ws_rast <- crop(subset(ws_rast,start_hr:end_hr),domain)
  #wd_rast <- crop(subset(wd_rast,start_hr:end_hr), domain)
  
  ws_dat <- data.frame(rasterToPoints(ws_rast)) %>% 
    melt(c('x', 'y')) %>% 
    mutate(date = ymd_hm(gsub("X", "", variable)),
           ws = value) %>% 
    select(-variable, -value)
  
  wd_dat <- data.frame(rasterToPoints(wd_rast)) %>% 
    melt(c('x', 'y')) %>% 
    mutate(date = ymd_hm(gsub("X", "", variable)),
           wd = value) %>% 
    select(-variable, -value)
  
  
  ws_wd_dis <- ws_dat %>% 
    left_join(wd_dat, by = c('date', 'x', 'y')) %>% 
    select(lon = x, lat = y, ws, wd, date) %>% 
    mutate(wd = (wd*pi)/180) ## convert to radians
  
  
  framez <- NROW(unique(ws_wd_dis$date))
  
  g1 <- ggplot(ws_wd_dis, 
               aes(x = lon, 
                   y = lat, 
                   fill = ws, 
                   angle = wd, 
                   radius = scales::rescale(ws, c(.2, .8)/2))) +
    geom_raster() +
    geom_spoke(arrow = arrow(angle = 20, length = unit(1.0, 'mm'))) + 
    scale_fill_distiller(palette = "RdYlGn") +  
    coord_equal(expand = 0) + 
    theme(legend.position = 'bottom', 
          legend.direction = 'horizontal',
          panel.background = element_rect(fill='white', colour='black'))+
    transition_manual(date) +
    labs(title = "Period: {current_frame}")
  
  # animate in a two step process:
  animate(g1, height = 1600, width =1200, fps = 2, nframes = framez)
  
}


u_v_plot <- function(u_rast = NA, v_rast = NA, domain, start_hr = 1, end_hr = 4){
  
  
  if(is.na(ws_rast)){
    
    windDir <-function(u,v){
      (270-atan2(u,v)*180/pi)%%360 
    }
    
    ra_ws = sqrt(u_rast^2 + v_rast^2)
    ra_wd = windDir(u = u_rast, v = v_rast)
    
    names(ra_ws) <- u_rast@data@names
    names(ra_wd) <- u_rast@data@names
    
    ws_rast <- crop(subset(ra_ws,start_hr:end_hr),domain)
    wd_rast <- crop(subset(ra_wd,start_hr:end_hr), domain)
    
    
    ws_dat <- data.frame(rasterToPoints(ws_rast)) %>% 
      melt(c('x', 'y')) %>% 
      mutate(date = ymd_hm(gsub("X", "", variable)),
             ws = value) %>% 
      select(-variable, -value)
    
    wd_dat <- data.frame(rasterToPoints(wd_rast)) %>% 
      melt(c('x', 'y')) %>% 
      mutate(date = ymd_hm(gsub("X", "", variable)),
             wd = value) %>% 
      select(-variable, -value)
    
    
    ws_wd_dis <- ws_dat %>% 
      left_join(wd_dat, by = c('date', 'x', 'y')) %>% 
      select(lon = x, lat = y, ws, wd, date) %>% 
      mutate(wd = (wd*pi)/180) ## convert to radians
    
  } else {
    
    ws_rast <- crop(subset(ws_rast,start_hr:end_hr),domain)
    wd_rast <- crop(subset(wd_rast,start_hr:end_hr), domain)
    
    ws_dat <- data.frame(rasterToPoints(ws_rast)) %>% 
      melt(c('x', 'y')) %>% 
      mutate(date = ymd_hm(gsub("X", "", variable)),
             ws = value) %>% 
      select(-variable, -value)
    
    wd_dat <- data.frame(rasterToPoints(wd_rast)) %>% 
      melt(c('x', 'y')) %>% 
      mutate(date = ymd_hm(gsub("X", "", variable)),
             wd = value) %>% 
      select(-variable, -value)
    
    ws_wd_dis <- ws_dat %>% 
      left_join(wd_dat, by = c('date', 'x', 'y')) %>% 
      select(lon = x, lat = y, ws, wd, date) %>% 
      mutate(wd = (wd*pi)/180) ## convert to radians
    
  }
  
  framez <- NROW(unique(ws_wd_dis$date))
  
  g1 <- ggplot(ws_wd_dis, 
               aes(x = lon, 
                   y = lat, 
                   fill = ws, 
                   angle = wd, 
                   radius = scales::rescale(ws, c(.2, .8)/2))) +
    geom_raster() +
    geom_spoke(arrow = arrow(angle = 20, length = unit(1.0, 'mm'))) + 
    scale_fill_distiller(palette = "RdYlGn") +  
    coord_equal(expand = 0) + 
    theme(legend.position = 'bottom', 
          legend.direction = 'horizontal',
          panel.background = element_rect(fill='white', colour='black'))+
    transition_manual(date) +
    labs(title = "Period: {current_frame}")
  
  # animate in a two step process:
  animate(g1, height = 1600, width =1200, fps = 2, nframes = framez)
  
}


regrid <- function(raster2change, desired_grid, statistic = 'mean'){
  # raster2change <- conc_dom
  # desired_grid <- new_grid
  # statistic <- 'mean'
  nlayerz <- raster2change@data@nlayers
  layer_stack <- list()
  for (n in 1:nlayerz){
    
    rl <- subset(raster2change, n)
    
    resamp <- exactextractr::exact_resample(x = rl, y = desired_grid, fun = statistic)
    nam <- as.character(n)
    layer_stack[[nam]] <- resamp
    print(paste0('processing layer ', n))
    
  }
  
  new_brick <- brick(layer_stack)
  names(new_brick) <- names(raster2change)
  
  return(new_brick)
  
}


meteo_plot <- function(raster_in, all_layers = TRUE, statistic = 'mean', variable = 'temper', average = TRUE, start_hr = 1, end_hr = 64){
  
  if(all_layers == TRUE){
    raster_in <- subset(raster_in, start_hr, end_hr)
    
    bg <- basemaps::basemap_raster(ext=subset(raster_in,1), map_service = "carto", map_type = "light")
    
    var_df <- data.frame(var = c('temper', 'wspd_surf', 'wdir_surf'),
                         units = c('K', 'm/s', 'degrees'))
    
    legend_title <- filter(var_df, var == variable)
    
    if(NROW(legend_title)<1){
      legend_title <- variable
    }
    
    d8s <- as.character(ymd_hm(gsub("\\.", "_", gsub("X", "", names(raster_in)))))
    
    bks <- seq(min(raster_in@data@min), max(raster_in@data@max), by = (max(raster_in@data@max)-min(raster_in@data@min))/20)
    
    pal_A <- pals::jet(NROW(bks))
    
    tm_1 <- tm_shape(bg)+
      tm_rgb()+
      tm_shape(raster_in) +
      tm_raster(palette = pal_A,alpha = 0.6, style = 'cont', breaks = bks, title = parse(text = legend_title))+
      tm_layout(legend.outside = FALSE, frame = FALSE, legend.outside.position = 'left', title.size = 1, panel.labels = d8s,
                panel.label.color = 'white', panel.label.bg.color = 'black',
                bg.color = "white")+
      tm_legend(position = c("left", "top"))+
      tm_facets(nrow = 1, ncol = 1)
    
    
    dir.create("plots/meteo", recursive = TRUE)
    tmap_animation(tm_1, filename = paste0("plots/meteo/", variable, "_", d8s[1], "_", d8s[NROW(d8s)], "_.gif"), delay = 60)
    #tmap_animation(tm_1, filename = paste0("C:/Users/kellyb/OneDrive - TNO/topas_plots/", start_date, "_", end_date, "_", lbl, ".gif"),width = 2100, height = 2200, dpi = 300, delay = 60)
    
  } else {
    
    d8s <- as.character(ymd_hm(gsub("\\.", "_", gsub("X", "", names(raster_in)))))
    
    raster_in <- calc(raster_in, statistic)
    
    bg <- basemaps::basemap_raster(ext=subset(raster_in,1), map_service = "carto", map_type = "light")
    
    var_df <- data.frame(var = c('temper', 'wspd_surf', 'wdir_surf'),
                         units = c('K', 'm/s', 'degrees'))
    
    legend_title <- filter(var_df, var == variable)
    
    d8s <- as.character(ymd_hm(gsub("\\.", "_", gsub("X", "", names(raster_in)))))
    
    bks <- seq(min(raster_in@data@min), max(raster_in@data@max), by = (max(raster_in@data@max)-min(raster_in@data@min))/20)
    
    pal_A <- pals::jet(NROW(bks))
    
    tm_1 <- tm_shape(bg)+
      tm_rgb()+
      tm_shape(raster_in) +
      tm_raster(palette = pal_A,alpha = 0.6, style = 'cont', breaks = bks, title = parse(text = legend_title))+
      tm_layout(legend.outside = FALSE, frame = FALSE, legend.outside.position = 'left', title.size = 1, panel.labels = paste(variable,d8s),
                panel.label.color = 'white', panel.label.bg.color = 'black',
                bg.color = "white")+
      tm_legend(position = c("left", "top"))
    
    tmap_save(tm_1, paste0("plots/meteo/", variable, "_", statistic, "_", d8s[1], "_", d8s[NROW(d8s)], "_.gif"))
    
  }
  
}


find_aq_sites <- function(raster_list, type_def = 'background', area_def = c('urban', 'rural', 'suburban'), domain, species){
  #raster_list <- zoom_dat
  #species = 'tpm25' 
  # data_brick <- conc_dom
  #type_def = c('background')
  #area_def = c('urban')
  # variable = 'tpm25'
  # d8_range <- ymd_hm(gsub("\\.", "_", gsub("X", "", names(data_brick))))
  # min_d8 <- min(d8_range)
  # max_d8 <- max(d8_range)
  ## get all sites from database
  all_sites_sf <- list()
  for (s in species){
    
    saq_s <- varient_df %>% filter(le_specs == s) %>% select(saq_specs)
    
    rstr_in <- brick(raster_list[grepl(s, names(raster_list))])
    
    sitez <- saqgetr::get_saq_sites()
    st <- unique(sitez$site_type)
    min_d8 <- ymd_hm(gsub("\\.", "_", gsub("X", "", names(rstr_in))))[1]
    max_d8 <- ymd_hm(gsub("\\.", "_", gsub("X", "", names(rstr_in))))[rstr_in@data@nlayers]
    
    
    ## get more detailed data
    processes <- get_saq_processes()
    
    ## combine and convert to sf
    site_sums_sf <- sitez %>% 
      left_join(processes, by = 'site') %>% 
      filter(!is.na(latitude)) %>%
      filter(date_start.y < min_d8 & date_end.y > max_d8) %>% 
      mutate(variable_long = tolower(variable_long)) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
      filter(site_type %in% type_def) %>% 
      filter(site_area %in% area_def) %>% 
      filter(variable == saq_s$saq_specs) %>% 
      filter(period == 'hour')
    
    
    #db <- crop(rstr_in, sub_domain)
    db <- rstr_in
    
    ## define extent of brick
    r1 <- st_as_sf(st_as_sfc(st_bbox(subset(db,1))))
    ## filter points in area
    site_sf <- site_sums_sf[r1,]
    
    all_sites_sf[[s]] <- site_sf
    
    print(paste0('all ', s, ' sites imported between ', min_d8, " & ", max_d8))
    
  }
  
  all_sf <- do.call(rbind, all_sites_sf)
  
  return(all_sf)
  
}


plot_obs <- function(sites_sf, variable, write_out = FALSE){
  
  sitez <- sites_sf
  ## set palette for site types
  
  pal_st <- colorFactor("Set1", reverse = FALSE, domain = sitez$site_type)
  pal_sa <- colorFactor("Set3", reverse = FALSE, domain = sitez$site_area)
  
  lat <- mean(st_coordinates(sitez)[,2])
  lon <- mean(st_coordinates(sitez)[,1])
  ##plot on a map
  
  
  ## create base base
  m <- leaflet() %>% 
    setView(lon, lat, zoom = 6) %>% 
    addProviderTiles('CartoDB.Positron')
  
  
  m <- m %>% addCircleMarkers(data = sitez, fillColor = ~pal_st(site_type), color = 'black', weight = 1,
                              opacity = 1.0, fillOpacity = 1.0,
                              popup = paste("site code:", sitez$site, "<br>",
                                            "site name:", sitez$site_name, "<br>",
                                            "site type:", sitez$site_type, "<br>",
                                            "date start:", sitez$date_start.y, "<br>",
                                            "date end:", sitez$date_end.y, "<br>",
                                            "data source:", sitez$data_source.y, "<br>",
                                            "sampling process:", sitez$sampling_process, "<br>",
                                            "sampling process:", sitez$sampling_point, "<br>",
                                            "observation count:", sitez$observation_count.y, "<br>",
                                            "sampling period:", sitez$period, "<br>"), group = 'site type')
  
  
  m <- m %>% addLegend("bottomleft", pal=pal_st, values=sitez$site_type, opacity=1, title = "site type",
                       group = "site type")
  
  m <- m %>% addCircleMarkers(data = sitez, fillColor = ~pal_sa(site_area), color = 'black', weight = 1,
                              opacity = 1.0, fillOpacity = 1.0,
                              popup = paste("site code:", sitez$site, "<br>",
                                            "site name:", sitez$site_name, "<br>",
                                            "site type:", sitez$site_type, "<br>",
                                            "date start:", sitez$date_start.y, "<br>",
                                            "date end:", sitez$date_end.y, "<br>",
                                            "data source:", sitez$data_source.y, "<br>",
                                            "sampling process:", sitez$sampling_process, "<br>",
                                            "sampling process:", sitez$sampling_point, "<br>",
                                            "observation count:", sitez$observation_count.y, "<br>",
                                            "sampling period:", sitez$period, "<br>"), group = 'site area')
  
  
  m <- m %>% addLegend("bottomleft", pal=pal_sa, values=sitez$site_area, opacity=1, title = "site area",
                       group = "site area")
  
  m <- m %>% addLayersControl(baseGroups = c('site type', 'site area'),
                              options = layersControlOptions(collapsed = FALSE))
  
  m
  
  if(write_out == TRUE){
    
    dir.create("plots/concs/", recursive = TRUE)
    ## use html widget saveWidget function to make a standalone html
    withr::with_dir('plots/concs/', saveWidget(m, file = paste0(variable, ".html")))  
  }
  
  return(m)
  
  ## error catch in case one group gets an error
  ##  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


mod_obs_combine <- function(raster_list = NULL, raster_brick = NULL, sites_sf, species, time_step){
  # raster_list = NULL
  # raster_brick <- all
  # sites_sf <- sites_in
  # species <- 'tpm25'
  # time_step <- 'day'
  
  mod_obs_spec <- list()
  for (s in species){
    
    saq_s <- varient_df %>% filter(le_specs == s) %>% select(saq_specs)
    
    if(!is.null(raster_list)){
      rstr_in <- brick(raster_list[grepl(s, names(raster_list))])
    }
    
    if(!is.null(raster_brick)){
      rstr_in <- raster_brick
    }
    if(time_step == 'hour'){
      min_d8 <- ymd_hm(gsub("\\.", "_", gsub("X", "", names(rstr_in))))[1]
      max_d8 <- ymd_hm(gsub("\\.", "_", gsub("X", "", names(rstr_in))))[rstr_in@data@nlayers]
    }
    
    if(time_step == 'day'){
      min_d8 <- ymd(gsub("\\.", "_", gsub("X", "", names(rstr_in))))[1]
      max_d8 <- ymd(gsub("\\.", "_", gsub("X", "", names(rstr_in))))[rstr_in@data@nlayers]
    }
    
    s_sf <- filter(sites_sf, variable == saq_s$saq_specs)
    
    site_dat <- get_saq_observations(site = unique(s_sf$site), start = year(min_d8), end = year(max_d8), variable = saq_s$saq_specs) %>% 
      filter(date >= min_d8 & date <= max_d8) %>% 
      select(date, site, obs = value) %>% 
      timeAverage(time_step, type = 'site')
    
    sites_in <- unique(s_sf$site)
    
    df <- data.frame(t(raster::extract(rstr_in, s_sf)))
    
    names(df) <- s_sf$site
    
    if(time_step == 'hour'){
      rast_d8s <- ymd_hm(gsub("X", "", row.names(df)))
    }
    
    if(time_step == 'day'){
      rast_d8s <- ymd(gsub("X", "", row.names(df)))
    }
    
    dat <- data.frame(date = rast_d8s, df)
    
    dat_species <- dat %>% 
      melt(c('date')) %>% 
      select(date, site = 'variable', mod = value) %>% 
      left_join(site_dat, by = c('date', 'site')) %>% 
      mutate(species = s)
    
    mod_obs_spec[[s]] <- dat_species
    
    print(s)
    
  }
  
  all_mod_obs <- do.call(rbind, mod_obs_spec)
  
  return(all_mod_obs)
  
}

#library(aermod)
library(tidyverse)
library(ecmwfr)
library(lubridate)
library(ncdf4)

#dataset_options <- c("reanalysis-era5-land", "reanalysis-era5-single-levels")

# get the variables from ecmwf for download
variablez <- c("2m_dewpoint_temperature",
               "2m_temperature",
               "soil_temperature_level_1",
               "snow_cover",
               "surface_solar_radiation_downwards",
               "10m_u_component_of_wind",
               "10m_v_component_of_wind",
               "surface_pressure")


# path <- "./"
# 
# yrz <- c(2022,2023)
# 
# vzs <- get_era_sl(variables = variablez, dataset = "reanalysis-era5-land", api_key = "c95fc5f7-4e83-4e08-88ff-f02d03baa66e", path_out = path, years = yrz)


get_era_sl <- function(lon,
                       lat,
                       dataset = "reanalysis-era5-single-levels",
                       variables = "2m_dewpoint_temperature",
                       years,
                       api_key,
                       path_out) {
  
  #' Download ERA5 reanalysis data for a single location
  #'
  #' @param lon Longitude of the location (default: -2.283846)
  #' @param lat Latitude of the location (default: 53.439953)
  #' @param dataset Character string specifying the dataset. Options:
  #'     - "reanalysis-era5-single-levels"
  #'     - "reanalysis-era5-land"
  #' @param variables A vector of variables to download. Options:
  #'      - "2m_dewpoint_temperature",
  #'      -"2m_temperature",
  #'      -"soil_temperature_level_1",
  #'      - "snow_cover",
  #'      -"surface_solar_radiation_downwards",
  #'      -"10m_u_component_of_wind",
  #'      -"10m_v_component_of_wind",
  #'      -"surface_pressure"
  #' @param years A vector of years for which data is requested
  #' @param api_key API key for authentication (must be provided)
  #' @param path_out Directory to save output files (default: "out/")
  #' @return A dataframe containing the combined data for all variables and years
  
  if (is.null(api_key)) {
    stop("API key must be provided.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(path_out)) {
    dir.create(path_out, recursive = TRUE)
  }
  
  # Set the API key
  wf_set_key(key = api_key)
  
  all_varz <- list()
  
  for (v in variables) {
    tryCatch({
      
      all_yrz <- list()
      
      for (y in years) {
        
        start_date <- as.Date(paste0(y, "-01-01"))
        end_date <- as.Date(paste0(y, "-12-31"))
        kr8_d8 <- seq(from = start_date, to = end_date, by = "day")
        
        request <- list(
          dataset_short_name = dataset,
          product_type = "reanalysis",
          variable = v,
          year = as.character(y),
          month = unique(month(kr8_d8)),
          day = sprintf("%02d", 1:31),
          time = sprintf("%02d:00", 0:23),
          area = paste0(lat, "/", lon, "/", lat, "/", lon),
          format = "netcdf",
          target = paste0(v, "_", y, ".nc")
        )
        
        print(paste("requesting", v, "for", y))
        
        tryCatch({
          file <- wf_request(request = request, transfer = TRUE, path = path_out)
        }, error = function(e) {
          message("Error requesting variable ", v, " for year ", y, ": ", conditionMessage(e))
        })
        
      }
      
    }, error = function(e) {
      message("Error processing variable ", v, ": ", conditionMessage(e))
    })
    
  }
  
}

create_df_from_nc = function(){
  
  nc_file <- paste0(path_out, v, "_", y, ".nc")
  
  nf <- ncdf4::nc_open(nc_file)
  
  thyme <- ncvar_get(nf, "valid_time")
  d8 <- lubridate::ymd(str_sub(nf$dim$valid_time$units, 15, -1)) + lubridate::seconds(thyme)
  v_nam <- names(nf$var)[3]
  
  df <- data.frame(date = d8, value = ncdf4::ncvar_get(nf, v_nam))
  names(df) <- c("date", v_nam)
  
  ncdf4::nc_close(nf)
  
  all_yrz[[as.character(y)]] <- df
  
  message("Processed ", v, " for year ", y)
  
  
  if (length(all_yrz) > 0) {
    all_yrz_df <- do.call(rbind, all_yrz)
    write.csv(all_yrz_df, file = paste0(path_out, v_nam, "_", y, ".csv"), row.names = FALSE)
    all_varz[[v]] <- all_yrz_df
  }
  
  
  all_variables <- Reduce(function(x, y) dplyr::left_join(x, y, by = "date"), all_varz)
  
  return(all_variables)
  
}


library(mapedit)
library(mapview)
library(sf)

generate_grid <- function(lon = -2.283846,
                          lat = 53.439953){
  
  pt = st_sfc(
    st_point(c(lon,lat)),
    crs = 4326
  )
  
  domain <- mapview(pt) %>%
    editMap(title = "use the rectangle tool to draw a domain for ECMWF to download and select - done")
  
  domain_bb <- st_bbox(domain$finished)
  
  min_lon <- domain_bb[1]
  max_lon <- domain_bb[3]
  min_lat <- domain_bb[2]
  max_lat <- domain_bb[4]
  
  ecmwf_format_bb <- paste0(min_lat, "/", min_lon, "/", max_lat, "/", max_lon)
  
  return(ecmwf_format_bb)
  
}

# domain is in format (min_lat/min_lon/max_lat/max_lon) use function generate grid to define

get_era_domain <- function(domain,
                           dataset = "reanalysis-era5-single-levels",
                           variables = "2m_temperature",
                           years = yrz,
                           api_key = NULL,
                           path_out = "out/") {
  
  #' Download ERA5 reanalysis data for a domain
  #'
  #' @param lon Longitude of the location (default: -2.283846)
  #' @param lat Latitude of the location (default: 53.439953)
  #' @param dataset Character string specifying the dataset. Options:
  #'     - "reanalysis-era5-single-levels"
  #'     - "reanalysis-era5-land"
  #' @param variables A vector of variables to download. Options:
  #'      - "2m_dewpoint_temperature",
  #'      -"2m_temperature",
  #'      -"soil_temperature_level_1",
  #'      - "snow_cover",
  #'      -"surface_solar_radiation_downwards",
  #'      -"10m_u_component_of_wind",
  #'      -"10m_v_component_of_wind",
  #'      -"surface_pressure"
  #' @param years A vector of years for which data is requested
  #' @param api_key API key for authentication (must be provided)
  #' @param path_out Directory to save output files (default: "out/")
  #' @return A dataframe containing the combined data for all variables and years
  
  if (is.null(api_key)) {
    stop("API key must be provided.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(path_out)) {
    dir.create(path_out, recursive = TRUE)
  }
  
  # Set the API key
  wf_set_key(key = api_key)
  
  for (v in variables) {
    tryCatch({
      
      all_yrz <- list()
      
      for (y in years) {
        
        start_date <- as.Date(paste0(y, "-01-01"))
        end_date <- as.Date(paste0(y, "-12-31"))
        kr8_d8 <- seq(from = start_date, to = end_date, by = "day")
        
        months <- sprintf("%02d",unique(month(kr8_d8)))
        
        for (m in months){
          
          request <- list(
            dataset_short_name = dataset,
            product_type = "reanalysis",
            variable = v,
            year = as.character(y),
            month = m,
            day = sprintf("%02d", 1:31),
            time = sprintf("%02d:00", 0:23),
            area = domain,
            format = "netcdf",
            target = paste0(path_out, v, "_",m, y, ".nc")
          )
          
          tryCatch({
            file <- wf_request(request = request, transfer = TRUE, path = path_out)
          }, error = function(e) {
            message("Error requesting variable ", v, " for month ",m, " of year ", y, ": ", conditionMessage(e))
          })
          
          nc_file <- paste0(path_out, v, "_",m, y, ".nc")
          if (!file.exists(nc_file)) {
            next
          }
          
        }
        
      }
      
      
    }, error = function(e) {
      message("Error processing variable ", v, ": ", conditionMessage(e))
    })
  }
  
  if (length(all_varz) == 0) {
    stop("No data was successfully downloaded.")
  }
  
  all_variables <- Reduce(function(x, y) dplyr::left_join(x, y, by = "date"), all_varz)
  
  return(all_variables)
}

#library(aermod)
library(tidyverse)
library(ecmwfr)
library(lubridate)
library(ncdf4)

#dataset_options <- c("reanalysis-era5-land", "reanalysis-era5-single-levels")

# single levels https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
# land https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview



get_era <- function(bb,
                    dataset = "land",
                    location,
                    variables,
                    years,
                    months = seq(1,12),
                    api_key,
                    path_out) {
  
  #' Download ERA5 reanalysis data for a single location
  #'
  #' @param lon Longitude of the location (default: -2.283846)
  #' @param lat Latitude of the location (default: 53.439953)
  #' @param dataset Character string specifying the dataset. Options:
  #'     - "reanalysis-era5-single-levels"
  #'     - "reanalysis-era5-land"
  #' @param variables A vector of variables to download. Options:
  #'      - "2m_dewpoint_temperature",
  #'      -"2m_temperature",
  #'      -"soil_temperature_level_1",
  #'      - "snow_cover",
  #'      -"surface_solar_radiation_downwards",
  #'      -"10m_u_component_of_wind",
  #'      -"10m_v_component_of_wind",
  #'      -"surface_pressure"
  #' @param years A vector of years for which data is requested
  #' @param api_key API key for authentication (must be provided)
  #' @param path_out Directory to save output files (default: "out/")
  #' @return A dataframe containing the combined data for all variables and years
  
  if (is.null(api_key)) {
    stop("API key must be provided.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(paste0(path_out,location))) {
    dir.create(paste0(path_out,location), recursive = TRUE)
  }
  
  if(dataset == "land"){
    ds = "reanalysis-era5-land"
  }
  if(dataset == "single-level"){
    ds = "reanalysis-era5-single-levels"
  }
  
  # Set the API key
  wf_set_key(key = api_key)
  
  for (y in years) {
    
    for (m in months){
      
      for (v in variables) {
        
        tryCatch({
          
          request <- list(
            dataset_short_name = ds,
            product_type = "reanalysis",
            variable = v,
            year = as.character(y),
            month = m,
            #month = unique(month(kr8_d8)),
            day = sprintf("%02d", 1:31),
            time = sprintf("%02d:00", 0:23),
            area = ecmwf_format_bb,
            format = "netcdf",
            download_format = "unarchived",
            target = paste0(v, "_",sprintf("%02d",m),"_", y, ".nc")
          )
          
          print(paste("requesting", v, "for ",y, " in ", sprintf("%02d",m)))
          
          tryCatch({
            file <- wf_request(request = request, transfer = TRUE, path = paste0(path_out,location, "/"))
          }, error = function(e) {
            message("Error requesting variable ", v, " for month ", m, " year ", y, ": ", conditionMessage(e))
          })
          
        }, error = function(e) {
          message("Error requesting variable ", v, " for month ", m, " year ", y, ": ", conditionMessage(e))
        })
        
      }
      
    }
    
    
    
  }
  
}

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
library(rAHNextract)

select <- dplyr::select

# allows user to select a custom area based on latitude and longitude coordinates
get_custom_domain = function(lat,lon){
  
  latitude <- lat
  longitude <- lon
  
  ##convert to sf point
  location <- st_point(c(longitude, latitude))
  location <- st_sfc(location, crs = 4326)
  
  
  # ##Plot the point sources around site to get information, once drawn click FINISH this creates a "simple features" frame (sf)
  modelled_area <- mapview(location, map.types = c("OpenStreetMap", "Esri.WorldTopoMap", "Esri.WorldImagery", "Esri.WorldShadedRelief")) %>%
    editMap(title = "Use the rectangle tool to draw the area to be modelled")
  
  ##Once finished create variable to be plotted
  modelled_area <- modelled_area$finished
  
  return(modelled_area)
  
}

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
  
  ahn <- ahn_area(
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


