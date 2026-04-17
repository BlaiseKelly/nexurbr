#' Download LSOA 2021 boundary geometries
#'
#' Retrieves Lower Layer Super Output Area (LSOA) 2021 boundary polygons from
#' either the \code{geographr} package or a GeoPackage hosted on GitHub
#' (originally sourced from the ONS Open Geography Portal).
#'
#' @param provider Character. One of \code{"geographr"} (default) or
#'   \code{"ons"}. Controls where boundaries are loaded from.
#' @param lsoa_code Unquoted column name for the LSOA code field when
#'   \code{provider = "ons"} (passed via tidy evaluation).
#' @param lsoa_name Unquoted column name for the LSOA name field when
#'   \code{provider = "ons"}.
#' @return An \code{sf} data frame with columns \code{lsoa21_code},
#'   \code{lsoa21_name}, and \code{geometry}.
#' @examples
#' \dontrun{
#' lsoa_sf <- get_lsoa21_boundaries(provider = "geographr")
#' lsoa_sf <- get_lsoa21_boundaries(provider = "ons",
#'                                   lsoa_code = LSOA21CD,
#'                                   lsoa_name = LSOA21NM)
#' }
#' @export
get_lsoa21_boundaries <- function(provider = c("geographr","ons"), lsoa_code, lsoa_name){

  if(provider == "geographr"){
    lsoa_geo = geographr::boundaries_lsoa21 |>
      dplyr::select(lsoa21_code,lsoa21_name,geometry)
  }

  if(provider == "ons"){

    # unreliable ONS server switched to github release
    lsoa_url = "https://github.com/BlaiseKelly/stats19_stats/releases/download/boundaries-v1.0/lsoa_boundaries.gpkg"

    # import and normalise the names
    lsoa_geo <- sf::st_read(lsoa_url) |>
      dplyr::select(lsoa21_code = {{lsoa_code}},lsoa21_name = {{lsoa_name}},geometry = geom)

    # make sure it knows the geometry column
    st_geometry(lsoa_geo) = lsoa_geo$geometry

  }

  return(lsoa_geo)

}
