#' Download GB LiDAR elevation data from Environment Agency
#'
#' Retrieves DSM or DTM raster data from the Environment Agency's LiDAR
#' composite catalogue for England via the \code{gblidar} package.
#'
#' @param domain An \code{sf} object defining the area of interest (any CRS;
#'   will be transformed to EPSG:27700 internally).
#' @param product Character. Either \code{"DSM"} or \code{"DTM"}.
#' @param product_type Character. Either \code{"elevation"} or \code{"hillshade"}.
#'   Currently unused but reserved for future support.
#' @param res Numeric. Resolution in metres (1 or 2). Currently unused;
#'   resolution is determined by the EA composite dataset.
#'
#' @return A \code{SpatRaster} of the requested elevation product.
#'
#' @examples
#' \dontrun{
#' domain <- get_custom_domain(latitude = 51.45, longitude = -2.59)
#' dsm <- get_gb_lidar(domain, product = "DSM")
#' dtm <- get_gb_lidar(domain, product = "DTM")
#' }
#'
#' @export
get_gb_lidar <- function(domain,
                         product = c("DSM", "DTM"),
                         product_type = c("elevation", "hillshade"),
                         res = c(1, 2)) {

  product <- match.arg(product)
  product_type <- match.arg(product_type)

  if (rlang::is_installed("terra")) {
    options(gblidar.out_raster_type = "SpatRaster")
  }

  options(gblidar.progress = FALSE)

  domain_os <- domain |>
    sf::st_transform(27700)

  r <- gblidar::eng_composite(domain_os, product = product)

  return(r)
}
