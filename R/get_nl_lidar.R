#' Download Dutch AHN LiDAR elevation data
#'
#' Retrieves DSM or DTM raster data from the Dutch AHN (Actueel Hoogtebestand
#' Nederland) Web Coverage Service for a given study area.
#'
#' @param domain An \code{sf} object defining the area of interest (any CRS;
#'   will be transformed to EPSG:28992 internally).
#' @param product Character. Either \code{"DSM"} or \code{"DTM"}.
#' @param res Numeric. Raster resolution in metres. AHN supports 0.5 m and 5 m.
#'
#' @return A \code{SpatRaster} (or raster object, depending on \code{rAHNextract}
#'   version) of the requested elevation product.
#'
#' @examples
#' \dontrun{
#' domain <- get_custom_domain()
#' dsm <- get_nl_lidar(domain, product = "DSM", res = 0.5)
#' dtm <- get_nl_lidar(domain, product = "DTM", res = 5)
#' }
#'
#' @export
get_nl_lidar <- function(domain, product = c("DSM", "DTM"), res = c(0.5, 5)) {

  product <- match.arg(product)
  res <- res[1]

  domain_rdnew <- domain |>
    sf::st_transform(28992) |>
    sf::st_bbox()

  ahn <- rAHNextract::ahn_area(
    bbox = domain_rdnew,
    resolution = res,
    dem = product
  )

  return(ahn)
}
