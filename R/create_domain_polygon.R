#' Create a rectangular study domain from a point
#'
#' Buffers a point location to produce a polygon domain, useful as an input
#' to data retrieval functions that require an area rather than a single
#' coordinate (e.g. ECMWF ERA5 requests, LiDAR downloads).
#'
#' @param latitude Numeric. Latitude in decimal degrees (default: 52.114936, Utrecht).
#' @param longitude Numeric. Longitude in decimal degrees (default: 5.026817, Utrecht).
#' @param buffer_dist Numeric. Buffer distance in metres (applied in EPSG:27700).
#'
#' @return An \code{sfc} polygon in EPSG:4326.
#'
#' @examples
#' domain <- create_domain_polygon()
#' domain <- create_domain_polygon(latitude = 51.45, longitude = -2.59, buffer_dist = 500)
#'
#' @export
create_domain_polygon = function(latitude = 52.114936, longitude = 5.026817, buffer_dist = 40) {

  utm = get_utm_epsg(longitude,latitude)

  pt = sf::st_sfc(
    sf::st_point(c(longitude, latitude)),
    crs = 4326
  ) |>
    sf::st_transform(utm)

  domain = pt |>
    sf::st_buffer(buffer_dist) |>
    sf::st_transform(4326)

  return(domain)
}
