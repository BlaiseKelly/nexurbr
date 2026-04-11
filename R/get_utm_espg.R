#' Determine the UTM EPSG code for a given location
#'
#' Calculates the appropriate UTM zone and returns the corresponding EPSG code
#' for any point on earth. Useful for dynamically selecting a local metric
#' projection for buffering, distance calculations, or area measurements.
#'
#' @param longitude Numeric. Longitude in decimal degrees.
#' @param latitude Numeric. Latitude in decimal degrees.
#'
#' @return Integer. The EPSG code for the UTM zone (326xx for northern
#'   hemisphere, 327xx for southern).
#'
#' @examples
#' get_utm_epsg(5.03, 52.11)   # Utrecht -> 32631
#' get_utm_epsg(-2.59, 51.45)  # Bristol -> 32630
#' get_utm_epsg(2.17, 41.39)   # Barcelona -> 32631
#'
#' @export
get_utm_epsg = function(longitude, latitude) {
  zone = floor((longitude + 180) / 6) + 1
  if (latitude >= 0) 32600 + zone else 32700 + zone
}
