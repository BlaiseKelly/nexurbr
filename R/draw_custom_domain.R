#' Draw a custom study area polygon interactively
#'
#' Opens an interactive map centred on the given coordinates, allowing the user
#' to draw a polygon (typically a rectangle) defining a study domain. Requires
#' an interactive R session.
#'
#' @param latitude Numeric. Centre latitude in decimal degrees (default: 52.114936, Utrecht).
#' @param longitude Numeric. Centre longitude in decimal degrees (default: 5.026817, Utrecht).
#'
#' @return An \code{sf} object containing the drawn polygon(s) in EPSG:4326.
#'
#' @examples
#' \dontrun{
#' domain <- draw_custom_domain()
#' domain <- draw_custom_domain(latitude = 51.4545, longitude = -2.5879)
#' }
#'
#' @export
draw_custom_domain <- function(latitude = 52.114936, longitude = 5.026817) {
  
  
  location <- sf::st_point(c(longitude, latitude))
  location <- sf::st_sfc(location, crs = 4326)
  
  
  modelled_area <- mapview::mapview(
    location,
    map.types = c("OpenStreetMap", "Esri.WorldTopoMap",
                  "Esri.WorldImagery", "Esri.WorldShadedRelief")
  ) %>%
    mapedit::editMap(title = "Use the rectangle tool to draw the area to be modelled")
  
  modelled_area <- modelled_area$finished
  
  return(modelled_area)
}