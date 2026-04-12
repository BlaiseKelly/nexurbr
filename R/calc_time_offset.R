#' Calculate UTC time offset for a location and set of dates
#'
#' Looks up the time zone for a given coordinate and returns the UTC offset
#' in hours for each date-time, accounting for daylight saving transitions.
#'
#' @param dates POSIXct vector. Date-times to calculate offsets for.
#' @param lat Numeric. Latitude in decimal degrees.
#' @param lon Numeric. Longitude in decimal degrees.
#'
#' @return Numeric vector of UTC offsets in hours (e.g. 1 for CET, 2 for CEST).
#'
#' @examples
#' dates = as.POSIXct(c("2024-06-15 12:00", "2024-12-15 12:00"), tz = "UTC")
#' calc_time_offset(dates, lat = 52.11, lon = 5.03)
#'
#' @export
calc_time_offset = function(dates, lat, lon) {
  tz = lutz::tz_lookup_coords(lat = lat, lon = lon, method = "fast", warn = FALSE)
  offset = lutz::tz_offset(dates, tz = tz)$utc_offset_h
  return(offset)
}
