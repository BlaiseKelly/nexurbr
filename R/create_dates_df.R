

#' Create a dataframe of date-time sequences with temporal attributes
#'
#' Generates a regular time series between two dates and adds columns for
#' hour, day of week, day of year, and month. Useful for joining to
#' temporally-varying emissions or meteorological profiles (e.g. AERMOD
#' variable emission factors).
#'
#' @param start_date Character. Start date-time string parseable by
#'   \code{\link[base]{as.POSIXct}} (e.g. \code{"2024-01-01"}).
#' @param end_date Character. End date-time string.
#' @param time_zone Character. Time zone string (e.g. \code{"UTC"},
#'   \code{"Europe/Amsterdam"}). Note: currently unused internally;
#'   all parsing is done in UTC.
#' @param timestep Character. Time step passed to \code{\link[base]{seq.POSIXt}}
#'   (e.g. \code{"hour"}, \code{"30 min"}).
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{date}{POSIXct timestamp.}
#'   \item{variable}{Hour label in format \code{X01}–\code{X24}.}
#'   \item{dow}{Lowercase abbreviated day of week (e.g. \code{"mon"}).}
#'   \item{hr}{Integer hour 1–24.}
#'   \item{doy}{Integer day of year.}
#'   \item{month}{Lowercase abbreviated month (e.g. \code{"jan"}).}
#' }
#'
#' @examples
#' dates = create_dates_df("2024-01-01", "2024-01-07", "UTC", "hour")
#' head(dates)
#'
#' @export
create_dates_df = function(start_date, end_date, time_zone, timestep) {

  start_d8 = as.POSIXct(start_date, tz = "UTC")
  end_d8 = as.POSIXct(end_date, tz = "UTC")

  d8_df = data.frame(
    date = seq(from = start_d8, to = end_d8, by = timestep)
  ) |>
    dplyr::mutate(variable = paste0("X", sprintf("%02d", lubridate::hour(date) + 1)))

  d8_df$dow = tolower(as.character(lubridate::wday(d8_df$date, label = TRUE)))
  d8_df$hr = lubridate::hour(d8_df$date) + 1
  d8_df$doy = lubridate::yday(d8_df$date)
  d8_df$month = tolower(as.character(lubridate::month(d8_df$date, label = TRUE)))

  return(d8_df)
}
