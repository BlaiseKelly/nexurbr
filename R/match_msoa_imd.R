#' Compute population-weighted IMD decile for MSOAs
#'
#' Aggregates LSOA-level Index of Multiple Deprivation (IMD 2025) data to MSOA
#' level using a population-weighted mean of IMD deciles. Downloads the IMD
#' GeoPackage from GitHub if not supplied.
#'
#' @param IMD_lsoa_data Optional \code{sf} data frame of LSOA-level IMD data.
#'   If \code{NULL} (the default), the function downloads it from GitHub.
#' @return A data frame with columns \code{MSOA21CD},
#'   \code{imd_weighted} (population-weighted mean IMD decile), and
#'   \code{n_lsoa} (number of LSOAs in each MSOA).
#' @examples
#' \dontrun{
#' msoa_imd <- match_msoa_imd()
#' msoa_imd <- match_msoa_imd(IMD_lsoa_data = my_imd_sf)
#' }
#' @export
match_msoa_imd = function(IMD_lsoa_data = NULL, keep_geometry = NULL){

  if(is.null(IMD_lsoa_data)){

    IMD_lsoa_data = st_read("https://github.com/BlaiseKelly/IMD/releases/download/LSOA_IMD2025/LSOA_IMD2025_WGS84_-4854136717238973930.gpkg")

  }

  lsoa_geo = get_lsoa21_boundaries(provider = "geographr") |>
    st_transform(27700)

  lsoa_cent = st_centroid(lsoa_geo)

  lsoa_pop = read.csv("https://github.com/BlaiseKelly/lsoa_ons_population/releases/download/v0.1.1/lsoa21_pop_tot_2011_2024.csv") |>
    select(lsoa21cd,pop = X2024)

  msoa_geo = st_read("https://github.com/BlaiseKelly/stats19_stats/releases/download/msoa_boundaries-v1.0/msoa.gpkg") |>
    st_transform(27700) |>
    select(MSOA21CD,geom)


  msoa_lsoa = st_join(msoa_geo,lsoa_cent) |>
    select(MSOA21CD,lsoa21_code)

  if(is.null(keep_geometry)){
    st_geometry(msoa_lsoa) = NULL
  }

  # read in hcl names
  msoa_nm = read.csv("https://houseofcommonslibrary.github.io/msoanames/MSOA-Names-2.2.csv") |>
    distinct(msoa21cd,msoa21hclnm, localauthorityname)

  msoa_imd <- IMD_2025 |>
    st_set_geometry(NULL) |>
    left_join(lsoa_pop, by = c("LSOA21CD" = "lsoa21cd")) |>
    left_join(msoa_lsoa, by = c("LSOA21CD" = "lsoa21_code")) |>
    group_by(MSOA21CD) %>%
    summarise(
      imd_weighted = weighted.mean(IMDDecil, w = pop, na.rm = TRUE),
      population_2024 = sum(pop),
      n_lsoa = n()
    ) |>
    left_join(msoa_nm, by = c("MSOA21CD" = "msoa21cd"))


  return(msoa_imd)

}
