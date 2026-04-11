
# get nl lidar using the WCS service
get_nl_lidar = function(domain,product = c("DSM", "DTM"),res = c(0.5,5)){

  domain_rdnew = domain |>
    st_transform(28992) |>
    st_bbox()

  ahn = rAHNextract::ahn_area(
    bbox = domain_rdnew,
    resolution = res,
    dem = product
  )

  return(ahn)

}


