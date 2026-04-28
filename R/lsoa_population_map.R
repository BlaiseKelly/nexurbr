
# binning population from lsoa data

population_map = function(lsoa_poly,
                          population,
                             fill_palette = "tol.incandescent",
                             backgroup_map = "CartoDB.Positron"){
  
  
  lsoa_in = IMD_poly |> 
    left_join(pop_want, by = c("LSOA21CD" = "LSOA.2021.Code"))
  
  
  lsoa_pop = readRDS("OUT/DAT/pop_age_sex.RDS")
  
  lsoa_score_pop = left_join(lsoa_scores,lsoa_pop, by = c("LSOA21CD" = "LSOA.2021.Code"))
  
  backgroup_map = "CartoDB.Positron"
  
  tmap_mode("view")
  # initial geometry extent to start map off
  tm1 = tm_shape(geo_buf) +
    tm_polygons(fill_alpha = 0,
                col_alpha = 0,
                popup.vars = FALSE,
                group = "LAs", 
                group.control = "hide")
  
  geo_df = expand.grid(age = unique(lsoa_pop$age_band),sex = unique(lsoa_pop$sex)) |> 
    mutate(age_sex = paste(sex,age))
  geo_type = geo_df$age_sex
  # create function to loop through
  tm1 = reduce(geo_type, function(map_obj, ct) {
    
    demo_df = filter(geo_df,age_sex == ct)
    
    pop_df = filter(lsoa_score_pop,age_band == demo_df$age & sex == demo_df$sex)
    
    map_obj = map_obj +
      tm_shape(pop_df) +
      tm_polygons(fill = "pop",
                  fill.scale = tm_scale_intervals(n = 10),
                  #fill.legend = tm_legend(show = FALSE),
                  lwd = 0.5,
                  group = ct,
                  group.control = "radio")
    
    map_obj
    
    
  }, .init = tm1)
  
  map_title = "18 and under population of LSOA areas split by age and sex"
  
  tm1 <- tm1 + 
    tm_basemap(backgroup_map, group.control = "none")+ 
    tm_title(map_title)+
    # tm_add_legend(
    #   type = "polygons",
    #   labels = names(qpal),
    #   fill = unname(qpal)
    # )+
    tm_view(control.collapse = FALSE,overlay.groups = geo_type[1])
  
  tm1
  
  
}
