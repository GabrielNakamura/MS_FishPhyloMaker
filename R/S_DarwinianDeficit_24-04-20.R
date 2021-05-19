
# libraries and data ------------------------------------------------------
library(rgdal)     
library(ggplot2)  
library(ggmap)   

phylo_drainages <- readRDS(here::here("output", "phylo_drainages.rds"))
shapefile <- rgdal::readOGR(here::here("data", "Basin042017_3119.shp"))
PD_deficit_all <- readRDS(here::here("output", "PD_deficit_all.rds"))


# formating shapefile for plot --------------------------------------------

for(i in 1:length(shapefile@polygons)){
  shapefile@polygons[[i]]@ID <- shapefile@data$BasinName[i]
}
shapefile_df <- fortify(shapefile)

PD_dataframe <- data.frame(deficit = unlist(PD_deficit_all), id = names(PD_deficit_all))
deficit_value <- unlist(lapply(shapefile_df$id, function(x) 
  PD_dataframe[which(x == PD_dataframe$id), "deficit"]
)
)
data_map_PDeficit <- data.frame(shapefile_df, deficit = deficit_value)


# plotting map ------------------------------------------------------------
wrld <- map_data("world")
w <- ggplot(wrld) +
  geom_polygon(aes(x = long, y = lat, group = group), alpha = 0.3) +
  geom_polygon(data = data_map_PDeficit, aes(x = long, y = lat, group = id, 
                                             fill = deficit),
                            colour = "black", size = .15) + scale_fill_viridis_c() +
  labs(x = "Longitude", y = "Latitude")
w

cowplot::save_plot(filename = here::here("output", "images", "world_map.png"), 
                   w, base_width = 10, base_height = 5)
