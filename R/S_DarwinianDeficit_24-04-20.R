
# libraries and data ------------------------------------------------------
library(rgdal)     # R wrapper around GDAL/OGR
library(ggplot2)   # for general plotting
library(ggmap)    # for fortifying shapefiles
phylo_drainages <- readRDS(here::here("output", "phylo_drainages.rds"))
shapefile <- rgdal::readOGR(here::here("data", "Basin042017_3119.shp"))

# Next the shapefile has to be converted to a dataframe for use in ggplot2
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

# Now the shapefile can be plotted as either a geom_path or a geom_polygon.
# Paths handle clipping better. Polygons can be filled.
# You need the aesthetics long, lat, and group.
p <- ggplot(data = data_map_PDeficit) +
  geom_polygon(aes(x = long, y = lat, group = id, fill = deficit),
               colour = "black")
quartz()
Map_darwinian_shortfall <- p + scale_fill_viridis_c()
pdf("Map_darwinian.pdf", width = 15, height = 5)
print(Map_darwinian_shortfall)
dev.off()
