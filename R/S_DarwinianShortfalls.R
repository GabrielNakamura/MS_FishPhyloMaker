
# libraries and data ------------------------------------------------------
library(ggplot2)  
library(tidyverse)
library(sf)
library(cartogram)
library(patchwork)

phylo_drainages <- readRDS(here::here("output", "phylo_drainages.rds"))
PD_deficit_all <- readRDS(here::here("output", "PD_deficit_all.rds"))
shapefile <- sf::read_sf(here::here("data"), as_tibble = T)

# formatting data ----------------------------------------------------------

PD_deficit_org <- data.frame(BasinName = shapefile$BasinName, 
                             PD_deficit = unlist(PD_deficit_all)[match(shapefile$BasinName, 
                                                                       names(unlist(PD_deficit_all)
                                                                             )
                                                                       )
                                                                 ]
                             )

# formating shapefiles for plot --------------------------------------------

sf_darwinian <- 
  shapefile %>%
  st_transform(crs = "+proj=robin") %>% 
  left_join(PD_deficit_org)


# plotting map ------------------------------------------------------------

map_PD_deficit_wrld <- ggplot() +
  geom_sf(data = sf_darwinian, aes(geometry = geometry, 
                                   fill = PD_deficit),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "BluYl", 
                                  direction = -1, 
                                  limits = c(0, 1),  ## max percent overall
                                  breaks = seq(0, 1, by = .1),
                                  labels = glue::glue("{seq(0, 1, by = 0.1)}")) +
  guides(fill = guide_colorbar(barheight = unit(5, units = "mm"),  
                               barwidth = unit(100, units = "mm"),
                               direction = "horizontal",
                               ticks.colour = "grey20",
                               title.position = "top",
                               label.position = "bottom",
                               title.hjust = 0.5)) +
  labs(fill = expression(PD[inserted]/PD[total])) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "grey70")) 
  


# saving plot -------------------------------------------------------------

ggsave(filename = here::here("output", "images", "world_map_DarwinianShort.png"), plot = map_PD_deficit_wrld,
       width = 10, height = 5, 
       dpi = 700)
ggsave(here::here("output", "images", "world_map_DarwinianShort.pdf"), 
       width = 10, height = 5, device = cairo_pdf)

