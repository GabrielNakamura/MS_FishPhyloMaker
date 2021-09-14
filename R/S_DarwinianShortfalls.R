
# libraries and data ------------------------------------------------------
library(ggplot2)  
library(tidyverse)
library(sf)
library(cartogram)
library(patchwork)

phylo_drainages <- readRDS(here::here("output", "phylo_drainages.rds"))
PD_deficit_all <- readRDS(here::here("output", "PD_deficit_all.rds"))
PD_deficit_order <- readRDS(here::here("output", "PD_deficit_perOrder.rds"))
shapefile <- sf::read_sf(here::here("data"), as_tibble = T)


# setting theme -----------------------------------------------------------

theme_update(
  plot.title = element_text(family = "Bitter", size = 32, hjust = 0.5),
  legend.position = "bottom",
  legend.title = element_text(family = "Bitter", color = "black", face = "bold", size = 12),
  legend.text = element_text(family = "Roboto Mono", color = "black", size = 10)
)

theme_set(theme_bw())

# formatting data deficit ----------------------------------------------------------

PD_deficit_all <- unlist(PD_deficit_all[match(rownames(PD_deficit_order), names(PD_deficit_all))])

df_deficit <- data.frame(BasinName = rownames(PD_deficit_order), PD_deficit_all, PD_deficit_order)



# formatting shapefiles for plot --------------------------------------------

sf_darwinian <- 
  shapefile %>%
  st_transform(crs = "+proj=robin") %>% 
  left_join(df_deficit)

sf_darwinian_rmzero <- sf_darwinian[- which(sf_darwinian$PD_deficit_all == 0.00000000), ]

# plotting map ------------------------------------------------------------

map_PD_deficit_wrld <- ggplot() +
  geom_sf(data = sf_darwinian_rmzero, aes(geometry = geometry, 
                                          fill = PD_deficit_all),
          color = "transparent", size = 0.1) +
  rcartocolor::scale_fill_carto_c(palette = "SunsetDark", 
                                  direction = 1, 
                                  limits = c(0, 1),  ## max percent overall
                                  breaks = seq(0, 1, by = .1),
                                  labels = glue::glue("{seq(0, 1, by = 0.1)}")) +
  guides(fill = guide_colorbar(barheight = unit(2.3, units = "mm"),  
                               barwidth = unit(100, units = "mm"),
                               direction = "horizontal",
                               ticks.colour = "grey20",
                               title.position = "top",
                               label.position = "bottom",
                               title.hjust = 0.5)) +
  labs(subtitle = "All species", fill = expression(PD[inserted]/PD[total])) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"),
        legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
        legend.text = element_text(family = "Times", color = "black", size = 10), 
        axis.text = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 11, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6))
        ) 


map_PD_deficit_perci <- ggplot() +
  geom_sf(data = sf_darwinian_rmzero, aes(geometry = geometry, 
                                          fill = deficitPerci),
          color = "transparent", size = 0.1) +
  labs(subtitle = "Perciformes") +
  rcartocolor::scale_fill_carto_c(palette = "SunsetDark", 
                                direction = 1, 
                                limits = c(0, 1),  ## max percent overall
                                breaks = seq(0, 1, by = .1),
                                labels = glue::glue("{seq(0, 1, by = 0.1)}")) +
  theme(legend.position = "none", panel.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.text = element_text(family = "Arial", color = "black", size = 10), 
        axis.text = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), 
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 11, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6))
        ) 

map_PD_deficit_chara <- ggplot() +
  geom_sf(data = sf_darwinian_rmzero, aes(geometry = geometry, 
                                          fill = deficitChar),
          color = "transparent", size = 0.1) +
  labs(subtitle = "Characiformes", fill = expression(PD[inserted]/PD[total])) +
  rcartocolor::scale_fill_carto_c(palette = "SunsetDark", 
                                  direction = 1, 
                                  limits = c(0, 1),  ## max percent overall
                                  breaks = seq(0, 1, by = .1),
                                  labels = glue::glue("{seq(0, 1, by = 0.1)}")) +
  theme(legend.position = "none", panel.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.text = element_text(family = "Arial", color = "black", size = 10), 
        axis.text = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), 
        plot.subtitle = element_text(family = "Arial",
                                     color = "black", 
                                     size = 11, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6))
  ) 

map_PD_deficit_cypri <- ggplot() +
  geom_sf(data = sf_darwinian_rmzero, aes(geometry = geometry, 
                                          fill = deficitCypri),
          color = "transparent", size = 0.1) +
  labs(subtitle = "Cypriniformes", fill = expression(PD[inserted]/PD[total])) +
  rcartocolor::scale_fill_carto_c(palette = "SunsetDark", 
                                  direction = 1, 
                                  limits = c(0, 1),  ## max percent overall
                                  breaks = seq(0, 1, by = .1),
                                  labels = glue::glue("{seq(0, 1, by = 0.1)}")) +
  theme(legend.position = "none", panel.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.text = element_text(family = "Arial", color = "black", size = 10), 
        axis.text = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), 
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black", 
                                     size = 11, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6))
  ) 

map_PD_deficit_siluri <- ggplot() +
  geom_sf(data = sf_darwinian_rmzero, aes(geometry = geometry, 
                                          fill = deficitSilu),
          color = "transparent", size = 0.1) +
  labs(subtitle = "Siluriformes", fill = expression(PD[inserted]/PD[total])) +
  rcartocolor::scale_fill_carto_c(palette = "SunsetDark", 
                                  direction = 1, 
                                  limits = c(0, 1),  ## max percent overall
                                  breaks = seq(0, 1, by = .1),
                                  labels = glue::glue("{seq(0, 1, by = 0.1)}")) +
  theme(legend.position = "none", panel.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.text = element_text(family = "Arial", color = "black", size = 10), 
        axis.text = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), 
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 11, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6))
  ) 

# full pannel maps ------------------------------------------------------------


full_map <- map_PD_deficit_wrld / (map_PD_deficit_chara | map_PD_deficit_cypri | map_PD_deficit_perci | map_PD_deficit_siluri) +
  plot_layout(heights = c(1, .75))

# saving plot -------------------------------------------------------------

ggsave(filename = here::here("output", "images", "world_map_DarwinianShort.png"), plot = full_map,
       width = 10, height = 5, 
       dpi = 500)
ggsave(here::here("output", "images", "world_map_DarwinianShort.pdf"), plot = full_map,
       width = 10, height = 5, device = cairo_pdf)

