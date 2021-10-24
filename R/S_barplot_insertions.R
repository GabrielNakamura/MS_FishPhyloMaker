
# reading data and libraries ----------------------------------------------

phylo_complete <- readRDS(here::here("output", "phylo_all.rds"))
insertions_all <- readRDS(here::here("output", "insertions_all.rds"))
list_spp_perBasin <- readRDS(here::here("output", "list_spp_perBasin.rds"))

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(patchwork)
library(ggplot2)

# setting theme -----------------------------------------------------------

theme_set(theme_bw())

# barplot -------------------------------------------------------------

# organizing data
Eco_basin <- lapply(1:length(list_spp_perBasin), function(x){
  drainage_basins[which(drainage_basins$X1.Basin.Name == names(list_spp_perBasin[x])), "X3.Ecoregion"]
})


# corrigir nomes para nomes válidos
data_basin_ecoregion <- data.frame(spp = unlist(list_spp_perBasin), 
                                   basin.name = names(unlist(list_spp_perBasin)),
                                   Ecoregion = unlist(lapply(1:length(Eco_basin), 
                                                             function(x) rep(Eco_basin[[x]],
                                                                             times = length(list_spp_perBasin[[x]])
                                                             )
                                   )
                                   )
)




data_basin_ecoregion$insertions <- insertions_all[match(data_basin_ecoregion$spp,insertions_all$s), "insertions"]
data_basin_ecoregion$insertions <- factor(data_basin_ecoregion$insertions, 
                                          levels = c("Present_in_Tree", 
                                                     "Congeneric_insertion",
                                                     "Family_insertion",
                                                     "Congeneric_Family_level", 
                                                     "Order_insertion", 
                                                     "Not_inserted"))
data_basin_ecoregion <- data_basin_ecoregion[-which(is.na(data_basin_ecoregion$insertions) == TRUE), ]


unique(data_basin_ecoregion[which(data_basin_ecoregion$insertions != "Present_in_Tree" & data_basin_ecoregion$insertions != "Not_inserted" & data_basin_ecoregion$Ecoregion == "Nearctic"), "spp"])
length(unique(data_basin_ecoregion[which(data_basin_ecoregion$insertions == "Congeneric_insertion" & data_basin_ecoregion$Ecoregion == "Nearctic"), "spp"]))

# Stacked + percent
barplot_insertion <- ggplot(data = data_basin_ecoregion, aes(x = Ecoregion, fill = insertions)) +
  geom_bar(na.rm = TRUE) +
  rcartocolor::scale_fill_carto_d(palette = "Safe", 
                                  labels = c("Present",
                                             "Congeneric",
                                             "Family", 
                                             "Congeneric Family",
                                             "Order", 
                                             "Not inserted")) +
  labs(y = "Total number of insertions", fill = "Type of insertion") +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "mm"),
        legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
        legend.text = element_text(family = "Times", color = "black", size = 12), 
        axis.text = element_text(family = "Times", color = "black", size = 12), panel.border = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 9, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6)
        )
  )
barplot_insertion

perc.insertions <- data_basin_ecoregion %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::count(insertions)

barplot_insertion_perc <- ggplot2::ggplot(data = perc.insertions, aes(y = n, x = Ecoregion)) +
  # coord_flip() +
  geom_bar(aes(fill = insertions), 
           position = position_fill(reverse = TRUE),
           stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.8) +
  labs(x = NULL, y = NULL) +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 10)) + # Adjust label text
  guides(fill = guide_legend(reverse = TRUE)) +
  rcartocolor::scale_fill_carto_d(palette = "Safe", 
                                  labels = c("Present",
                                             "Congeneric",
                                             "Family", 
                                             "Congeneric Family",
                                             "Order", 
                                             "Not inserted")) +
  labs(y = "Proportion of insertions", fill = "Type of insertion") +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "mm"),
        legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
        legend.text = element_text(family = "Times", color = "black", size = 12), 
        axis.text = element_text(family = "Times", color = "black", size = 12), panel.border = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 9, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6)
        )
  )
barplot_insertion_perc

# plotting tree -----------------------------------------------------------

df_insertion <- insertions_all[- which(insertions_all$insertions == "Not_inserted"), ]
df_insertion_org <- df_insertion[match(df_insertion$s, phylo_complete$tip.label), ]
df_insertion_tree <- data.frame(df_insertion_org$insertions)
rownames(df_insertion_tree) <- df_insertion_org$s
orders <- names(sort(table(insertions_all$o), decreasing = TRUE)[1:7])
nodedf <- data.frame(nodes = unlist(lapply(orders, 
                                           function(x) phytools::findMRCA(tree = phylo_complete, 
                                                                          tips = insertions_all[which(insertions_all$o == x), "s"])
)
))
names_df_order <- names(sort(table(insertions_all$o), decreasing = TRUE)[1:7])
nodedf$Orders <- names_df_order

# adjusting the order of insertions
insertions_all$insertions <- factor(insertions_all$insertions, 
                                   levels = c("Present_in_Tree", 
                                              "Congeneric_insertion",
                                              "Family_insertion",
                                              "Congeneric_Family_level", 
                                              "Order_insertion", 
                                              "Not_inserted")
                                   )



# plotting tree -----------------------------------------------------------

phylo <- ggtree::ggtree(phylo_complete, layout = "circular") + 
  geom_hilight(data = nodedf, mapping = aes(node = nodes), extendto = 400,
               alpha = 0.3, fill = "grey", color = "grey60",
               size = 0.05) +
  geom_cladelab(data = nodedf, 
                mapping=aes(node = nodes, 
                            label = Orders),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=FALSE, 
                fontsize=2.4,
                fontface="italic"
  ) +
  geom_fruit(data = insertions_all, geom = geom_tile,
             mapping = aes(y = s, x = insertions, fill = insertions),
             color = "grey50", offset = 0.1, size = 0.02, 
             pwidth = 0.3, stat = "identity") +
  rcartocolor::scale_fill_carto_d(palette = "Safe", 
                                   labels = c("Present",
                                              "Congeneric",
                                              "Family", 
                                              "Congeneric Family",
                                              "Order")) +
  labs(subtitle = "", fill = "Insertions") +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"),
        legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
        legend.text = element_text(family = "Times", color = "black", size = 12), 
        axis.text = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 9, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6)
                                     )
        )



# saving figures ----------------------------------------------------------

ggsave(filename = here::here("output", "images", "phylogeny_complete.png"), plot = phylo,
       width = 7, height = 8, 
       dpi = 500)
ggsave(here::here("output", "images", "barplot_insertions_perc.png"), plot = barplot_insertion_perc,
       width = 7, height = 7, dpi = 500)

ggsave(here::here("output", "images", "barplot_insertions.png"), plot = barplot_insertion,
       width = 7, height = 7, dpi = 500)



