# read data ---------------------------------------------------------------
raw_data <- read.csv(here::here("data", "osm-raw-data.csv"), sep=";")
data <- raw_data[-which(duplicated(raw_data$SpecCode)==TRUE),]
data$Genus.species <- gsub("[.]","_",data$Genus.species)
species_list <- data$Genus.species

# read results ------------------------------------------------------------

res_afrotropics <- readRDS(here::here("output", "phylo_afrotropics.rds"))
res_indomalay <- readRDS(here::here("output", "phylo_indomalay.rds"))
res_neartic <- readRDS(here::here("output", "phylo_neartic.rds"))
res_neotropic <- readRDS(here::here("output", "phylo_neotropic.rds"))

res_afrotropics$Insertions_data$basin  <- "Afrotropic"
res_indomalay$Insertions_data$basin  <- "Indomalay"
res_neartic$Insertions_data$basin  <- "Neartic"
res_neotropic$Insertions_data$basin  <- "Neotropic"

new_data <- rbind(res_afrotropics$Insertions_data, res_indomalay$Insertions_data, 
                  res_neartic$Insertions_data, res_neotropic$Insertions_data)
new_data$insertions <- as.factor(new_data$insertions)
new_data$insertions2 <- factor(new_data$insertions, levels = c("Present_in_Tree", "Congeneric_insertion", "Congeneric_insertion_roundFamily", "Family_insertion", "Order_insertion"))
head(new_data)

# Stacked + percent
library(ggplot2)
ggplot(new_data, aes(x = insertions2, group = basin)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat = "count") +
  labs(y = "Percent", fill = "Insertion type") +
  facet_wrap(~basin) +
  scale_y_continuous(labels = scales::percent)  +
  scale_x_discrete(labels = c("Present in Tree", 
                              "Congeneric insertion",
                              "Family insertion",
                              "Congeneric at Family",
                              "Order insertion")) +
  rcartocolor::scale_fill_carto_d(palette = "SunsetDark", 
                                  direction = 1) +
  theme(legend.position = "none", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
