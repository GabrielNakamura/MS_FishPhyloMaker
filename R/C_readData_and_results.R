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

