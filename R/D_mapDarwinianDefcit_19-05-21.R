
# reading data and functions ----------------------------------------------
# data
osm_raw_data <- read.csv(here::here("data", "osm-raw-data.csv"), sep = ";")
drainage_basins <- read.csv(here::here("data", "Drainage_Basins_Table.csv"), sep = ";")
occ_drainage <- read.csv(here::here("data", "Occurrence_Table.csv"), sep = ";")
occ_drainage$X2.Species.Name.in.Source <- gsub("[.]", "_", occ_drainage$X2.Species.Name.in.Source)
occ_drainage$X6.Fishbase.Valid.Species.Name <- gsub("[.]", "_", occ_drainage$X6.Fishbase.Valid.Species.Name)

# function
source(here::here("R", "functions", "F_darwinianDeficit.R"))

# libraries
library(FishPhyloMaker)

names_spp_valid <- unique(occ_drainage$X6.Fishbase.Valid.Species.Name)

# obtaining phylogeny  ----------------------------------------------------
taxon_data <- FishPhyloMaker::FishTaxaMaker(data = names_spp_valid, allow.manual.insert = TRUE)
Pimelodidae
Siluriformes 
Petromyzontidae
Petromyzontiformes
Nemacheilidae
Cypriniformes
Trichomycteridae
Siluriformes
Procatopodidae
Cyprinodontiformes


# correcting families and orders ------------------------------------------
FishPhyloMaker_data <- taxon_data$Taxon_data_FishPhyloMaker
FishPhyloMaker_data <- FishPhyloMaker_data[- which(is.na(unique(FishPhyloMaker_data$s)) == TRUE), ]
FishPhyloMaker_data[which(FishPhyloMaker_data$f == "Cichlidae"), "o"] <- "Cichliformes"
FishPhyloMaker_data[which(FishPhyloMaker_data$o == "Scorpaeniformes"), "o"] <- "Perciformes"

system.time(phylo_drainages <- FishPhyloMaker(data = FishPhyloMaker_data,
                                              insert.base.node = TRUE,
                                              return.insertions = TRUE,
                                              progress.bar = TRUE))
saveRDS(object = phylo_drainages, file = here::here("output", "phylo_drainages.rds"))


PD_deficit_all <- lapply(phy_perBasin, function(x) PD_defict(phylo = x, 
                                                             data = phylo_drainages$Insertions_data, 
                                                             level = "Congeneric_insertion")
)

saveRDS(PD_deficit_all, file = here::here("output", "PD_deficit_all.rds"))

