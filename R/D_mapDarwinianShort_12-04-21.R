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

insertions_allDrainages <- phylo_drainages$Insertions_data
phylo_allDrainages <- phylo_drainages$Phylogeny
basin_names <- unique(occ_drainage$X1.Basin.Name)
list_spp_perBasin <- lapply(basin_names, function(x) occ_drainage[which(x == occ_drainage$X1.Basin.Name), "X6.Fishbase.Valid.Species.Name"])
names(list_spp_perBasin) <- unique(occ_drainage$X1.Basin.Name)
phy_perBasin <- lapply(list_spp_perBasin, function(x){
  rm_tip <- ape::drop.tip(phylo_drainages$Phylogeny, tip = x)$tip.label
  ape::drop.tip(phylo_drainages$Phylogeny, tip = rm_tip)
})

PD_deficit_all <- lapply(phy_perBasin, function(x) PD_defict(phylo = x, 
                                           data = insertions_allDrainages, 
                                           level = "Congeneric_insertion")
       )


