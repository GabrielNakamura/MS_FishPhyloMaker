
# reading data and functions ----------------------------------------------
phylo_all <- readRDS(here::here("output", "phylo_all.rds"))
insertions_all <- readRDS(here::here("output", "insertions_all.rds"))
drainage_basins <- read.csv(here::here("data", "Drainage_Basins_Table.csv"), sep = ";")
occ_drainage <- read.csv(here::here("data", "Occurrence_Table.csv"), sep = ";")
occ_drainage$X2.Species.Name.in.Source <- gsub("[.]", "_", occ_drainage$X2.Species.Name.in.Source)
occ_drainage$X6.Fishbase.Valid.Species.Name <- gsub("[.]", "_", occ_drainage$X6.Fishbase.Valid.Species.Name)

# libraries
library(FishPhyloMaker)

basin_names <- unique(occ_drainage$X1.Basin.Name) # extracting basin names
list_spp_perBasin <- lapply(basin_names, 
                            function(x) occ_drainage[which(x == occ_drainage$X1.Basin.Name)
                                                     , "X6.Fishbase.Valid.Species.Name"]) # valid species names occurrence

names(list_spp_perBasin) <- unique(occ_drainage$X1.Basin.Name) # basin names

# prunning phylogeny for each basin
phy_perBasin <- lapply(list_spp_perBasin, function(x){
  rm_tip <- ape::drop.tip(phylo_all, tip = x)$tip.label
  ape::drop.tip(phylo_all, tip = rm_tip)
})

PD_deficit_all <- lapply(phy_perBasin, function(x) FishPhyloMaker::PD_deficit(phylo = x, 
                                                                              data = insertions_all, 
                                                                              level = c("Congeneric_insertion",
                                                                                        "Family_insertion")
                                                                              )
)

# deficit by orders -------------------------------------------------------

phy_perBasin <- lapply(list_spp_perBasin, function(x){
  rm_tip <- ape::drop.tip(phylo_all, tip = x)$tip.label
  ape::drop.tip(phylo_all, tip = rm_tip)
})

phy_perBasin_rmphy <- phy_perBasin[-match(names(unlist(lapply(phy_perBasin, function(x) which(length(x) == 0)))), 
                                          names(phy_perBasin))]
## siluriformes

phy_perBasin_Siluri <- lapply(phy_perBasin_rmphy, function(x) ape::drop.tip(x, 
                                                                            ape::drop.tip(x, 
                                                                                          tip = dplyr::filter(phylo_drainages$Insertions_data, 
                                                                                                              o == "Siluriformes")[, "s"])$tip.label)
                              )
PD_deficit_Siluri <- unlist(lapply(phy_perBasin_Siluri, function(x) FishPhyloMaker::PD_defict(phylo = x, 
                                                                     data = phylo_drainages$Insertions_data, 
                                                                     level = c("Congeneric_insertion",
                                                                               "Family_insertion"))
))

## Characiformes

phy_perBasin_Characi <- lapply(phy_perBasin_rmphy, function(x) ape::drop.tip(x, 
                                                                            ape::drop.tip(x, 
                                                                                          tip = dplyr::filter(phylo_drainages$Insertions_data, 
                                                                                                              o == "Characiformes")[, "s"])$tip.label)
)
PD_deficit_Characi <- unlist(lapply(phy_perBasin_Characi, function(x) FishPhyloMaker::PD_defict(phylo = x, 
                                                                                              data = phylo_drainages$Insertions_data, 
                                                                                              level = c("Congeneric_insertion",
                                                                                                        "Family_insertion"))
))

## Cipryniformes

phy_perBasin_Cyprini <- lapply(phy_perBasin_rmphy, function(x) ape::drop.tip(x, 
                                                                             ape::drop.tip(x, 
                                                                                           tip = dplyr::filter(phylo_drainages$Insertions_data, 
                                                                                                               o == "Cypriniformes")[, "s"])$tip.label)
)
PD_deficit_Cyprini <- unlist(lapply(phy_perBasin_Cyprini, function(x) FishPhyloMaker::PD_defict(phylo = x, 
                                                                                               data = phylo_drainages$Insertions_data, 
                                                                                               level = c("Congeneric_insertion",
                                                                                                         "Family_insertion"))
))

## Perciformes

phy_perBasin_Perci <- lapply(phy_perBasin_rmphy, function(x) ape::drop.tip(x, 
                                                                             ape::drop.tip(x, 
                                                                                           tip = dplyr::filter(phylo_drainages$Insertions_data, 
                                                                                                               o == "Perciformes")[, "s"])$tip.label)
)
PD_deficit_Perci <- unlist(lapply(phy_perBasin_Perci, function(x) FishPhyloMaker::PD_defict(phylo = x, 
                                                                                               data = phylo_drainages$Insertions_data, 
                                                                                               level = c("Congeneric_insertion",
                                                                                                         "Family_insertion"))
))

## data frame deficit per order

PD_deficit_perOrder <- data.frame(deficitChar = PD_deficit_Characi,
                               deficitPerci = PD_deficit_Perci, 
                               deficitCypri = PD_deficit_Cyprini, 
                               deficitSilu = PD_deficit_Siluri)

# saving objects ----------------------------------------------------------

saveRDS(object = phylo_drainages, file = here::here("output", "phylo_drainages.rds"))
saveRDS(object = PD_deficit_perOrder, file = here::here("output", "PD_deficit_perOrder.rds"))
saveRDS(PD_deficit_all, file = here::here("output", "PD_deficit_all.rds"))

