
# read libraries ---------------------------------------------------------------

library(FishPhyloMaker)

# analyzing all data ------------------------------------------------------

species_list <- unique(occ_drainage$X6.Fishbase.Valid.Species.Name)

all_taxa_names <- FishTaxaMaker(data = species_list, allow.manual.insert = TRUE)
Cyprinidae
Loricariidae
Pseudopimelodidae
Petromyzontidae
Nemacheilidae
Loricariidae
Trichomycteridae
Procatopodidae

phylo_all_spp <- FishPhyloMaker(data = all_taxa_names$Taxon_data_FishPhyloMaker, 
                                insert.base.node = TRUE, 
                                return.insertions = TRUE)

phylo_all <- phylo_all_spp$Phylogeny
insertions <- phylo_all_spp$Insertions


# separating ecoregions ---------------------------------------------------

ecoregion <- names(table(occ_drainage$X3_Ecoregion))
list_ecoregions <- lapply(ecoregion, function(x) data[which(data$X3_Ecoregion == x), ])
names(list_ecoregions) <- ecoregion
species_list_ecoregions <- lapply(list_ecoregions, function(x) x$Genus.species)
names(species_list_ecoregions) <- ecoregion


# phylogenies per ecoregion ------------------------------------------------------


# saving results ----------------------------------------------------------
saveRDS(object = phylo_all, file = here::here("output", "phylo_all.rds"))
saveRDS(object = insertions, file = here::here("output", "insertions_all.rds"))

