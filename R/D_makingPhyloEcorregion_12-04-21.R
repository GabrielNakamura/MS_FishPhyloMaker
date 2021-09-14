
# read libraries ---------------------------------------------------------------

library(FishPhyloMaker)


# analyzing all data ------------------------------------------------------

all_taxa_names <- FishTaxaMaker(data = species_list, allow.manual.insert = TRUE)

all_taxa_names <- FishTaxaMaker(data = species_list, allow.manual.insert = TRUE)

phylo_all_spp <- FishPhyloMaker(data = all_taxa_names$Taxon_data_FishPhyloMaker, 
                                insert.base.node = TRUE, 
                                return.insertions = TRUE)


phylo_all <- phylo_all_spp$Phylogeny

# separating ecoregions ---------------------------------------------------
ecoregion <- names(table(data$X3_Ecoregion))
list_ecoregions <- lapply(ecoregion, function(x) data[which(data$X3_Ecoregion == x), ])
names(list_ecoregions) <- ecoregion
species_list_ecoregions <- lapply(list_ecoregions, function(x) x$Genus.species)
names(species_list_ecoregions) <- ecoregion


# making phylogenies ------------------------------------------------------

taxa_afrotropics <- FishPhyloMaker::FishTaxaMaker(data = species_list_ecoregions$Afrotropic, allow.manual.insert = TRUE)
taxa_afrotropics$Taxon_data_FishPhyloMaker[, "o"][which(taxa_afrotropics$Taxon_data_FishPhyloMaker$f == "Cichlidae")] <- "Cichliformes"
res_afrotropics <- FishPhyloMaker(data = taxa_afrotropics$Taxon_data_FishPhyloMaker, insert.base.node = TRUE, return.insertions = TRUE)


taxa_australasia <- FishTaxaMaker(data = species_list_ecoregions$Australasia)
taxa_australasia$Taxon_data_FishPhyloMaker[, "o"][which(taxa_australasia$Taxon_data_FishPhyloMaker$o == "Scorpaeniformes")] <- "Perciformes"
res_australasia <- FishPhyloMaker(data = taxa_australasia$Taxon_data_FishPhyloMaker, 
                                  return.insertions = TRUE, 
                                  insert.base.node = TRUE)


taxa_indomalay <- FishTaxaMaker(data = species_list_ecoregions$`Indo-Malay`)
res_indomalay <- FishPhyloMaker(data = taxa_indomalay$Taxon_data_FishPhyloMaker, 
                                return.insertions = TRUE, 
                                insert.base.node = TRUE)


taxa_neartic <- FishTaxaMaker(data = species_list_ecoregions$Nearctic)
taxa_neartic$Taxon_data_FishPhyloMaker[which(taxa_neartic$Taxon_data_FishPhyloMaker$o == "Scorpaeniformes"), "o"] <- "Perciformes"
taxa_neartic$Taxon_data_FishPhyloMaker[which(taxa_neartic$Taxon_data_FishPhyloMaker$f == "Cichlidae"), "o"] <- "Cichliformes"
res_neartic <- FishPhyloMaker(data = taxa_neartic$Taxon_data_FishPhyloMaker, return.insertions = TRUE, insert.base.node = TRUE)

taxa_neotropic <- FishTaxaMaker(data = species_list_ecoregions$Neotropic)
taxa_neotropic$Taxon_data_FishPhyloMaker[which(taxa_neotropic$Taxon_data_FishPhyloMaker$f == "Cichlidae"), "o"] <- "Cichliformes"
res_neotropic <- FishPhyloMaker(data = taxa_neotropic$Taxon_data_FishPhyloMaker, 
                                return.insertions = TRUE, 
                                insert.base.node = TRUE)

taxa_oceania <- FishTaxaMaker(data = species_list_ecoregions$Oceania)
res_oceania <- FishPhyloMaker(data = taxa_oceania$Taxon_data_FishPhyloMaker, 
                              return.insertions = TRUE, 
                              insert.base.node = TRUE)

taxa_paleartic <- FishTaxaMaker(data = species_list_ecoregions$Palearctic)
taxa_paleartic$Taxon_data_FishPhyloMaker[which(taxa_paleartic$Taxon_data_FishPhyloMaker$f == "Cichlidae"), "o"] <- "Cichliformes"
taxa_paleartic$Taxon_data_FishPhyloMaker[which(taxa_paleartic$Taxon_data_FishPhyloMaker$o == "Scorpaeniformes"), "o"] <- "Perciformes"
res_paleartic <- FishPhyloMaker(data = taxa_paleartic$Taxon_data_FishPhyloMaker, 
                                return.insertions = TRUE, 
                                insert.base.node = TRUE)

# saving results ----------------------------------------------------------
saveRDS(object = phylo_all, file = here::here("output", "phylo_all.rds"))
saveRDS(object = res_afrotropics, file = here::here("output", "phylo_afrotropics.rds"))
saveRDS(res_paleartic, file = here::here("output", "phylo_paleartic.rds"))
saveRDS(res_oceania, file = here::here("output", "phylo_oceania.rds"))
saveRDS(object = res_neotropic, file = here::here("output", "phylo_neotropic.rds"))
saveRDS(object = res_neartic, file = here::here("output", "phylo_neartic.rds"))
saveRDS(object = res_neotropic, file = here::here("output", "phylo_neotropic.rds"))
saveRDS(object = res_indomalay, file = here::here("output", "phylo_indomalay.rds"))
saveRDS(object = res_australasia, file = here::here("output", "phylo_australasia.rds"))


