
# read libraries ---------------------------------------------------------------
library(FishPhyloMaker)


# analyzing all data ------------------------------------------------------

taxa_all <- FishPhyloMaker::FishTaxaMaker(data = species_list, allow.manual.insert = TRUE)
phylo_all <- FishPhyloMaker::FishPhyloMaker(data = taxa_all$Taxon_data_FishPhyloMaker, 
                                            insert.base.node = TRUE, 
                                            return.insertions = TRUE,
                                            progress.bar = TRUE)

# separating ecoregions ---------------------------------------------------
ecoregion <- names(table(data$X3_Ecoregion))
list_ecoregions <- lapply(ecoregion, function(x) data[which(data$X3_Ecoregion == x), ])
names(list_ecoregions) <- ecoregion
species_list_ecoregions <- lapply(list_ecoregions, function(x) x$Genus.species)
names(species_list_ecoregions) <- ecoregion

taxa_afrotropics <- FishPhyloMaker::FishTaxaMaker(data = species_list_ecoregions$Afrotropic, allow.manual.insert = TRUE)
taxa_afrotropics$Taxon_data_FishPhyloMaker[, "o"][which(taxa_afrotropics$Taxon_data_FishPhyloMaker$f == "Cichlidae")] <- "Cichliformes"
res_afrotropics <- FishPhyloMaker(data = taxa_afrotropics$Taxon_data_FishPhyloMaker, insert.base.node = TRUE, return.insertions = TRUE)
saveRDS(object = res_afrotropics, file = here::here("output", "phylo_afrotropics.rds"))


taxa_australasia <- FishTaxaMaker(data = species_list_ecoregions$Australasia)
taxa_australasia$Taxon_data_FishPhyloMaker[, "o"][which(taxa_australasia$Taxon_data_FishPhyloMaker$o == "Scorpaeniformes")] <- "Perciformes"
res_australasia <- FishPhyloMaker(data = taxa_australasia$Taxon_data_FishPhyloMaker, 
                                  return.insertions = TRUE, 
                                  insert.base.node = TRUE)
saveRDS(object = res_australasia, file = here::here("output", "phylo_australasia.rds"))


taxa_indomalay <- FishTaxaMaker(data = species_list_ecoregions$`Indo-Malay`)
res_indomalay <- FishPhyloMaker(data = taxa_indomalay$Taxon_data_FishPhyloMaker, 
                                return.insertions = TRUE, 
                                insert.base.node = TRUE)
saveRDS(object = res_indomalay, file = here::here("output", "phylo_indomalay.rds"))


taxa_neartic <- FishTaxaMaker(data = species_list_ecoregions$Nearctic)
taxa_neartic$Taxon_data_FishPhyloMaker[which(taxa_neartic$Taxon_data_FishPhyloMaker$o == "Scorpaeniformes"), "o"] <- "Perciformes"
taxa_neartic$Taxon_data_FishPhyloMaker[which(taxa_neartic$Taxon_data_FishPhyloMaker$f == "Cichlidae"), "o"] <- "Cichliformes"
res_neartic <- FishPhyloMaker(data = taxa_neartic$Taxon_data_FishPhyloMaker, return.insertions = TRUE, insert.base.node = TRUE)
saveRDS(object = res_neartic, file = here::here("output", "phylo_neartic.rds"))

taxa_neotropic <- FishTaxaMaker(data = species_list_ecoregions$Neotropic)
taxa_neotropic$Taxon_data_FishPhyloMaker[which(taxa_neotropic$Taxon_data_FishPhyloMaker$f == "Cichlidae"), "o"] <- "Cichliformes"
res_neotropic <- FishPhyloMaker(data = taxa_neotropic$Taxon_data_FishPhyloMaker, 
                                return.insertions = TRUE, 
                                insert.base.node = TRUE)
saveRDS(object = res_neotropic, file = here::here("phylo_neotropic.rds"))

taxa_oceania <- FishTaxaMaker(data = species_list_ecoregions$Oceania)
res_oceania <- FishPhyloMaker(data = taxa_oceania$Taxon_data_FishPhyloMaker, 
                              return.insertions = TRUE, 
                              insert.base.node = TRUE)
saveRDS(res_oceania, file = here::here("phylo_oceania.rds"))

taxa_paleartic <- FishTaxaMaker(data = species_list_ecoregions$Palearctic)
taxa_paleartic$Taxon_data_FishPhyloMaker[which(taxa_paleartic$Taxon_data_FishPhyloMaker$f == "Cichlidae"), "o"] <- "Cichliformes"
taxa_paleartic$Taxon_data_FishPhyloMaker[which(taxa_paleartic$Taxon_data_FishPhyloMaker$o == "Scorpaeniformes"), "o"] <- "Perciformes"
res_paleartic <- FishPhyloMaker(data = taxa_paleartic$Taxon_data_FishPhyloMaker, 
                                return.insertions = TRUE, 
                                insert.base.node = TRUE)
saveRDS(res_paleartic, file = here::here("phylo_paleartic.rds"))

