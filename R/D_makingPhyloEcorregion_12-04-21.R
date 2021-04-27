
# read data and libraries ---------------------------------------------------------------
library(FishPhyloMaker)
raw_data <- read.csv(here::here("data", "osm-raw-data.csv"), sep=";")
data <- raw_data[-which(duplicated(raw_data$SpecCode)==TRUE),]
data$Genus.species <- gsub("[.]","_",data$Genus.species)
species_list <- data$Genus.species


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

# calculating PD ----------------------------------------------------------

library(picante)

exclude <- drop.tip(phylo_afrotropics$Phylogeny, tip = match(phylo_afrotropics$Insertions_data[which(phylo_afrotropics$Insertions_data$insertions ==
                                                "Present_in_Tree"), "s"],
      phylo_afrotropics$Phylogeny$tip.label
      ))$tip.label
phylo_afrotropics_present <- drop.tip(phylo_afrotropics$Phylogeny, tip = exclude)
PD_present_afrotropic <- sum(phylo_afrotropics_present$edge.length)

congeneric_exclude<- drop.tip(phylo_afrotropics$Phylogeny, tip = match(phylo_afrotropics$Insertions_data[which(phylo_afrotropics$Insertions_data$insertions ==
                                                                                            "Congeneric_insertion"), "s"],
                                                  phylo_afrotropics$Phylogeny$tip.label
))$tip.label
phylo_afrotropics_congeneric <- drop.tip(phylo_afrotropics$Phylogeny, congeneric_exclude)
PD_congeneric_afrotropic <- sum(phylo_afrotropics_congeneric$edge.length)
PD_total_afrotropic <- sum(PD_present_afrotropic, PD_congeneric_afrotropic)
PD_congeneric_afrotropic/PD_total_afrotropic
PD_present_afrotropic/PD_total_afrotropic


defict_neotropic <- PD_defict(phylo = phylo_neotropic$Phylogeny, data = phylo_neotropic$Insertions_data, level = "Congeneric_insertion")
defict_neartic <- PD_defict(phylo = phylo_neartic$Phylogeny, data = phylo_neartic$Insertions_data, level = "Congeneric_insertion")
defict_afrotropic <- PD_defict(phylo = phylo_afrotropics$Phylogeny, data = phylo_afrotropics$Insertions_data, level = "Congeneric_insertion")
defict_indomalay <- PD_defict(phylo = phylo_indomalay$Phylogeny, data = phylo_indomalay$Insertions_data, level = "Congeneric_insertion")

