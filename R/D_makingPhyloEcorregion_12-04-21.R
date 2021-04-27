
# read data and libraries ---------------------------------------------------------------
devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main", force = TRUE)
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

taxa_afrotropics <- FishTaxaMaker(data = species_list_ecoregions$Afrotropic)
res_afrotropics <- FishPhyloMaker(data = taxa_afrotropics, return.insertions = TRUE)
Nematabramis
Cnesterodon
Cyprinidae
Heterochromis
Nematabramis Acapoeta
Heterochromis
Petrocephalus
Distichodontidae
Heterochromis
Cichlidae
Caprichromis
Distichodontidae
Channidae
Hemistichodus
Redigobius
Poeciliidae
Horabagrus
Cichlidae
Corematodus
Amphiliidae

Pseudopimelodidae
Siluriformes
Hiodontidae
Gonorynchiformes
saveRDS(object = res_afrotropics, file = here::here("phylo_afrotropics.rds"))


taxa_australasia <- FishTaxaMaker(data = species_list_ecoregions$Australasia)
res_australasia <- FishPhyloMaker(data = taxa_australasia, return.insertions = TRUE)
Tandanus
Pandaka
Soleidae
Pseudogobiopsis Synechogobius
Glyptauchen

Chirocentridae
Clupeiformes
Clupeiformes
saveRDS(object = res_australasia, file = here::here("phylo_australasia.rds"))


taxa_indomalay <- FishTaxaMaker(data = species_list_ecoregions$`Indo-Malay`)
res_indomalay <- FishPhyloMaker(data = taxa_indomalay, return.insertions = TRUE)
Soleidae
Sinogastromyzon
Cyprinidae
Genyonemus Seriphus
Soleidae
Parambassis
Syngnathidae
Pseudogobiopsis
Neosilurus
Platyallabes
Doryrhamphus Microphis
Siluridae
Cyprinidae
Pseudogobiopsis
Phallostethidae
Schilbeidae
Aspidoparia
Sillaginodes
Siluridae
Sagamia
Cyprinidae

Clupeiformes
Chirocentridae
Denticipitidae Chirocentridae
Xiphiidae
Perciformes
Clupeiformes
saveRDS(object = res_indomalay, file = here::here("phylo_indomalay.rds"))


taxa_neartic <- FishTaxaMaker(data = species_list_ecoregions$Nearctic)
res_neartic <- FishPhyloMaker(data = taxa_neartic, return.insertions = TRUE)
saveRDS(object = res_neartic, file = here::here("phylo_neartic.rds"))

taxa_neotropic <- FishTaxaMaker(data = species_list_ecoregions$Neotropic)
res_neotropic <- FishPhyloMaker(data = taxa_neotropic, return.insertions = TRUE)
Schizodon
Paravandellia
Ctenoluciidae
Compsura Serrapinnus
Curimatidae
Gymnotidae
Xenodexia
Lamontichthys
Characidae
Lamontichthys Fonchiiichthys
Anostomus
Anostomus
Plagioscion
Trichomycteridae
Lophiosilurus
Characidae
Fluviphylax
Trichomycteridae
Pinguipedidae
Anostomidae
Physopyxis
Hoplocharax
Anostomidae
Monopterus
Auchenipteridae
saveRDS(object = res_neotropic, file = here::here("phylo_neotropic.rds"))

taxa_oceania <- FishTaxaMaker(data = species_list_ecoregions$Oceania)
res_oceania <- FishPhyloMaker(data = taxa_oceania, return.insertions = TRUE)
saveRDS(res_oceania, file = here::here("phylo_oceania.rds"))

taxa_paleartic <- FishTaxaMaker(data = species_list_ecoregions$Palearctic)
res_paleartic <- FishPhyloMaker(data = taxa_paleartic, return.insertions = TRUE)
Cabdio
Securicula
Gobiidae
Aspidoparia

Clupeiformes
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

