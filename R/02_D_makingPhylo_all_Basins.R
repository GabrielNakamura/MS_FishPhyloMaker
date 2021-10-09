
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


# manual inspection of taxa table -----------------------------------------

all_taxa_names$Taxon_data_FishPhyloMaker$o <- gsub("/.*", 
                                                   replacement = "", 
                                                   all_taxa_names$Taxon_data_FishPhyloMaker$o) # removing slash in Perciformes

phylo_all_spp <- FishPhyloMaker(data = all_taxa_names$Taxon_data_FishPhyloMaker, 
                                insert.base.node = TRUE, 
                                return.insertions = TRUE)

phylo_all <- phylo_all_spp$Phylogeny
insertions <- phylo_all_spp$Insertions


# saving results ----------------------------------------------------------
saveRDS(object = phylo_all, file = here::here("output", "phylo_all.rds"))
saveRDS(object = insertions, file = here::here("output", "insertions_all.rds"))
saveRDS(object = all_taxa_names, file = here::here("output", "all_taxa_names.rds"))

# Exploratory analysis of insertions -----------------------------------------------------------------

length(which(insertions$insertions != c("Present_in_Tree", "Not_inserted"))) # species inserted in phylogeny
length(which(insertions$insertions == "Present_in_Tree"))
length(which(insertions$insertions ==  "Not_inserted")) # species not inserted in the phylogenetic tree
length(which(insertions$insertions ==  "Not_inserted"))/dim(insertions)[1]
sum(is.na(match(all_taxa_names$All_info_fishbase$user_spp,
        gsub(" ", "_", all_taxa_names$All_info_fishbase$valid_names), nomatch = NA))) # total number of species not inserted in the phylogenetic tree

sum(is.na(match(all_taxa_names$All_info_fishbase$user_spp,
                gsub(" ", "_", all_taxa_names$All_info_fishbase$valid_names), nomatch = NA))) - length(all_taxa_names$All_info_fishbase$user_spp) # number of species inserted 
