# sensitivity analysis 

# reading data and packages -----------------------------------------------

library(fishtree)
library(rfishbase)
source(here::here("R", "functions", "function_evalPhyloMaker.R"))

# percentage of insertions per ecoregion ----------------------------------

unlist(lapply(lapply(lapply(unique(data_basin_ecoregion$Ecoregion), 
                            function(x) which(data_basin_ecoregion$Ecoregion == x &
                                                data_basin_ecoregion$insertions != "Present_in_Tree")), function(y){
                                                  unique(data_basin_ecoregion[y,"spp"])
                                                }), function(z) length(z)/length(phylo_complete$tip.label)))


# sensitivity analysis ----------------------------------------------------


