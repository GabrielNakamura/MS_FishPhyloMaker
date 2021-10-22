# sensitivity analysis 

# reading data and packages -----------------------------------------------

library(fishtree)
library(rfishbase)
library(parallel)
source(here::here("R", "functions", "function_evalPhyloMaker.R"))

# percentage of insertions per ecoregion ----------------------------------

unlist(lapply(lapply(lapply(unique(data_basin_ecoregion$Ecoregion), 
                            function(x) which(data_basin_ecoregion$Ecoregion == x &
                                                data_basin_ecoregion$insertions != "Present_in_Tree")), function(y){
                                                  unique(data_basin_ecoregion[y,"spp"])
                                                }), function(z) length(z)/length(phylo_complete$tip.label)))


# sensitivity analysis ----------------------------------------------------

n.cluster <- parallel::detectCores()
parallel <- parallel::makeCluster(n.cluster, type = "PSOCK")
res_01 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                    probs = 0.1, 
                    insert.base.node = TRUE, 
                    return.insertions = TRUE, 
                    progress.bar = TRUE)
parallel::stopCluster(parallel)
res_02 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                              probs = 0.2, 
                              insert.base.node = TRUE, 
                              return.insertions = TRUE, 
                              progress.bar = TRUE)
parallel::stopCluster(parallel)
res_03 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                              probs = 0.3, 
                              insert.base.node = TRUE, 
                              return.insertions = TRUE, 
                              progress.bar = TRUE)
parallel::stopCluster(parallel)

saveRDS(object = res_01, file = here::here("output", "res_010_20runs.rds"))
