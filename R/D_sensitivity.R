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

system.time(res_01 <- parallel::parLapply(parallel, 1:30, fun = eval_PhyloMaker,
                    probs = 0.1, 
                    insert.base.node = TRUE, 
                    return.insertions = TRUE, 
                    progress.bar = TRUE))
parallel::stopCluster(parallel)

res_0.1 <- lapply(rep(0.1, 50), function(x) eval_PhyloMaker(probs = x, 
                                                 insert.base.node = TRUE, 
                                                 return.insertions = TRUE, 
                                                 progress.bar = TRUE)
       )

res_0.2 <- lapply(rep(0.2, 50), function(x) eval_PhyloMaker(probs = x, 
                                                            insert.base.node = TRUE, 
                                                            return.insertions = TRUE, 
                                                            progress.bar = TRUE)
)

res_0.3 <- lapply(rep(0.3, 50), function(x) eval_PhyloMaker(probs = x, 
                                                            insert.base.node = TRUE, 
                                                            return.insertions = TRUE, 
                                                            progress.bar = TRUE)
)

saveRDS(object = res_0.1, file = here::here("output", "res_010_50runs.rds"))
