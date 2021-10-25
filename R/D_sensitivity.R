# sensitivity analysis 

# reading data and packages -----------------------------------------------

library(fishtree)
library(rfishbase)
library(parallel)
source(here::here("R", "functions", "function_evalPhyloMaker.R"))

# preparing backbone tree -------------------------------------------------

progress.bar <- TRUE
source(here::here("R","functions", "internal_filter_rank.R"))
source(here::here("R", "functions", "internal_treedata_modif.R"))

fishbasedata <- as.data.frame(data.frame(rfishbase::load_taxa()))
tree_complete <- fishtree::fishtree_phylogeny()
data <- FishPhyloMaker::FishTaxaMaker(data = tree_complete$tip.label, allow.manual.insert = FALSE)
data <- data$All_info_fishbase[, c("user_spp", "Family", "Order")]
colnames(data) <- c("s", "f", "o")
data$o <- gsub("/.*", 
               replacement = "", 
               data$o)
rank_order <- as.character(unique(data$o))
rank_family <- as.character(unique(data$f))
spp <- as.character(data$s)
all_families <- unique(unlist(lapply(rank_order, function(x) {
  fishbasedata[which(x == fishbasedata$Order), 5]
})))
families_in_orders <- suppressWarnings(all_families[which(unique(data$f) != 
                                                            all_families)])
families_order_and_data <- unique(c(rank_family, families_in_orders))
list_family <- vector(mode = "list", length = length(families_order_and_data))
if (progress.bar == TRUE) {
  pb_download_family <- progress::progress_bar$new(format = "  downloading families [:bar] :percent", 
                                                   total = length(families_order_and_data), clear = FALSE, 
                                                   width = 60, current = "<", incomplete = ">", 
                                                   complete = ">")
}
for (i in 1:length(families_order_and_data)) {
  list_family[[i]] <- tryCatch(paste(fishtree::fishtree_phylogeny(rank = families_order_and_data[i], 
                                                                  type = "chronogram_mrca")$tip.label), error = function(e) paste(families_order_and_data[i]))
  if (progress.bar == TRUE) {
    pb_download_family$tick()
  }
}
names(list_family) <- families_order_and_data
monotipic_family <- names(unlist(lapply(list_family, 
                                        function(x) which(length(x) == 1))))
list_monotipic <- vector(mode = "list", length = length(monotipic_family))
for (i in 1:length(monotipic_family)) {
  list_monotipic[[i]] <- tryCatch(fishtree::fishtree_taxonomy(rank = monotipic_family[i])[[1]]$taxonomy[[9]], 
                                  error = function(e) paste("not.found", "_", monotipic_family[i], 
                                                            sep = ""))
}
orders_to_add <- unique(unlist(list_monotipic[-which(sub("_.*", 
                                                         "", unlist(list_monotipic)) == "not.found")]))
differences_orders_toadd <- setdiff(rank_order, orders_to_add)
if (length(differences_orders_toadd) >= 1) {
  all_orders_include <- c(differences_orders_toadd, 
                          orders_to_add)
}
all_orders_include <- unique(c(rank_order, unique(orders_to_add)))
list_order <- vector(mode = "list", length = length(all_orders_include))
if (progress.bar == TRUE) {
  pb_download_order <- progress::progress_bar$new(format = "  downloading orders [:bar] :percent", 
                                                  total = length(all_orders_include), clear = FALSE, 
                                                  width = 60, current = "<", incomplete = ">", 
                                                  complete = ">")
}
for (i in 1:length(all_orders_include)) {
  list_order[[i]] <- tryCatch(paste(fishtree::fishtree_phylogeny(rank = all_orders_include[i], 
                                                                 type = "chronogram_mrca")$tip.label), error = function(e) paste(print(all_orders_include[i])))
  if (progress.bar == TRUE) {
    pb_download_order$tick()
  }
}
names(list_order) <- all_orders_include
phylo_order <- filter_rank(order = list_order)
phylo_order <- ape::makeNodeLabel(phy = phylo_order)
order_rm_list <- names(unlist(lapply(list_order, function(x) which(length(x) == 
                                                                     1))))
list_order <- list_order[-match(order_rm_list, names(list_order))]
list_non_monotipic <- list_family[setdiff(names(list_family), 
                                          monotipic_family)]
if (progress.bar == TRUE) {
  pb_names_family <- progress::progress_bar$new(format = "  Naming family nodes [:bar] :percent", 
                                                total = length(list_non_monotipic), clear = FALSE, 
                                                width = 60, current = "<", incomplete = ">", 
                                                complete = ">")
}
for (i in 1:length(list_non_monotipic)) {
  phylo_order <- ape::makeNodeLabel(phylo_order, "u", 
                                    nodeList = list(Fam_name = list_non_monotipic[[i]]))
  phylo_order$node.label[which(phylo_order$node.label == 
                                 "Fam_name")] <- paste(names(list_non_monotipic)[i])
  if (progress.bar == TRUE) {
    pb_names_family$tick()
  }
}
families_in_tree <- families_order_and_data[which(!is.na(match(families_order_and_data, 
                                                               phylo_order$node.label)) == T)]
families_monotipic_notfound <- setdiff(monotipic_family, 
                                       families_in_tree)
for (i in 1:length(families_monotipic_notfound)) {
  spp_tmp <- tryCatch(fishtree::fishtree_taxonomy(rank = families_monotipic_notfound[i])[[1]]$species, 
                      error = function(e) paste("not.found", "_", families_monotipic_notfound[i], 
                                                sep = ""))
  spp_tmp <- gsub("\\ ", "_", spp_tmp)
  list_family[which(families_monotipic_notfound[i] == 
                      names(list_family))] <- list(spp_tmp)
}
families_not_found_fishtree <- names(unlist(lapply(lapply(list_family, 
                                                          function(x) {
                                                            sub("_.*", "", x)
                                                          }), function(y) which(length(y) == 1) & which(y == 
                                                                                                          "not.found"))))
list_family_tobeaddnames <- list_family[-match(families_not_found_fishtree, 
                                               names(list_family))]
family_no_spp_in_tree <- names(unlist(lapply(lapply(list_family_tobeaddnames, 
                                                    function(x) {
                                                      sum(!is.na(match(x, phylo_order$tip.label)))
                                                    }), function(y) which(y == 0))))
list_family_tobeaddnames <- list_family_tobeaddnames[-match(family_no_spp_in_tree, 
                                                            names(list_family_tobeaddnames))]
if (length(list_family_tobeaddnames) > 0) {
  if (progress.bar == TRUE) {
    pb_names_family_toadd <- progress::progress_bar$new(format = "  Naming monotipic family nodes [:bar] :percent", 
                                                        total = length(list_family_tobeaddnames), clear = FALSE, 
                                                        width = 60, current = "<", incomplete = ">", 
                                                        complete = ">")
  }
  phylo_order <- phytools::force.ultrametric(phylo_order)
  for (i in 1:length(list_family_tobeaddnames)) {
    na_check <- sum(!is.na(match(list_family_tobeaddnames[[i]], 
                                 phylo_order$tip.label)))
    if (na_check == 1) {
      spp_singleton <- unlist(list(list_family_tobeaddnames[[i]][!is.na(match(list_family_tobeaddnames[[i]], 
                                                                              phylo_order$tip.label))]))
      spp_singleton_add <- paste(sub("_.*", "", spp_singleton), 
                                 "_", "singleton", sep = "")
      phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                    species = spp_singleton_add)
      list_family_tobeaddnames[i] <- list(c(spp_singleton, 
                                            spp_singleton_add))
    }
    phylo_order <- ape::makeNodeLabel(phylo_order, 
                                      "u", nodeList = list(Fam_name = list_family_tobeaddnames[[i]]))
    phylo_order$node.label[which(phylo_order$node.label == 
                                   "Fam_name")] <- paste(names(list_family_tobeaddnames)[i])
    if (progress.bar == TRUE) {
      pb_names_family_toadd$tick()
    }
  }
}


# sensitivity analysis ----------------------------------------------------

n.cluster <- parallel::detectCores()
parallel <- parallel::makeCluster(n.cluster, type = "PSOCK")
res_01 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                    probs = 0.1, 
                    tree_complete = phylo_order,
                    insert.base.node = TRUE, 
                    return.insertions = TRUE, 
                    progress.bar = TRUE)
res_015 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                              probs = 0.15, 
                              tree_complete = phylo_order,
                              insert.base.node = TRUE, 
                              return.insertions = TRUE, 
                              progress.bar = TRUE)
res_02 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                              probs = 0.2, 
                              tree_complete = phylo_order,
                              insert.base.node = TRUE, 
                              return.insertions = TRUE, 
                              progress.bar = TRUE)
res_025 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                              probs = 0.25, 
                              tree_complete = phylo_order,
                              insert.base.node = TRUE, 
                              return.insertions = TRUE, 
                              progress.bar = TRUE)
res_03 <- parallel::parLapply(parallel, 1:50, fun = eval_PhyloMaker,
                              probs = 0.3, 
                              tree_complete = phylo_order,
                              insert.base.node = TRUE, 
                              return.insertions = TRUE, 
                              progress.bar = TRUE)
parallel::stopCluster(parallel)

saveRDS(object = res_01, file = here::here("output", "res_010_50runs.rds"))
saveRDS(object = res_015, file = here::here("output", "res_015_50runs.rds"))
saveRDS(object = res_02, file = here::here("output", "res_020_50runs.rds"))
saveRDS(object = res_025, file = here::here("output", "res_025_50runs.rds"))
saveRDS(object = res_030, file = here::here("output", "res_030_50runs.rds"))
