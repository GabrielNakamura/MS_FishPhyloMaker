res_01 <- readRDS(here::here("output", "res_010_20runs.rds"))
res_02 <- readRDS(here::here("output", "res_020_50runs.rds"))
res_03 <- readRDS(here::here("output", "res_03_50runs.rds"))


complete_phylo <- fishtree::fishtree_phylogeny()
prunedtree1 <- ape::drop.tip(phy = list_res$Phylogeny, tip = list_res$samp_spp)
tree_rm <- ape::drop.tip(prunedtree1, tip = list_res$Insertions_data$s)
tree_simul <- ape::drop.tip(phy = prunedtree1, tip = tree_rm$tip.label)
tree_rm_original <- ape::drop.tip(phy = complete_phylo, tip = list_res$samp_spp)
tree_original <- ape::drop.tip(phy = complete_phylo, tip = tree_rm_original$tip.label)
genus_rm <- unlist(strsplit(x = list_res$Insertions_data$s[which(is.na(match(list_res$Insertions_data$s, list_res$Phylogeny$tip.label)) == TRUE)],
                            split = "_.*"))
tip_rm_2 <- tree_original$tip.label[match(genus_rm, 
                                           unlist(strsplit(tree_original$tip.label, split = "_.*")))]

phylo_obs <- ape::drop.tip(tree_original, tip = tip_rm_2)


cophenetic_simul <- cophenetic(tree_simul)
cophenetic_original <- cophenetic(phylo_obs)
vegan::mantel(cophenetic_original, cophenetic_simul, permutations = 0)

