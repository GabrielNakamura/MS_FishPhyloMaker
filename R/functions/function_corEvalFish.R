correlation_eval <- 
  function(list_tree){
    complete_phylo <- fishtree::fishtree_phylogeny()
    tree_simul_rm <- ape::drop.tip(phy = list_tree$Phylogeny, tip = list_tree$Insertions_data$s)
    tree_simul <- ape::drop.tip(phy = list_tree$Phylogeny,
                                tip = tree_simul_rm$tip.label)
    spp_sing <- list_tree$samp_spp[which(unlist(lapply(strsplit(list_tree$samp_spp, split = "*._"), function(x) x[2]))
                             == "singleton")]
    spp_substitute_simul <- complete_phylo$tip.label[match(unlist(strsplit(spp_sing, split = "_.*")), 
                                                           unlist(strsplit(complete_phylo$tip.label, split = "_.*")))]
    list_tree$samp_spp[match(spp_sing, list_tree$samp_spp)] <- spp_substitute_simul
    tree_original_rm <- ape::drop.tip(phy = complete_phylo, tip = list_tree$samp_spp)
    tree_original <- ape::drop.tip(phy = complete_phylo, tip = tree_original_rm$tip.label)
    genus_rm <- unlist(strsplit(x = list_tree$Insertions_data$s[which(is.na(match(list_tree$Insertions_data$s, 
                                                                                  list_tree$Phylogeny$tip.label)) == TRUE)],
                                split = "_.*"))
    tip_rm_2 <- tree_original$tip.label[match(genus_rm, 
                                              unlist(strsplit(tree_original$tip.label, split = "_.*")))]
    phylo_obs <- ape::drop.tip(tree_original, tip = tip_rm_2)
    names_remove_simul <- unlist(strsplit(names(which(table(list_tree$samp_spp) > 1)), split = "_.*"))
    if(length(phylo_obs$tip.label) != length(tree_simul$tip.label)){
      genus_duplicated <- which(table(list_tree$samp_spp) > 1)
      names_remove_simul2 <- tree_simul$tip.label[match(unlist(strsplit(names(genus_duplicated), split = "_.*")), 
                                                        unlist(strsplit(tree_simul$tip.label, split = "_.*")))]
      tree_simul <- ape::drop.tip(phy = tree_simul, tip = names_remove_simul2)
      phylo_obs <- ape::drop.tip(phy = phylo_obs, tip = names(genus_duplicated))
    }
    cophenetic_simul <- cophenetic(tree_simul)
    cophenetic_original <- cophenetic(phylo_obs)
    distinctness_treeSimul <- phyloregion::evol_distinct(tree = tree_simul, type = c("equal.splits", "fair.proportion"))
    distinctness_treeObs <- phyloregion::evol_distinct(tree = phylo_obs, type = c("equal.splits", "fair.proportion"))
    distinctness_all <- data.frame(distinctness_obs = distinctness_treeObs,
                                   distinctness_simul = distinctness_treeSimul, 
                                   row.names = names(distinctness_treeObs)) 
    cor_distinctness <- cor(distinctness_treeSimul, distinctness_treeObs)
    cor_cophenetic <- cor(as.vector(cophenetic_original), as.vector(cophenetic_simul))
    correlations <- c(cor_distinctness, cor_cophenetic)
    names(correlations) <- c("cor_distinctness", "cor_cophenetic")
    list_res <- vector(mode = "list", length = 4)
    list_res$cophenetic_simuldist <- cophenetic_simul
    list_res$cophenetic_obsdist <- cophenetic_original
    list_res$distinctness <- distinctness_all
    list_res$correlation <- correlations
    return(list_res)
  }

