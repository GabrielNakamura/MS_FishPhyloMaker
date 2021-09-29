# sensibility analysis 


# reading data and packages -----------------------------------------------

library(fishtree)
library(FishPhyloMaker)

library(fishtree)
tree_complete <- fishtree_phylogeny()
probs <- seq(0.1, 0.7, 0.05)
samp <- ceiling(length(tree_complete$tip.label) * probs)
spp_samp <- lapply(probs, 
                   function(x) tree_complete$tip.label[sample(1:length(tree_complete$tip.label),
                                                              size = length(tree_complete$tip.label)*x)])
simul_unknown_tree <- lapply(spp_samp, function(x) ape::drop.tip(phy = tree_complete, tip = x))
names(simul_unknown_tree) <- paste("samp", probs, sep = "_")

lapply(simul_unk)
simul_unknown_spp <- lapply(spp_samp, function(x) {
  gsub("_.*", "_notFound", x)
}
)
