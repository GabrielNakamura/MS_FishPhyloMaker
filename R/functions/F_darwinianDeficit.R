PD_defict <- function(phylo, data, level = "Congeneric_insertion"){
  names_exclude <- phylo$tip.label[na.omit(match(data[which(data$insertions == "Present_in_Tree"), "s"], 
                                phylo$tip.label))]
  
  if(all(!is.na(match( phylo$tip.label, names_exclude))) == TRUE){
    PD_present <- sum(phylo$edge.length)
    Darwinian_defict <- 0
  } else{
    exclude <- ape::drop.tip(phylo, tip = names_exclude)$tip.label
    phylo_present <- ape::drop.tip(phylo, tip = exclude)
    PD_present <- sum(phylo_present$edge.length)
    if(length(level) == 1){
      level_exclude <- ape::drop.tip(phylo, tip = phylo$tip.label[na.omit(match(data[which(data$insertions == level), "s"], 
                                                                                phylo$tip.label))])$tip.label
    }
    if(length(level) > 1){
      level_exclude <- ape::drop.tip(phylo, data$insertions[!is.na(match(data$insertions, level)), "s"])$tip.label 
    }
    phylo_level <- ape::drop.tip(phylo, level_exclude)
    PD_level <- sum(phylo_level$edge.length)
    PD_total <- sum(PD_present, PD_level)
    Darwinian_defict <- PD_level/PD_total
  }
  return(Darwinian_defict)
}
