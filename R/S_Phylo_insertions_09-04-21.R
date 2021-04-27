
# reading phylogenies -----------------------------------------------------
phylo_afrotropics <- readRDS(file = "phylo_afrotropics.rds")
phylo_australasia <- readRDS(file = "phylo_australasia.rds")
phylo_neartic <- readRDS(file = "phylo_neartic.rds")
phylo_indomalay <- readRDS(file = "phylo_indomalay.rds")
phylo_oceania <- readRDS(file = "phylo_oceania.rds")
phylo_paleartic <- readRDS(file = "phylo_paleartic.rds")
phylo_neotropic <- readRDS(file = "phylo_neotropic.rds")

# plotting figures --------------------------------------------------------
library(phytools)
library(ggtree)
library(ggplot2)
library(cowplot)
tree_afrotropics <- phylo_afrotropics$Phylogeny
insertions_afrotropics <- phylo_afrotropics$Insertions_data[-which("Present_in_Tree" == 
                                                           phylo_afrotropics$Insertions_data[, "insertions"]), ]
insertions_afrotropics <- phylo_afrotropics$Insertions_data

p.base_afrotropics <- ggtree(tree_afrotropics, layout = "circular", size = .3)  %<+% insertions_afrotropics +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) 


p.full_afrotropics <- p.base_afrotropics +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = T,
                        labels = c("Congeneric", "Congeneric F", "Family", "Order", "Presente in tree"))


tree_australasia <- phylo_australasia$Phylogeny
#insertions_australasia <- phylo_australasia$Insertions_data[-which("Present_in_Tree" == 
#                                                                     phylo_australasia$Insertions_data[, "insertions"]), ]
insertions_australasia <- phylo_australasia$Insertions_data

p.base_australasia <- ggtree(tree_australasia, layout = "circular", size = .3)  %<+% insertions_australasia +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) 


p.full_australasia <- p.base_australasia +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = T,
                        labels = c("Congeneric", "Family", "Order", "Presente in tree"))


tree_indomalay <- phylo_indomalay$Phylogeny
#insertions_australasia <- phylo_australasia$Insertions_data[-which("Present_in_Tree" == 
#                                                                     phylo_australasia$Insertions_data[, "insertions"]), ]
insertions_indomalay <- phylo_indomalay$Insertions_data

p.base_indomalay <- ggtree(tree_indomalay, layout = "circular", size = .3)  %<+% insertions_indomalay +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) 


p.full_indomalay <- p.base_indomalay +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = T,
                        labels = c("Congeneric", "Congeneric F", "Family", "Order", "Presente in tree"))

tree_neartic <- phylo_neartic$Phylogeny
#insertions_australasia <- phylo_australasia$Insertions_data[-which("Present_in_Tree" == 
#                                                                     phylo_australasia$Insertions_data[, "insertions"]), ]
insertions_neartic <- phylo_neartic$Insertions_data

p.base_neartic <- ggtree(tree_neartic, layout = "circular", size = .3)  %<+% insertions_neartic +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) 


p.full_neartic <- p.base_neartic +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = T,
                        labels = c("Congeneric", "Presente in tree"))


tree_paleartic <- phylo_paleartic$Phylogeny
#insertions_australasia <- phylo_australasia$Insertions_data[-which("Present_in_Tree" == 
#                                                                     phylo_australasia$Insertions_data[, "insertions"]), ]
insertions_paleartic <- phylo_paleartic$Insertions_data

p.base_paleartic <- ggtree(tree_paleartic, layout = "circular", size = .3)  %<+% insertions_paleartic +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) 


p.full_paleartic <- p.base_paleartic +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = T,
                        labels = c("Congeneric", "Family", "Not inserted", "Order", "Present in tree"))


tree_neotropic <- phylo_neotropic$Phylogeny
#insertions_australasia <- phylo_australasia$Insertions_data[-which("Present_in_Tree" == 
#                                                                     phylo_australasia$Insertions_data[, "insertions"]), ]
insertions_neotropic <- phylo_neotropic$Insertions_data

p.base_neotropic <- ggtree(tree_neotropic, layout = "circular", size = .3)  %<+% insertions_neotropic +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) 


p.full_neotropic <- p.base_neotropic +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = T,
                        labels = c("Congeneric", "Congeneric F", "Family", "Order", "Present in tree"))


quartz()
plot_grid(p.full_afrotropics, p.full_australasia, p.full_indomalay, p.full_neartic, p.full_neotropic, p.full_paleartic)

