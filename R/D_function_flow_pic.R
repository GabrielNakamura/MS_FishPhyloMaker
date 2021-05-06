# Script to genetare the figures that represent the step-by-step of the 
# FishPhyloMaker funtion

library(phytools)
library(geiger)
library(ggtree)
library(ggimage)
library(cowplot)
library(viridis)

## TREE gnerated by the PhyloFishMaker function
tree.full <- ape::read.tree(here::here("data", "spp_example.txt"))

# finding the phylopic images for family nodes
fam.pic <- c("Cichlidae", "Loricariidae", 
             "Characidae")

f.pic <- ggimage::phylopic_uid(fam.pic)

#finding the node position for plot the phylopics
phytools::plotTree(tree = tree.full, node.numbers = T)
# 19 = Characidae
# 15 = Loricariidae
# 12 = Cichlidae

# create a data frame with the informations
phylopic_info <- data.frame(node = c(12, 15, 19),
                            phylopic = c("335e6e09-38a2-4a56-9e63-335dfc9e722a",
                                         "ae1d3ca5-f6d9-49dd-a431-e1bea5a3bdf3",
                                         "7be4fe69-28df-4a1f-99a0-666aafba3950"),
                            family = c("Cichlidae", "Loricariidae",
                                       "Characidae"))

# create a list specifying the type of insertions this species are 
# not in the phylogenetic tree, and will be included
# CG = congeneric insertion
# FN = insertion on the family node
# CGF = congeneric insertion, after the insertion of a species in family round
# NG = specifying a genus to place the species that will be add
# BG = specifying two genus names to insert the species between

insertions.types <- list(CG = "Gymnotus_sp1",
                         FN = "Rineloricaria_sp1",
                         CGF = "Rineloricaria_sp2", 
                         NG = "Rineloricaria_sp1",
                         BG = "Rineloricaria_sp1")

tree.a <- groupOTU(tree.full, insertions.types) 
tree.a$tip.label <- gsub("_", " ", tree.a$tip.label)

# plot for the complete tree
plot.a <- ggtree(tree.a) %<+% phylopic_info + 
  geom_tiplab(size = 3) + labs(title = "Backbone tree") +
  theme(legend.position = c(0,1)) +
  geom_nodelab(aes(image = phylopic), geom = "phylopic",
               color = 'gray40', size = .1) +
  geom_tiplab(size = 3, show.legend = FALSE) +
  xlim(0, 370) + geom_treescale(x = 0, width = 20, fontsize = 3.5, offset = .15)

# initiate the insertion procedure:
# first round: congeneric
tree.insert <- add.species.to.genus(tree = tree.full, species = "Gymnotus_sp1")

tree.b <- groupOTU(tree.insert, insertions.types)
tree.b$tip.label <- gsub("_", " ", tree.b$tip.label)

plot.b <- ggtree(tree.b, aes(color = group, linetype = group))  +
  scale_color_manual(values = c("black", "#33CCCC"), breaks = "CG") +
  scale_linetype_manual(values = c("solid", "dashed"), breaks = "CG") +
  labs(colour = "Insertion type", linetype = "Insertion type", 
       title = "(i) Congeneric insertion") +
  geom_tiplab(size = 3, show.legend = FALSE) + xlim(0, 370) +
  theme(legend.position = "none")


# second round - option 1: add into family node, as polytomy
tree.insert1<- bind.tip(tree.insert, "Rineloricaria_sp1", 
                        where = fastMRCA(tree.insert,"Loricaria_luciae", "Ancistrus_reisi"), 
                        position = 0)

tree.op3 <- groupOTU(tree.insert1, insertions.types)
tree.op3$tip.label <- gsub("_", " ", tree.op3$tip.label)

plot.op3 <- ggtree(tree.op3, aes(color = group, linetype = group))  +
  scale_color_manual(values = c("black", "#33CCCC", "black"), breaks = c("CG", "FN")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed"), breaks = c("CG", "FN")) +
  labs(colour = "Insertion type", linetype = "Insertion type", 
       title = "(ii) Family-level insertions", subtitle = "Option 3 - Family node") +
  geom_tiplab(size = 3, show.legend = F) + 
  xlim(0, 370) +
  theme(legend.position = "none")


# Second round - option 2: specifying a relationship with some genus
pos.genus <- which(c(tree.insert$tip.label, 
                     tree.insert$node.label) == "Loricaria_luciae") 

posit.genus <- tree.insert$edge.length[sapply(pos.genus, function(x, y) 
  which(y == x), y = tree.insert$edge[, 2])]

tree.insert2 <- bind.tip(tree.insert, "Rineloricaria_sp1",
                         where = pos.genus, position = posit.genus/2) #

tree.op1 <- groupOTU(tree.insert2, insertions.types)
tree.op1$tip.label <- gsub("_", " ", tree.op1$tip.label)


plot.op1 <- ggtree(tree.op1, aes(color = group, linetype = group))  +
  scale_color_manual(values = c("black", "#33CCCC", "black", "black"), 
                     breaks = c ("CG", "FN", "NG")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed", "dashed"), 
                        breaks = c ("CG", "FN", "NG")) +
  labs(colour = "Insertion type", linetype = "Insertion type", 
       subtitle = "Option 1 - Near to a genus", 
       title = "(ii) Family-level insertions") +
  geom_tiplab(size = 3, show.legend = F) + 
  xlim(0, 370) +
  theme(legend.position = "none")


# second round - option 3: specifying two genus
tree.insert3 <- bind.tip(tree.insert, "Rineloricaria_sp1", 
                         where = fastMRCA(tree.insert, "Curculionichthys_insperatus", "Loricaria_luciae"), 
                         position = 0)

tree.op2 <- groupOTU(tree.insert3, insertions.types)
tree.op2$tip.label <- gsub("_", " ", tree.op2$tip.label)

plot.op2 <- ggtree(tree.op2, aes(color = group, linetype = group))  +
  scale_color_manual(values = c("black", "#33CCCC", "black", "black","black"), 
                     breaks = c ("CG", "FN", "NG", "BG")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed", "dashed", "dashed"), 
                        breaks = c ("CG", "FN", "NG", "BG")) +
  labs(colour = "Insertion type", linetype = "Insertion type",
       title = "(ii) Family-level insertions",
       subtitle = "Option 2 - Between two genera") +
  geom_tiplab(size = 3, show.legend = FALSE) + 
  xlim(0, 370) +
  theme(legend.position = "none")

# third round: congeneric in the family round
tree.insert4 <- add.species.to.genus(tree = tree.insert3, species = "Rineloricaria_sp2")

tree.d <- groupOTU(tree.insert4, insertions.types)
tree.d$tip.label <- gsub("_", " ", tree.d$tip.label)

plot.d <- ggtree(tree.d, aes(color = group, linetype = group))  +
  scale_color_manual(values = c("black", "black", "black", "#33CCCC", "black", "black"), 
                     breaks = c ("CG", "FN", "NG", "BG", "CGF")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed", "dashed", "dashed", "dashed"), 
                        breaks = c ("CG", "FN", "NG", "BG", "CGF")) +
  labs(colour = "Insertion type", linetype = "Insertion type", 
       title = "(iii) Congeneric insertion at family-level") +
  geom_tiplab(size = 3, show.legend = F) + 
  xlim(0, 370) +
  theme(legend.position = "none")

# Finaly, we can cut the tree according with your species pool
spp.pool <- c("Cichlasoma_paranaense", "Ancistrus_reisi", "Ancistrus_pirareta", 
              "Rineloricaria_sp1", "Rineloricaria_sp2", "Hypostomus_iheringii",
              "Curculionichthys_insperatus", "Bryconamericus_exodon", "Astyanax_lacustris",
              "Gymnotus_inaequilabiatus", "Gymnotus_sp1")

spp.pool <- setNames(spp.pool, spp.pool)

tree.cutted <- geiger::treedata(phy = tree.insert4, spp.pool)$phy
tree.cutted <- groupOTU(tree.cutted, insertions.types)
tree.cutted$tip.label <- gsub("_", " ", tree.cutted$tip.label)

plot.e <- ggtree(tree.cutted, aes(linetype = group, color = group))  +
  labs(title = "Final tree", linetype = "Insertion type", color = "Insertion type") +
  scale_color_manual(values = c("black", viridis(5)), 
                     breaks = c("CG", "FN", "NG", "BG", "CGF")) +
  geom_treescale(x = 0, width = 20, fontsize = 3.5, offset = .15) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed", 
                                   "dashed", "dashed", "dashed"),
                        breaks = c("CG", "FN", "NG", "BG", "CGF")) +
  geom_tiplab(size = 3, show.legend = FALSE) +
  # guides(linetype = guide_legend(override.aes = list(size = .8))) +
  xlim(0, 370) + 
  theme(legend.position = "none")


# saving all figures together
squematic.fig <- cowplot::plot_grid(NULL, NULL,
                                    plot.op1,
                                    NULL, NULL,
                                    plot.a, plot.b,
                                    plot.op2, 
                                    plot.d, plot.e,
                                    NULL, NULL, 
                                    plot.op3, nrow = 3 , ncol = 5,
                                    scale = 0.9)

cowplot::save_plot(here::here("output", "images", "Squematic_fig.png"), 
                   squematic.fig, 
                   base_height = 15, 
                   base_width = 30, 
                   ncol = 1, nrow = 1)
