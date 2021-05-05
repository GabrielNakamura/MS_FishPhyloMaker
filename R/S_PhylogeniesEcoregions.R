# loading packages -----------------------------------------------
library(ggtree)
library(ggplot2)
library(cowplot)


# number of insertions ----------------------------------------------------
percentage_insertion_afrotropics <- round(100*(length(unlist(lapply(c("Congeneric_insertion", "Congeneric_insertion_roundFamily", "Family_insertion", "Order_insertion"), function(x){
  which(x == res_afrotropics$Insertions_data$insertions)
})))/dim(res_afrotropics$Insertions_data)[1]), 1)

percentage_insertion_neotropic <- round(100*(length(unlist(lapply(c("Congeneric_insertion", "Congeneric_insertion_roundFamily", "Family_insertion", "Order_insertion"), function(x){
  which(x == res_neotropic$Insertions_data$insertions)
})))/dim(res_neotropic$Insertions_data)[1]), 1)

percentage_insertion_indomalay <- round(100*(length(unlist(lapply(c("Congeneric_insertion", "Congeneric_insertion_roundFamily", "Family_insertion", "Order_insertion"), function(x){
  which(x == res_indomalay$Insertions_data$insertions)
})))/dim(res_indomalay$Insertions_data)[1]), 1)

percentage_insertion_neartic <- round(100*(length(unlist(lapply(c("Congeneric_insertion", "Congeneric_insertion_roundFamily", "Family_insertion", "Order_insertion"), function(x){
  which(x == res_neartic$Insertions_data$insertions)
})))/dim(res_neartic$Insertions_data)[1]), 1)

# plotting trees and insertions -------------------------------------------
cols <- c(Present_in_Tree = 'lightblue', 
          Congeneric_insertion = 'yellow2', 
          Congeneric_insertion_roundFamily = 'orange', 
          Family_insertion = "red",
          Order_insertion = "purple")

p.base_afrotropics <- ggtree(res_afrotropics$Phylogeny, layout = "circular", size = .3)  %<+% res_afrotropics$Insertions_data +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) + #  plot the scale bar
  annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar
phylo_afrotropics <- p.base_afrotropics + 
  geom_tippoint(aes(color = insertions), size = .5, alpha = .8) +
  scale_color_manual(values = cols, name = "Insertions", labels = c("Congeneric", 
                                                                    "Congeneric in Family",
                                                                    "Family", 
                                                                    "Order",
                                                                    "Present in tree")
  )
p.base_indomalay <- ggtree(res_indomalay$Phylogeny, layout = "circular", size = .3)  %<+% res_indomalay$Insertions_data +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) + #  plot the scale bar
  annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar
phylo_indomalay <- p.base_indomalay + 
  geom_tippoint(aes(color = insertions), size = .5, alpha = .8) +
  scale_color_manual(values = cols, name = "Insertions", labels = c("Congeneric", 
                                                                    "Congeneric in Family",
                                                                    "Family", 
                                                                    "Order",
                                                                    "Present in tree")
  )
p.base_neartic <- ggtree(res_neartic$Phylogeny, layout = "circular", size = .3)  %<+% res_neartic$Insertions_data +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) + #  plot the scale bar
  annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar
phylo_neartic <- p.base_neartic + 
  geom_tippoint(aes(color = insertions), size = .5, alpha = .8) +
  scale_color_manual(values = cols, name = "Insertions", labels = c("Congeneric", 
                                                                    "Congeneric in Family",
                                                                    "Family", 
                                                                    "Order",
                                                                    "Present in tree")
                     )
p.base_neotropic <- ggtree(res_neotropic$Phylogeny, layout = "circular", size = .3)  %<+% res_neotropic$Insertions_data +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) + #  plot the scale bar
  annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar
phylo_neotropic <- p.base_neotropic + 
  geom_tippoint(aes(color = insertions), size = .5, alpha = .8) +
  scale_color_manual(values = cols, name = "Insertions", labels = c("Congeneric", 
                                                                    "Congeneric in Family",
                                                                    "Family", 
                                                                    "Order",
                                                                    "Present in tree"))

pgrid <- plot_grid(
  phylo_afrotropics + theme(legend.position="none"),
  phylo_indomalay + theme(legend.position="none"),
  phylo_neartic + theme(legend.position="none"),
  phylo_neotropic + theme(legend.position="none"), 
  labels = c(paste('Afrotropic', paste(percentage_insertion_afrotropics, "%", sep = ""), sep = " "),
             paste('Indo-Malay', paste(percentage_insertion_indomalay, "%", sep = ""), sep = " "),
             paste("Neartic", paste(percentage_insertion_neartic, "%", sep = ""), sep = " "),
             paste("Neotropic", paste(percentage_insertion_neotropic, "%", sep = ""), sep = " ")),
  nrow = 2
)

# extract the legend from one of the plots
legend <- get_legend(
  phylo_afrotropics + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 2)))
)

phylo_all <- plot_grid(pgrid, 
                       legend, ncol = 1, rel_heights = c(1, .1)
)


# saving phylo image ------------------------------------------------------
data_saving <- date()
data_edit <- unlist(strsplit(data_saving, split = " "))[c(1, 2, 3, 5)]
name_plot <- paste("Fig_phyloEcoregions", 
                   data_edit[1],
                   data_edit[2], 
                   data_edit[3],
                   data_edit[4],
                   sep = "_")
name_plot_final <- paste(name_plot, "pdf", sep = ".")
pdf(file = here::here("output", "images", name_plot_final), width = 10, height = 10)
print(phylo_all)
dev.off()
