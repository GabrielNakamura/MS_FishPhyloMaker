# plotting sensitivity analysis


# reading data ------------------------------------------------------------

df_sensitivity <- readRDS(file = here::here("output", "df_sensitivity_cor.rds"))
library(tidyverse)
library(cowplot)
library(ggforce)

# plot --------------------------------------------------------------------


theme_set(theme_bw())
insertionPercentagelab <- glue::glue("{seq(10, 60, by = 5)}")

distinctness <- ggplot(data = df_sensitivity, 
                       aes(x = percent_insert, y = cor_distinctness, fill = cor_distinctness)) +
  geom_boxplot() +
  geom_jitter(color="gray", size=0.6, alpha=0.6) +
  scale_x_discrete(labels= insertionPercentagelab) +
  scale_y_continuous(limits = c(.4, 1)) +
  labs(x = "Percentage of inserted species (%)", 
       y = "Correlation coefficient", 
       tag = "A") +
  theme(legend.position = "", panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "mm"),
        legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
        legend.text = element_text(family = "Times", color = "black", size = 12), 
        axis.text = element_text(family = "Times", color = "black", size = 12), panel.border = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 9, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6)
        )
        )


coph <- ggplot(data = df_sensitivity, 
               aes(x = percent_insert, y = cor_cophenetic)) +
  geom_boxplot() +
  geom_jitter(color="gray", size=0.8, alpha=0.6) +
  scale_x_discrete(labels= insertionPercentagelab) +
  labs(x = "Percentage of inserted species (%)",
       y = "Correlation coefficient",
       tag = "B") +
  theme(legend.position = "", panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "mm"),
        legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
        legend.text = element_text(family = "Times", color = "black", size = 12), 
        axis.text = element_text(family = "Times", color = "black", size = 12), 
        panel.border = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(family = "Arial", 
                                     color = "black",
                                     size = 9, 
                                     hjust = 0.5, 
                                     margin = margin(b = 6)
        )
  )

sensitivity_plot <- cowplot::plot_grid(distinctness, coph, nrow = 2)



# saving graphic ----------------------------------------------------------

ggsave(filename = here::here("output", "images", "sensitivity_plot.png"), plot = sensitivity_plot, 
       width = 7, height = 9,
       dpi = 500)
