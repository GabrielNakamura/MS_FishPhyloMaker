# plotting sensitivty analysis


# reading data ------------------------------------------------------------

df_sensitivity <- readRDS(file = here::here("output", "df_sensitivity_cor.rds"))
library(tidyverse)


# plot --------------------------------------------------------------------


theme_set(theme_bw())
insertionPercentagelabx <- glue::glue("{seq(10, 25, by = 5)}")

distinctness <- ggplot(data = df_sensitivity, 
                       aes(x = percent_insert, y = cor_distinctness, fill = cor_distinctness)) +
  geom_boxplot() +
  geom_jitter(color="gray", size=0.8, alpha=0.6) +
  scale_x_discrete(labels= insertionPercentagelab) +
  scale_y_continuous(name = "Correlation coefficient",
                     limits = c(0, 1)) +
  ggtitle("A") +
  xlab("Percentage of inserted species") +
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
               aes(x = percent_insert, y = cor_cophenetic, fill = cor_cophenetic)) +
  geom_boxplot() +
  geom_jitter(color="gray", size=0.8, alpha=0.6) +
  scale_x_discrete(labels= insertionPercentagelab) +
  scale_y_continuous(name = "Correlation coefficient",
                     limits = c(0, 1)) +
  ggtitle("B") +
  xlab("Percentage of inserted species") +
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

library(patchwork)
sensitivity_plot <- distinctness|coph

ggsave(filename = here::here("output", "images", "sensitivity_plot.png"), plot = sensitivity_plot, 
       width = 10, height = 7,
       dpi = 500)
