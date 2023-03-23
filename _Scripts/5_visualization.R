###################################################
### Plots & Visualizations for Force Experiment ###
###################################################

# Author: Devan Pancura


### Import required packages ###

library(ggplot2)
library(ez)


### Define plot themes ###

plot_theme <- theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    legend.key = element_blank(),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 22),
    axis.text = element_text(color = "black", size = 18),
    axis.title.y = element_text(margin = margin(r = 0.5, unit = "cm")),
    axis.title.x = element_text(margin = margin(t = 0.5, unit = "cm")),
    legend.title = element_text(size = 22),
    axis.ticks = element_line(color = "black", linewidth = 1),
    legend.text = element_text(size = 22),
    legend.spacing.x = unit(0.5, "cm"),
    legend.box.margin = margin(l = 1, unit = "cm"),
  )


# Plot some descriptives about force trial data

ggplot(force_data, aes(x = target)) +
  geom_bar(color = "navyblue", fill = "navyblue", width = 3) +
  labs(x = "Force (% MVC)", y = "Number of Successful Trials") +
  geom_text(
    aes(label = after_stat(count)), stat = "count",
    size = 8, vjust = 1.5, colour = "white"
  ) +
  plot_theme

ggplot(accuracy, aes(x = target, y = accuracy)) +
  geom_point(color = "navyblue", size = 2) +
  labs(x = "Force (% MVC)", y = "Accuracy (%)") +
  geom_errorbar(
    aes(ymin = accuracy - sd, ymax = accuracy + sd),
    linewidth = 0.3, width = 1, col = "navyblue"
  ) +
  plot_theme

ggplot(avg_accuracy, aes(x = target, y = avg_acc)) +
  geom_point(size = 4.5, col = "navyblue") +
  geom_errorbar(
    aes(ymin = avg_acc - sd, ymax = avg_acc + sd),
    linewidth = 0.5, col = "navyblue"
  ) +
  xlab("Force (% MVC)") +
  ylab("Accuracy (%)") +
  plot_theme
