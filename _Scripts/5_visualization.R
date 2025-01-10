###################################################
### Plots & Visualizations for Force Experiment ###
###################################################

# Author: Devan Pancura


### Import required packages ###

library(dplyr)
library(bayestestR)
library(ggplot2)


### Plot force accuracy by target force

force_acc <- force_means_pp %>%
  group_by(target) %>%
  summarize(
    acc_pct = mean(target_achieved) * 100
  )

ggplot(force_acc, aes(x = as.factor(target), y = acc_pct)) +
  geom_col(color = "lightskyblue3", fill = "lightskyblue3") +
  labs(x = "Force (% MVC)", y = "Target Achieved (%)") +
  geom_text(
    aes(label = paste0(sprintf(acc_pct, fmt = "%#.1f"), "%")),
    size = 4, nudge_y = -5, colour = "white"
  ) +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 5))
  ) +
  theme_classic()



### Plot the patterns of results from the models ###

pd <- position_dodge(0.55)

median_eti <- function(emm, .width = 0.9) {
  # Computes ETIs and merges then on to the EMM table
  emm_eti <- as_tibble(eti(emm, ci = .width)) %>%
    separate("Parameter", c("target", "order"), extra = "merge")
  left_join(as_tibble(emm), emm_eti, by = c("target", "order"))
}


# MEP amplitudes for motor execution trials

ggplot(median_eti(emm_phys), aes(x = target, y = emmean, color = order)) +
  geom_point(position = pd, size = 2.2) +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high), width = 0,
    linewidth = 0.8, position = pd
  ) +
  xlab("Target Force (% MVC)") +
  ylab("MEP Amplitude (Standardized)") +
  labs(color = "Order") +
  scale_color_manual(values = c("lightskyblue3", "lightskyblue4")) +
  scale_y_continuous(breaks = c(-1.5, -1.0, -0.5, 0, 0.5, 1.0)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 5))
  )


# MEP ampitudes for motor imagery trials

ggplot(median_eti(emm_mi), aes(x = target, y = emmean, color = order)) +
  geom_point(position = pd, size = 2.2) +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high), width = 0,
    linewidth = 0.8, position = pd
  ) +
  xlab("Target Force (% MVC)") +
  ylab("MEP Amplitude (Standardized)") +
  labs(color = "Order") +
  scale_color_manual(values = c("lightskyblue3", "lightskyblue4")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 5))
  )
